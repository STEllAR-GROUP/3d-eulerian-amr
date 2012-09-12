#include "../defs.h"
#ifdef BINARY
#include "binary_driver.h"
#include <stdlib.h>
#include <omp.h>
#include <time.h>

#include "initialize.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "binary.h"
#include "../state.h"
#include <omp.h>

void BinaryDriver::set_output_frequency_by_time(Real t) {
	int_output_freq = -1;
	time_output_freq = t;
}

void BinaryDriver::reset_time() {
	this->ostep_cnt = 0;
	root->set_time(0.0);
	step_cnt = 0;
}

void BinaryDriver::set_verbosity(int a) {
	assert( a >= 0 );
	verbosity = a;
}

void BinaryDriver::set_output_frequency_by_step(int a) {
	assert( a >= 1 );
	int_output_freq = a;
}

void BinaryDriver::step_to_time(Real tmax) {
	root->set_time(0.0);
	clock_t start_time = clock();
	bool last_step = false;
	bool do_output;
	Real tleft, dt, last_dt;
	Real next_output_time;
	dt = 0.0;
	Vector<State, 4> state_sum;
	if (get_time() == 0.0 && int_output_freq != 0) {
#ifndef READ_FROM_FILE
		int iters = 0;
		iters = 0;
		root->exec_function(initialize);
		root->enforce_boundaries();
		root->inject_from_children();
		dt = root->max_dt();
		dt *= min(1.0, MAXINITDT);
		int last_max_level;
		int last_max_level2;
		do {
			last_max_level2 = last_max_level;
			last_max_level = root->max_level();
			this->step(dt, true);
			//			root->check_for_refine();
			if (verbosity >= 1) {
				state_sum = root->state_sum();
				printf("%4i %10.3e %10.3e %4i %4i ", step_cnt, dt, get_time(), root->ngrids(), root->max_level());
				/*		printf("%10.4e %10.4e %10.4e \n", state_sum[1].rho() / state_sum[0].rho(), state_sum[2].rho() / state_sum[0].rho(),
				 state_sum[3].rho() / state_sum[0].rho());*/
				printf("\n");
			}
			root->set_time(0.0);
			root->exec_function(initialize);
			root->enforce_boundaries();
			root->inject_from_children();
			dt = min(root->max_dt(), MAXINITDT);
			step_cnt = 0;
			iters++;
		}while (root->max_level() < GridNode::max_refine_level || iters <= 2);
#else
		FILE* fp = fopen(READ_FROM_FILE, "rb");
		assert(fp);
		root->read(fp);
		fclose(fp);
#endif
		printf("%e\n", Binary::Omega);
#ifndef SCF_CODE
		this->set_output_frequency_by_time(
				2.0 * M_PI / Binary::Omega0 / Real(FRAME_RATE));
#endif
		ostep_cnt = int(root->get_time() / time_output_freq);
		root->pot_from_grid();
		solve();
		root->GridNode::output("X", nint(get_time() / time_output_freq));
	}
	tmax *= 2.0 * M_PI / Binary::Omega0;
	tmax += get_time();
#ifdef SCF_CODE
	tmax = 1.0 + 100;
#endif
	GridNode::set_time_normal(2.0 * M_PI / Binary::Omega0);
	do {
		Real start_time = omp_get_wtime();
		dt = root->max_dt();
		if (get_time() == 0.0) {
			dt = MAXINITDT * dt;
		} else if (step_cnt != 0) {
			dt = min(last_dt * MAXDTINC, dt);
		}
		last_dt = dt;
		do_output = false;
		tleft = tmax - get_time();
		if (int_output_freq < 0) {
			next_output_time = Real(ostep_cnt + 1) * time_output_freq;
			if (dt + get_time() >= next_output_time) {
				dt = next_output_time - get_time();
				ostep_cnt++;
				do_output = true;
			} else {
				dt = min(
						dt,
						(next_output_time - get_time()) / int(
								(next_output_time - get_time()) / dt + 1.0));
			}
		} else if (int_output_freq != 0) {
			do_output = bool(step_cnt % int_output_freq == 0);
		}
		//		printf( "%e\n", dt );
		if (step_cnt % CHKPT_FREQ == 0) {
			FILE * fp;
			if (step_cnt % (2 * CHKPT_FREQ) == 0) {
				fp = fopen("X.hello.chk", "wb");
			} else {
				fp = fopen("X.goodbye.chk", "wb");
			}
			root->write(fp);
			fclose(fp);
		}
		if (tleft <= dt) {
			dt = tleft;
			last_step = true;
			if (int_output_freq < 0 && next_output_time == tmax) {
				do_output = true;
			}
		}
		if (verbosity >= 50) {
			root->output_tree();
		}
		if (verbosity >= 1) {
			printf("\n");
		}
		dt = this->step(dt);
		if (verbosity >= 1) {
			state_sum = root->state_sum() + root->flow_off();
			printf("%4i %.3e %.3e % 5i %5i %5i %6.2f", step_cnt, get_time(),
					dt, root->ngrids(), root->max_level(),
					root->fine_point_count(), omp_get_wtime() - start_time);
			state_sums_out(stdout, state_sum[0]);
			FILE* fp = fopen("ts.dat", "at");
			fprintf(fp, "%e %e\n", get_time(), dt);
			fclose(fp);
		}
		if (do_output) {
			root->GridNode::output("X", nint(get_time() / time_output_freq));
			printf("\n --- output ---");
		}
		_3Vec com;
		Vector<Real, 4> m1;
		Vector<Real, 4> m2;
		m1 = root->mass_sum(0);
		m2 = root->mass_sum(1);
		for (int i = 0; i < 3; i++) {
			com[i] = (m1[0] * m1[i + 1] + m2[0] * m2[i + 1]) / (m1[0] + m2[0]);
		}
		FILE* f = fopen("com.dat", "at");
		fprintf(f, "%e %e %e\n", root->get_time(), com[0], com[1]);
		fclose(f);

#ifdef CENTER_OF_MASS_CORRECTION
		Real dR, dRdot;
		const Real w0 = 4.0 * Binary::Omega;
		dR = sqrt(com[0] * com[0] + com[1] * com[1]);
		dRdot = (dR - sqrt(
				last_com[0] * last_com[0] + last_com[1] * last_com[1])) / dt;
		if (last_dx != 0.0) {
			State::fR = -w0 * (w0 * dR + 2.0 * dRdot);
			State::ftheta = atan2(com[1], com[0]);
			FILE* f = fopen("com_corr.dat", "at");
			fprintf(f, "%0.8e %0.8e %0.8e %0.8e %0.8e\n", get_time(), dR,
					dRdot, State::fR, State::ftheta);
			fclose(f);
		}
#endif
#ifdef DYNAMIC_OMEGA
		Real domega = (m1[2] - m2[2]) / (m1[1] - m2[1]);
		if (last_dx != 0.0) {
			domega -= (last_dy / last_dx);
			domega /= dt;
			domega += Binary::Omega * (m1[2] - m2[2]) / (m1[1] - m2[1]);
			Binary::Omega += domega;
			FILE* f = fopen("domega.dat", "at");
			Binary::phase += (Binary::Omega - Binary::Omega0) * dt;
			fprintf(f, "%0.8e %0.8e %0.8e %0.8e %0.8e\n", get_time(),
					Binary::Omega, domega, Binary::Omega0, Binary::phase);
			fclose(f);
		}
		Binary::phase += (Binary::Omega - Binary::Omega0) * dt;
#endif
		last_com = com;
		last_dy = m1[2] - m2[2];
		last_dx = m1[1] - m2[1];

	} while (!last_step);
	printf("time = %e\n", (clock() - start_time) / Real(CLOCKS_PER_SEC));
}

void BinaryDriver::set_output_off() {
	int_output_freq = 0;
}

Real BinaryDriver::get_time() const {
	return root->get_time();
}

BinaryDriver::BinaryDriver(Binary* _root) {
	ChildIndex i;
	root = _root;
	step_cnt = 0;
	ostep_cnt = 0;
	set_output_off();
#ifndef NDEBUG
	set_verbosity(DEBUG_V);
#else
	set_verbosity(DEBUG_V);
#endif
	hydrotime = poissontime = 1.0;
	last_dy = 0.0;
	last_dx = 0.0;
}

void BinaryDriver::state_sums_out(FILE* fp, const State& s) {
	//	Real half_phi = root->pot_sum();
	//	for (int i = 0; i < STATE_NF; i++) {
	//	fprintf(fp, "%12.5e ", s[i]);
	//	}
	/*hydrotime = omp_get_wtime() - poissontime - hydrotime;
	 Real ttime = (poissontime + hydrotime) * 0.01;
	 fprintf(fp, "%5.2f ", poissontime + hydrotime);
	 poissontime = 0.0;
	 hydrotime = omp_get_wtime();*/
	FILE* f;
	State s0;
	f = fopen("sums.dat", "at");
	fprintf(f, "%e ", root->get_time());
	for (int i = 0; i < STATE_NF; i++) {
		fprintf(f, "%.14e ", s[i]);
	}
	fprintf(f, "\n");
	fclose(f);
	s0 = root->flow_off();
	f = fopen("flow.dat", "at");
	fprintf(f, "%e ", root->get_time());
	for (int i = 0; i < STATE_NF; i++) {
		fprintf(f, "%e ", s0[i]);
	}
	fprintf(f, "\n");
	fclose(f);
	s0 = root->state_max();
	State::rho_floor = max(s0[0], s0[1]) / 1.0e+14;
	f = fopen("max.dat", "at");
	fprintf(f, "%e ", root->get_time());
	for (int i = 0; i < STATE_NF; i++) {
		fprintf(f, "%e ", s0[i]);
	}
	fprintf(f, "\n ");
	fclose(f);

	s0 = root->state_min();
	f = fopen("min.dat", "at");
	fprintf(f, "%e ", root->get_time());
	for (int i = 0; i < STATE_NF; i++) {
		fprintf(f, "%e ", s0[i]);
	}
	fprintf(f, "\n");
	fclose(f);

}

void BinaryDriver::solve(Real ell_toler) {
	Real start_time = omp_get_wtime();
	Real error, dlzerror;
	int iter;
	iter = 0;
	if (DEBUG_V >= 40) {
		printf("Computing Poisson boundary conditions\n");
	}
	root->compute_source_sums();
	root->compute_physical_boundaries(root);
	if (DEBUG_V >= 40) {
		printf("Solving Poisson Equation\n");
	}
	root->compute_forces();
	//	error = sqrt(root->dphil2());
	Real last_dlzerror;
	dlzerror = 1.0e+99;
	do {
		if (iter > MAX_POISSON_ITER) {
			printf("Max Poisson iterations\n");
			//	root->output(nint(get_time() / time_output_freq));
			abort();
		}
		error = root->vcycle();
		last_dlzerror = dlzerror;
		dlzerror = root->sum_dlz();
		iter++;
		FILE* fp = fopen("poisson.dat", "at");
		fprintf(fp, "%5i %e %e %e\n", iter, error, dlzerror,
				fabs(log(fabs(dlzerror / last_dlzerror))));
		fclose(fp);
#ifdef SCF_CODE
		if (error < ELL_TOLER * 1000.0 && iter >= 5) {
			break;
		}
	}while (error > ell_toler || (iter < 5));
#else
	} while (error > ell_toler);
#endif
	root->pot_to_grid();
	root->inject_from_children();
	poissontime += omp_get_wtime() - start_time;

}

void BinaryDriver::sub_step(Real dt, Real beta, int rk) {
	FILE* fp = fopen("poisson.dat", "at");
	fprintf(fp, "%i\n", rk);
	fclose(fp);
	if (DEBUG_V >= 40) {
		printf("Doing Hydro Step\n");
	}
	root->enforce_boundaries();
	root->clear_difs();
	root->compute_x_flux();
	root->adjust_x_flux();
	root->sum_x_difs();
	root->compute_y_flux();
	root->adjust_y_flux();
	root->sum_y_difs();
	root->compute_z_flux();
	root->adjust_z_flux();
	root->sum_z_difs();

	root->error_from_parent();
	root->add_difs(dt, beta);
#ifndef SCF_CODE
	if (rk == 1) {
		root->store_phi();
		solve();
	} else if (rk == 2) {
		root->restore_half_phi();
		solve();
	} else if (rk == 3) {
		root->restore_phi();
		solve();
	}
	//	printf("%e\n", root->domega());
#else
	solve();
#endif
}

Real BinaryDriver::step(Real dt, bool test_step) {
	assert( dt > 0.0 );
	bool resolve = false;
	if (step_cnt <= 0) {
#ifdef READ_FROM_FILE
		root->pot_from_grid();
#endif
		solve();
		root->store_phi();
		FILE* fp = fopen("X.start.chk", "wb");
		root->write(fp);
		fclose(fp);
	}
	Binary::binary_integrals_t B;
#ifndef SCF_CODE
	root->M1M2data(&B);
	FILE* fp = fopen("lobe.dat", "at");
	fprintf(
			fp,
			"%e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e %0.8e\n",
			get_time()/*1*/, B.q/*2*/, B.a/*3*/, B.jorb/*4*/,
			B.total_energy/*5*/, B.kinetic/*6*/, B.m1/*7*/, B.js1/*8*/,
			B.I1/*9*/, B.m2/*10*/, B.js2/*11*/, B.I2/*12*/, B.mc/*13*/,
			B.jc/*14*/, B.Ic/*15*/, B.V1/*16*/, B.V2/*17*/);
	fclose(fp);
#endif
	root->store();
#ifdef SCF_CODE
	if (!test_step) {
		this->sub_step(dt, 1.0, 1);
	}
#else
	REDO: this->sub_step(dt, 1.0, 1);
	this->sub_step(dt, 0.25, 2);
	this->sub_step(dt, 2.0 / 3.0, 3);
	Real test_dt = root->max_dt();
	if (test_dt < dt / MAXDTINC) {
		printf("Timestep of %e is too big, trying %e\n", dt, dt / MAXDTINC);
		dt = dt / MAXDTINC;
		root->restore();
		goto REDO;
	}
#endif
	root->enforce_boundaries();
	if (step_cnt % 8 == 0 || test_step) {
		State smax = root->state_max();
		Real mind = max(smax.rho(0), smax.rho(1));
		mind /= DYNAMIC_RANGE;
		//	printf("\nSetting minimum refinement density to %e\n", mind);
		GridNode::min_refine_density = mind;
#ifndef FIXED_MESH
#ifndef SCF_CODE
		root->add_pot_et();
#endif
		bool clear_flags;
		clear_flags = (step_cnt % 500) == 0;
		if (root->check_for_refine(clear_flags)) {
			if (!test_step) {
				solve(ELL_TOLER / 10.0);
				root->store_phi();
				resolve = true;
				Vector<Real, 4> m1 = root->mass_sum(0);
				Vector<Real, 4> m2 = root->mass_sum(1);
				Vector<Real, 4> l1 = root->find_l1(m1, m2);
				root->mark_lobes(l1[0], l1[1], l1[2], m1[1], m2[1], m1[2],
						m2[2]);
			}
		}
#ifndef SCF_CODE
		root->sub_pot_et();
#endif
#endif

	}
	step_cnt++;
	root->set_time(get_time() + dt);
	static Real omega_p1, K2, K1, virial;
#ifdef SCF_CODE
	if (get_time() > 10.0 ) {
		if (fabs(log(K1 / K2)) < 1.0e-5) {
			this->reset_time();
			Binary::phase = 0.0;
			Binary::Omega0 = Binary::Omega;
			this->set_output_frequency_by_time(2.0 * M_PI / Binary::Omega0 / Real(FRAME_RATE));
			char buffer[256];
			FILE* fp = fopen("X.start.chk", "wb");
			root->write(fp);
			fclose(fp);
			sprintf(buffer, "cp X.start.chk X.%i.%.2lf.%.2lf.%.2lf.chk\n", GridNode::max_refine_level, Binary::q,
					Binary::fill_factor, Binary::a);
			system(buffer);
			system("mkdir scf_data\n");
			system("mv *.chk ./scf_data\n");
			system("mv *.gz ./scf_data\n");
			system("mv *.dat ./scf_data\n");
			//			printf("\nSCF complete - Starting hydro-code\n");
			printf("\nSCF complete\n");
			exit(1);
		}
	}
	if (get_time() > 0.0 && !test_step) {
		Real phi1max, phi2min;
		Vector<Real, 4> m1 = root->mass_sum(0);
		Vector<Real, 4> m2 = root->mass_sum(1);
		Real phi0_1_np1;
		Vector<Real, 4> l1;
		l1 = root->find_l1(m1, m2);
		omega_p1 = root->next_omega();
		Binary::l1x = l1[1];
		root->mark_lobes(l1[0], l1[1], 0.0,m1[1], m2[1],0.0,0.0);
		phi2min = root->find_phimin(1.0e-6, 2);
		Binary::phi0_2 = l1[0] + (phi2min - l1[0]) * (1.0 - Binary::fill_factor);
		if (Binary::phi0_1 > l1[0]) {
			phi1max = root->find_phimax(1.0e-6, 1);
			Binary::phi0_1 = phi1max;
			printf("%e\n", phi1max);
		}
		K1 = root->find_K(Binary::phi0_1, 1);
		K2 = root->find_K(Binary::phi0_2, 2);
		Real eps = 1.0 + 1.0e-1;
		Binary::phi0_1 *= pow(K1 / K2, dt);
		Binary::Omega *= pow(omega_p1 / Binary::Omega, dt);
		Binary::K1 = K1;
		Binary::K2 = K2;
		Real d0 = 1.01;
		Real xcom = (m1[1] * m1[0] + m2[1] * m2[0]) / (m1[0] + m2[0]);
		Real ycom = (m1[2] * m1[0] + m2[2] * m2[0]) / (m1[0] + m2[0]);
		_3Vec origin = root->get_origin();
		origin[0] += 1.5 * xcom;
		origin[1] += 1.5 * ycom;
		FILE* fp = fopen("origin.dat", "at");
		fprintf(fp, "%e %e %e %e %e %e\n", root->get_time(), origin[0], origin[1], origin[2], xcom, ycom);
		fclose(fp);
		root->set_origin(origin);
		printf("t=%e  m1=%e m2=%e Omega=%e xcom=%e 	Kerr=%e	 ", root->get_time(), m1[0], m2[0], Binary::Omega,
				xcom, fabs(log(K1 / K2)));
	}
#endif
	FILE* f = fopen("lz.dat", "at");
	fprintf(f, "%e %e %e\n", root->get_time(), root->sum_lz(), root->sum_dlz());
	fclose(f);
	return dt;
}

BinaryDriver::~BinaryDriver() {
}
#endif
