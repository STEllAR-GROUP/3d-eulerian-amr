#include "../defs.h"
#ifdef SINGLE
#define USE_HYDRODRIVER
#endif
#ifdef EULER
#define USE_HYDRODRIVER
#endif

#ifdef USE_HYDRODRIVER
#include "hydro_driver.h"
#include "initialize.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

HydroDriver::HydroDriver(GridNode* a) {
	ChildIndex i;
	root = a;
	step_cnt = 0;
	ostep_cnt = 0;
	set_output_off();
	set_rk(RK_ORDER);
#ifndef NDEBUG
	set_verbosity(DEBUG_V);
#else
	set_verbosity(DEBUG_V);
#endif
}

void HydroDriver::set_rk(int a) {
	assert( a >= 1 );assert( a <= 3 );
	rk = a;
}

void HydroDriver::set_output_frequency_by_time(Real t) {
	int_output_freq = -1;
	time_output_freq = t;
}

void HydroDriver::reset_time() {
	this->ostep_cnt = 0;
	root->set_time(0.0);
	step_cnt = 0;
}

void HydroDriver::set_verbosity(int a) {
	assert( a >= 0 );
	verbosity = a;
}

void HydroDriver::set_output_frequency_by_step(int a) {
	assert( a >= 1 );
	int_output_freq = a;
}

int HydroDriver::get_rk() const {
	return rk;
}

void HydroDriver::state_sums_out(FILE* fp, const State& s) {
	for (int i = 0; i < STATE_NF; i++) {
		fprintf(fp, "%10.3e ", s[i]);
	}
}

HydroDriver::~HydroDriver() {
}

void HydroDriver::step_to_time(Real tmax) {
	root->set_time(0.0);
	clock_t start_time = clock();
	bool last_step = false;
	bool do_output;
	Real tleft, dt, last_dt;
	Real next_output_time;
	dt = 0.0;
	Vector<State, 4> state_sum;
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
		if (verbosity >= 1) {
			state_sum = root->state_sum();
			printf("%4i %10.3e %10.3e %4i %4i ", step_cnt, dt, get_time(), root->ngrids(), root->max_level());
			printf("\n");
		}
		root->set_time(0.0);
		root->exec_function(initialize);
		root->enforce_boundaries();
		root->inject_from_children();
		dt = min(root->max_dt(), MAXINITDT);
		step_cnt = 0;
		iters++;
	} while (root->max_level() < GridNode::max_refine_level || iters <= 2);
	root->output("X", 0.0);
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
				dt = min(dt, (next_output_time - get_time()) / int((next_output_time - get_time()) / dt + 1.0));
			}
		} else if (int_output_freq != 0) {
			do_output = bool(step_cnt % int_output_freq == 0);
		}
		//		printf( "%e\n", dt );
		if (step_cnt % CHKPT_FREQ == 0) {
			FILE* fp;
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
		this->step(dt);
		if (verbosity >= 1) {
			state_sum = root->state_sum() + root->flow_off();
			printf("%4i %.3e %.3e % 5i %5i %5i %6.2f", step_cnt, get_time(), dt, root->ngrids(), root->max_level(),
					root->fine_point_count(), omp_get_wtime() - start_time);
			state_sums_out(stdout, state_sum[0]);
			FILE* fp = fopen("ts.dat", "at");
			fprintf(fp, "%e %e\n", get_time(), dt);
			fclose(fp);
		}
		if (do_output) {
			root->output("X", nint(get_time() / time_output_freq));
			printf("\n --- output ---");
		}
	} while (!last_step);
	printf("time = %e\n", (clock() - start_time) / Real(CLOCKS_PER_SEC));
}

void HydroDriver::step(Real dt, bool) {
	assert( dt > 0.0 );
	root->store();
	if (rk == 1) {
		this->sub_step(dt, 1.0, 1);
		root->inject_from_children();
	} else if (rk == 2) {
		this->sub_step(dt, 1.0, 1);
		root->inject_from_children();
		this->sub_step(dt, 0.5, 2);
		root->inject_from_children();
	} else {
		this->sub_step(dt, 1.0, 1);
		root->inject_from_children();
		this->sub_step(dt, 0.25, 2);
		root->inject_from_children();
		this->sub_step(dt, 2.0 / 3.0, 3);
		root->inject_from_children();
	}
	root->enforce_boundaries();
	root->check_for_refine();
	step_cnt++;
	root->set_time(get_time() + dt);
}

void HydroDriver::step() {
	step(root->max_dt());
}

void HydroDriver::sub_step(Real dt, Real beta, int) {
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
}

Real HydroDriver::get_time() const {
	return root->get_time();
}

void HydroDriver::set_output_off() {
	int_output_freq = 0;
}

#endif
