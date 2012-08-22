#ifndef SINGLE_DRIVER_H_
#define SINGLE_DRIVER_H_
#ifdef SINGLE

#include "single.h"
#include "euler/hydro_driver.h"

class SingleDriver : public HydroDriver {
private:
	Single* root;
public:
	SingleDriver(Single* s) : HydroDriver(s) {
		root=s;
	}
	virtual ~SingleDriver() {

	}
	virtual void state_sums_out(const State& s);
	void solve(Real ell_toler) {
		int iter;
		root->compute_source_sums();
		root->compute_physical_boundaries(root);
		root->compute_forces();
		Real error;
		iter = 0;
		do {
			if (iter > MAX_POISSON_ITER) {
				printf("Max Poisson iterations\n");
				//	root->output(nint(get_time() / time_output_freq));
				abort();
			}
			error = root->vcycle();
			iter++;
			FILE* fp = fopen("poisson.dat", "at");
			fprintf(fp, "%5i %e\n", iter, error );
			fclose(fp);
		}
		while (error > ell_toler);
		root->pot_to_grid();
		root->inject_from_children();
	}
	virtual void sub_step(Real dt, Real beta, int rk) {
		FILE* fp = fopen("poisson.dat", "at");
		fprintf(fp, "%i\n", rk);
		fclose(fp);
		if (DEBUG_V >= 40) {
			printf("Doing Hydro Step\n");
		}
		root->enforce_boundaries();
		root->clear_difs();
		root->compute_x_flux(false);
		root->adjust_x_flux();

		root->sum_x_difs(false);
		root->compute_y_flux(false);
		root->adjust_y_flux();
		root->sum_y_difs(false);
		root->compute_z_flux(false);
		root->adjust_z_flux();
		root->sum_z_difs(false);

		root->compute_x_flux(true);
		root->adjust_x_flux();
		root->sum_x_difs(true);
		root->compute_y_flux(true);
		root->adjust_y_flux();
		root->sum_y_difs(true);

		root->error_from_parent();
		root->add_difs(dt, beta);
		if (rk == 1) {
			root->store_phi();
			solve(ELL_TOLER);
		} else if (rk == 2) {
			root->restore_half_phi();
			solve(ELL_TOLER);
		} else if (rk == 3) {
			root->restore_phi();
			solve(ELL_TOLER);
		}
		//	printf("%e\n", root->domega());
	}
	virtual void step(Real dt, bool test_step) {
		assert( dt > 0.0 );
		if (step_cnt <= 0) {
#ifdef READ_FROM_FILE
			root->pot_from_grid();
#endif
			solve(ELL_TOLER);
			root->store_phi();
			FILE* fp = fopen("X.start.chk", "wb");
			root->write(fp);
			fclose(fp);
		}
		root->store();
		this->sub_step(dt, 1.0, 1);
		this->sub_step(dt, 0.25, 2);
		this->sub_step(dt, 2.0 / 3.0, 3);
		Real test_dt = root->max_dt();
		root->enforce_boundaries();
		if (step_cnt % 8 == 0 || test_step) {
			State smax = root->state_max();
			Real mind = smax.get_rho(0);
			mind /= DYNAMIC_RANGE;
			GridNode::min_refine_density = mind;
			root->add_pot_et();
			bool clear_flags;
			clear_flags = (step_cnt % 64) == 0;
			clear_flags = true;
			if (root->check_for_refine(clear_flags)) {
				solve(ELL_TOLER / 10.0);
				root->store_phi();
			}
			root->sub_pot_et();
		}
		step_cnt++;
		root->set_time(get_time() + dt);
		state_sums_out((root->state_sum() + root->flow_off())[0]);
	}
};

#endif /* HYDRO_DRIVER_H_ */
#endif
