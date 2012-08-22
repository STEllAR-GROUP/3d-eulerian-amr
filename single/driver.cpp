#include "../defs.h"
#ifdef SINGLE
#include "driver.h"
#include "initialize.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>



void SingleDriver::state_sums_out(const State& s) {
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
#ifndef SINGLE
#ifdef POLYTROPIC
	State::rho_floor = max(s0[0], s0[1]) / 1.0e+10;
#else
       State::rho_floor = max(s0[0], s0[1]) / 1.0e+17;
  #endif
#endif
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

#endif
