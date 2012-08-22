#include <fenv.h>
#include <stdio.h>
#include "program.h"
#include "single/single.h"
#include "single/driver.h"
#include "binary/binary_driver.h"
#include "euler/hydro_driver.h"
#include <stdlib.h>
#include <omp.h>

int Program::run(int argc, char* argv[]) {
	GridNode::max_refine_level = atoi(argv[1]);
	try {
#ifdef BINARY
		Real norb;
#ifdef SCF_CODE
			Real R2, q, a;
			if (argc < 5) {
				printf("Command format is scf <maxlevel> <mass ratio> <r2> <fill factor>\n");
				abort();
			}
			GridNode::max_refine_level = atoi(argv[1]);
			R2 = atof(argv[3]);
			Binary::q = q = atof(argv[2]);
			a = R2 / (0.49 * pow(q, 2.0 / 3.0) / (0.60 * pow(q, 2.0 / 3.0) + log(1.0 + pow(q, 1.0 / 3.0))));
			Binary::a = a;
			Binary::R2 = R2;
			printf("a=%e\n", Binary::a);
			Binary::fill_factor = atof(argv[4]);
			Binary::R1 = Binary::a / 8.0;
			Binary::Omega = sqrt((Binary::q + 1.0) / pow(Binary::a, 3.0));
			Binary::Omega0 = Binary::Omega;
#else
			if (argc < 3) {
				printf("Command format is amr <maxlevel> <norb>\n");
				abort();
			}
			GridNode::max_refine_level = atoi(argv[1]);
			norb = atof(argv[2]);
#endif
		assert( (GNX - 2*BW) % 4 == 0 );
		omp_set_dynamic(1);
		assert( GNX % 2 == 0 );

		Binary* root = new Binary;
		root->init();
		BinaryDriver* driver = new BinaryDriver(root);
		driver->set_output_frequency_by_time(OUTPUT_TIME_FREQ);
		driver->step_to_time(norb);
		delete root;
#endif

#ifdef SINGLE
		Single* root = new Single;
		root->init();
		SingleDriver* driver = new SingleDriver(root);
		driver->set_output_frequency_by_time(OUTPUT_TIME_FREQ);
		driver->set_rk(RK_ORDER);
		driver->step_to_time(TIME_MAX);
		delete root;
#endif

#ifdef EULER
		GridNode* root = new GridNode;
		root->init();
		HydroDriver* driver = new HydroDriver(root);
		driver->set_output_frequency_by_time(OUTPUT_TIME_FREQ);
		driver->set_rk(RK_ORDER);
		driver->step_to_time(TIME_MAX);
		delete root;
#endif

	} catch (...) {
		printf("Error 5\n");
		abort();
	}
	return 0;
}

Program::Program() {
	printf("Starting\n");
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#endif
}

Program::~Program() {
#ifndef NDEBUG
	printf("Terminating\n");
#endif
}
