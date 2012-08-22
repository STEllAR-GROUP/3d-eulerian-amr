#ifndef OPTIONS_____H
#define OPTIONS_____H

#define BINARY

#ifdef SINGLE
#include "./single/single_defs.h"
#endif

#ifdef BINARY
#include "./binary/binary_defs.h"
#endif

#ifdef EULER
#include "./euler/euler_defs.h"
#endif

#ifndef OMP_SCHEDULE
#define OMP_SCHEDULE        static
#endif
#define DEBUG_V 30
#define MINLEVEL           	(1)
#define MIN_GRID_LIFE       (RK_ORDER*(GNX-2*BW))
#define MAX_POISSON_ITER    10000

#endif

#ifndef GRID_DIM
#define GRID_DIM 1.0
#endif

