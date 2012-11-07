#ifndef OPTIONS_SCF_
#define OPTIONS_SCF_
#include "../real.h"
//#define SCF_CODE
#define ZTWD

#define CENTER_OF_MASS_CORRECTION
#define DYNAMIC_OMEGA
#define MIRROR_REFINE_Y


#define SAVE_RECONSTRUCTIONS
#define MAX_POINTS (1024*1024)
#define DYNAMIC_RANGE 1.0e+6
#ifndef SCF_CODE
#define READ_FROM_FILE "X.start.chk"
//#define READ_SILO
#endif
//#define BENCH 64
//#define PAR0
#define GRID_DIM            (1.0)
#define GRAVITY_EULER_STATE
#define USE_POISSON_DRIVE
//#define USE_HYDRO_DRIVE
#define MINMOD_THETA        1.3
#define Z_REFLECT
//#define POISSON_TEST
#define FRAME_RATE 100
#define PPM
#define BW 					(3)
#define EULER_GAMMA         (5.0/3.0)
#ifdef ZTWD
#define GNX 				(12+2*BW)
#else
#define GNX 				(12+2*BW)
#endif
#define RK_ORDER          	(3)
#define MAXLEVEL           	(4)
#define TIME_MAX           	(1000000.0)
//#define OUTPUT_STEP_FREQ    1
#define MAX_GRIDS           1024
#define MIN_GRIDS            64
//#define OUTPUT_TIME_FREQ   	(2.0*M_PI/1.097069/100.0)
#define OUTPUT_TIME_FREQ   	(1.0)
#define DEN_REF
#define REL_REF
#define VIRIAL_TOLERANCE 1.0e-6

#define KL  KR
#define KR 0.130821412355133


#define MAXDTINC           (1.25)
#ifdef SCF_CODE
#define GRID_CFL_FACTOR    0.3
#define CHKPT_FREQ 10
#else
#define CHKPT_FREQ         (32)
#define GRID_CFL_FACTOR    0.3
#endif


#define MAXINITDT          (5.0e-3)
#define ELL_TOLER          (1.0e-7)
#define ELL_TOLER2          (1.0)


void init_scf(const char*);
Real scf_phi(Real, Real, Real, Real);
Real scf_rho(Real, Real, Real, Real);
Real scf_et(Real, Real, Real, Real);
void scf_com(Real, Real, Real);


#endif
