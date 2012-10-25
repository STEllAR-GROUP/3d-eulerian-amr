#ifndef ___HELMHOLTZ___C
#define ___HELMHOLTZ___C

typedef struct {
	double p, e, cs, cv, abar, zbar, rho, T, s, dsdt;
} eos_t;

void helmholtz_eos( eos_t* );
void helmholtz_compute_T( eos_t*, bool=false );
void helmholtz_ztwd( double*, double*, double, double );
void helmholtz_set_cgs_units( double cm, double g, double s, double T );

#define ZTWD

#endif
