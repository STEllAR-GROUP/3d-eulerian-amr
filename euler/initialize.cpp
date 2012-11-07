
#include "../defs.h"
#ifdef EULER
#include "initialize.h"
#include "../grid/grid_node.h"
#include "state.h"
#include "math.h"

void initialize(GridNode* g) {
	State U = Vector<Real, STATE_NF> (0.0);
	Real r;
	const Real rho0 = 1.0;
	const Real ei0 = 1.0;
	const Real tau0 = pow(ei0, 1.0 / State::gamma);
	const Real rho1 = 1.0e-10;
	const Real ei1 = 1.0e-10;
	const Real tau1 = pow(ei1, 1.0 / State::gamma);


	const Real kappa = 1.0;
	const Real G = 1.0;
	const Real M_c = 2e-2;

	const Real eps = 0.4;
	const Real R_outer = 1.0747e-4;
	const Real R_inner = r_outer*(1.0-eps)/(1.0+eps);
	const Real C_2 = sqrt(2.0*G*M_c*R_inner*R_outer/(R_inner+R_outer));
	const Real C_1 = pow(0.5*(C_2/R_inner),2.0)-G*M_c/R_inner;

	for (int k = 0; k < GNX; k++) {
		for (int j = 0; j < GNX; j++) {
			for (int i = 0; i < GNX; i++) {
			  r = sqrt(g->xc(i) * g->xc(i) + g->yc(j) * g->yc(j));
			  if ( (r_inner <= r) && (r_outer >= r) ) {
			    U.set_rho( (0.5/kappa)*(C_1+G*M_c/r-0.5*pow((C_2/r),2.0) );
			    U.set_et(ei0);
			    U.set_tau(tau0);
			  } else {
			    U.set_rho(rho1);
			    U.set_et(ei1);
			    U.set_tau(tau1);
			  }
			  (*g)(i, j, k) = U;
			}
		}
	}
}
#endif
