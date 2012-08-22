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
	for (int k = 0; k < GNX; k++) {
		for (int j = 0; j < GNX; j++) {
			for (int i = 0; i < GNX; i++) {
				r = sqrt(g->xc(i) * g->xc(i) + g->yc(j) * g->yc(j) + g->zc(k) * g->zc(k));
				if (r < max(g->get_dx(), 0.005)) {
					U.set_rho(rho0);
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
