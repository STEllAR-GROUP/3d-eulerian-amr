#include "../defs.h"
#ifdef GR
#include "initialize.h"
#include "../grid/grid_node.h"
#include "state.h"
#include "math.h"
#include "../binary/lane_emden.h"

void initialize(GridNode* g) {
	State::turn_fluid_on();
	State U = Vector<Real, STATE_NF>(0.0);
	Real r, d, r1, d0, alpha;
	alpha = 1.0 / 15.0;
	d0 = pow(alpha * alpha / POLY_K / 5.0, -3.0);
	State::d_floor = d0 / 1.0e+10;
	for (int k = 0; k < GNX; k++) {
		for (int j = 0; j < GNX; j++) {
			for (int i = 0; i < GNX; i++) {
				U = Vector<Real, STATE_NF>(0.0);
				r = sqrt(g->xc(i) * g->xc(i) + g->yc(j) * g->yc(j) + g->zc(k) * g->zc(k));
				r1 = r / alpha;
				if (r1 < 3.7) {
					d = lane_emden(r1);
					U.d() = max(d0 * d, State::d_floor);
				} else {
					U.d() = State::d_floor;
				}
				U.tau() = POLY_K * pow(U.d(), State::fgamma) / (State::fgamma - 1.0);
				_3Vec X;
				X[0] = g->xc(i);
				X[1] = g->yc(j);
				X[2] = g->xc(k);
				U.floor(X);
				(*g)(i, j, k) = U;
			}
		}
	}
	State::turn_fluid_off();
}
#endif

