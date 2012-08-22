#include "../defs.h"
#ifdef SINGLE
#include "initialize.h"
#include "indexer_2d.h"
#include <math.h>
extern "C" {
#include "../binary/lane_emden.h"
}
#include "../oct_node/oct_node.h"
#include "../grid/grid_node.h"
#include "../binary/binary.h"

void initialize(GridNode* g) {
	Real R1 = 0.25;
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	const Real BinaryK = pow(R1 / 3.65375, 2) * 4.0 * M_PI / (2.5) * pow(1.0, 1.0 / 3.0);
	const Real E1 = sqrt((4.0 * 4.0 * atan(1.0) * pow(1.0, 1.0 - 1.0 / 1.5)) / ((1.5 + 1.0) * BinaryK));
	for (int index = 0; index <= indexer.max_index(); index++) {
		State U;
		Real x, y, z;
		U.set_rho(0.0);
		U.set_et(0.0);
		U.set_sz(0.0);
		int M0;
		if (g->is_finest()) {
			M0 = g->get_dx() / R1 + 2.0;
		} else {
			M0 = 2;
		}
		Real rho, r, x1, y1;
		x1 = 0.0;
		y1 = 0.0;
		int i0, j0, k0, j, k;
		j = indexer.x(index);
		k = indexer.y(index);
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i0,j0,k0,rho,r,x,y,z)
		for (int i = 0; i < GNX; i++) {
			y = minmod(g->yf(j) - y1, g->yf(j + 1) - y1);
			z = minmod(g->zf(k), g->zf(k + 1));
			x = minmod(g->xf(i) - x1, g->xf(i + 1) - x1);
			r = sqrt(x * x + y * y + z * z);
			r *= E1;
			rho = 0.0;
			if (r < 3.7) {
				for (k0 = 0; k0 < M0; k0++) {
					z = g->zf(k) + (Real(k0) + 0.5) * g->get_dx() / Real(M0);
					for (j0 = 0; j0 < M0; j0++) {
						y = g->yf(j) + (Real(j0) + 0.5) * g->get_dx() / Real(M0);
						for (i0 = 0; i0 < M0; i0++) {
							x = g->xf(i) + (Real(i0) + 0.5) * g->get_dx() / Real(M0);
							r = sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1) + z * z);
							r *= E1;
							if (r < 3.7) {
								rho += pow(lane_emden(r), 1.5) / pow(M0, 3);
							}
						}
					}
				}
			}
			x = g->xc(i);
			y = g->yc(j);
			z = g->zc(k);
			Real ei = BinaryK * pow(rho, 1.0 + 1.0 / 1.5) / (State::gamma - 1.0);
			rho = max(rho, State::rho_floor);
			(*g)(i, j, k).set_rho(rho, 0);
			(*g)(i, j, k).set_et(ei);
			(*g)(i, j, k).set_sx(0.0);
			(*g)(i, j, k).set_sy(0.0);
			(*g)(i, j, k).set_sz(0.0);
			_3Vec X;
			X[0] = x;
			X[1] = y;
			X[2] = z;
			(*g)(i, j, k).set_tau(pow(ei, 1.0 / State::gamma));
			(*g)(i, j, k).floor(X);
		}
	}
}
#endif
