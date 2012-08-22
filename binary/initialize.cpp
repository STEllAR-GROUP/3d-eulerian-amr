#include "../defs.h"
#ifdef BINARY
#include "initialize.h"
#include "indexer_2d.h"
#include <math.h>
extern "C" {
#include "lane_emden.h"
}
#include "../oct_node/oct_node.h"
#include "../grid/grid_node.h"
#include "../binary/binary.h"

void initialize(GridNode* g) {
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	for (int index = 0; index <= indexer.max_index(); index++) {
		const Real rho_c2 = Binary::q*Binary::M1 * (pow(3.65375 / Binary::R2, 3.0) / (2.71406 * 4.0 * M_PI));
		const Real rho_c1 = rho_c2 / pow(Binary::q, 2.0 * Binary::n / (3.0 - Binary::n));
		const Real BinaryK = pow(Binary::R2 / 3.65375, 2) * 4.0 * M_PI / (2.5) * pow(rho_c2, 1.0 / 3.0);
		Binary::K1 = BinaryK;
		Binary::K2 = BinaryK;
		const int j = indexer.x(index);
		const int k = indexer.y(index);
		const Real pi = 4.0 * atan(1.0);
		const Real E1 = sqrt((4.0 * pi * pow(rho_c1, 1.0 - 1.0 / Binary::n)) / ((Binary::n + 1.0) * BinaryK));
		const Real E2 = sqrt((4.0 * pi * pow(rho_c2, 1.0 - 1.0 / Binary::n)) / ((Binary::n + 1.0) * BinaryK));
		const Real x1 = Binary::a * Binary::q / (Binary::q + 1.0);
		const Real x2 = -Binary::a / (Binary::q + 1.0);
		//	printf("%e %e %e %e %e\n", rho_c1, rho_c2, BinaryK, E1, E2);
		State U;
		U.set_rho(0.0);
		U.set_sz(0.0);
		int M0;
		if (g->is_finest()) {
			M0 = g->get_dx() / Binary::R1 + 2.0;
		} else {
			M0 = 2;
		}
		Real r1, r2, x, y, z;
		Real dx = 0.0;
		Real dy = 0.0;
		Real rho1, rho2, rho, sx, sy;
		int i0, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i0,j0,k0,rho1,rho2,r1,r2,x,y,z,sx,sy,rho)
		for (int i = 0; i < GNX; i++) {
			y = minmod(g->yf(j), g->yf(j + 1));
			z = minmod(g->zf(k), g->zf(k + 1));
			x = minmod(g->xf(i) - x1, g->xf(i + 1) - x1);
			r1 = sqrt(x * x + y * y + z * z);
			r1 *= E1;
			x = minmod(g->xf(i) - x2, g->xf(i + 1) - x2);
			r2 = sqrt(x * x + y * y + z * z);
			r2 *= E2;
			rho1 = rho2 = 0.0;
			if (r1 < 3.7 || r2 < 3.7) {
				for (k0 = 0; k0 < M0; k0++) {
					z = g->zf(k) + (Real(k0) + 0.5) * g->get_dx() / Real(M0);
					for (j0 = 0; j0 < M0; j0++) {
						y = g->yf(j) + (Real(j0) + 0.5) * g->get_dx() / Real(M0);
						for (i0 = 0; i0 < M0; i0++) {
							x = g->xf(i) + (Real(i0) + 0.5) * g->get_dx() / Real(M0);
							r1 = sqrt((x - x1) * (x - x1) + y * y + z * z);
							r1 *= E1;
							r2 = sqrt((x - x2) * (x - x2) + y * y + z * z);
							r2 *= E2;
							if (r1 < 3.7) {
								rho1 += pow(lane_emden(r1), Binary::n) / pow(M0, 3) * rho_c1;
							}
							if (r2 < 3.7) {
								rho2 += pow(lane_emden(r2), Binary::n) / pow(M0, 3) * rho_c2;
							}
						}
					}
				}
			}
			x = g->xc(i);
			y = g->yc(j);
			z = g->zc(k);
			rho1 = max(rho1, State::rho_floor / 2.0);
			rho2 = max(rho2, State::rho_floor / 2.0);
			rho = rho1 + rho2;
			(*g)(i,j,k) = Vector<Real,STATE_NF>(0.0);
			(*g)(i, j, k).set_rho(rho1, 0);
			(*g)(i, j, k).set_rho(rho2, 1);
			(*g)(i, j, k).set_pot(0.0);
			_3Vec X;
			X[0] = x;
			X[1] = y;
			X[2] = z;
			(*g)(i, j, k).floor(X);
		}
	}
}
#endif
