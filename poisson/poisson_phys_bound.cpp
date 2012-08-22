#include "poisson_phys_bound.h"
#include "indexer_2d.h"
#include "poisson.h"
#include <stdlib.h>

PoissonPhysBound::PoissonPhysBound() :
	GridPhysBound(), PoissonInterface() {
}

void PoissonPhysBound::set_sibling(OctFace f, OctNode* sib0) {
	GridPhysBound::set_sibling(f, sib0);
	const Poisson* sib = static_cast<const Poisson*> (sib0);
	const Indexer2d indexer(0, PNX - 1, 0, PNX - 1);
	int j, k;
	reflect = false;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		switch (f) {
		case XU:
#ifdef X_REFLECT
			reflect = true;
#endif
			direction = 0;
			x[j][k] = sib->pxc(0);
			y[j][k] = sib->pyc(j);
			z[j][k] = sib->pzc(k);
			break;
		case YU:
#ifdef Y_REFLECT
			reflect = true;
#endif
			direction = 1;
			x[j][k] = sib->pxc(j);
			y[j][k] = sib->pyc(0);
			z[j][k] = sib->pzc(k);
			break;
		case ZU:
#ifdef Z_REFLECT
			reflect = true;
#endif
			direction = 2;
			x[j][k] = sib->pxc(j);
			y[j][k] = sib->pyc(k);
			z[j][k] = sib->pzc(0);
			break;
		case XL:
			direction = 0;
			x[j][k] = sib->pxc(PNX - 1);
			y[j][k] = sib->pyc(j);
			z[j][k] = sib->pzc(k);
			break;
		case YL:
			direction = 1;
			x[j][k] = sib->pxc(j);
			y[j][k] = sib->pyc(PNX - 1);
			z[j][k] = sib->pzc(k);
			break;
		case ZL:
			direction = 2;
			x[j][k] = sib->pxc(j);
			y[j][k] = sib->pyc(k);
			z[j][k] = sib->pzc(PNX - 1);
			break;
		}
	}
}

void PoissonPhysBound::compute_boundaries(const Poisson* root) {
	const Indexer2d indexer(0, PNX - 1, 0, PNX - 1);
	int j, k;
	if (!reflect) {
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k)
		for (int index = 0; index <= indexer.max_index(); index++) {
			j = indexer.x(index);
			k = indexer.y(index);
			phi[j][k] = root->compute_phi(x[j][k], y[j][k], z[j][k]);
		}
	}
}

PoissonPhysBound::~PoissonPhysBound() {
}

Real PoissonPhysBound::get_phi(int i, int j, int k) const {
	if (!reflect) {
		if (direction == 0) {
			return phi[j][k];
		} else if (direction == 1) {
			return phi[i][k];
		} else {
			assert(direction==2);
			return phi[i][j];
		}
	} else {
		if (direction == 0) {
			if (i == 1) {
				return static_cast<const Poisson*> (sibling)->get_phi(PNX - 2, j, k);
			} else if (i == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_phi(1, j, k);
			} else {
				assert(false);
			}
		} else if (direction == 1) {
			if (j == 1) {
				return static_cast<const Poisson*> (sibling)->get_phi(i, PNX - 2, k);
			} else if (j == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_phi(i, 1, k);
			} else {
				assert(false);
			}
		} else if (direction == 2) {
			if (k == 1) {
				return static_cast<const Poisson*> (sibling)->get_phi(i, j, PNX - 2);
			} else if (k == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_phi(i, j, 1);
			} else {
				assert(false);
			}
		}
		printf("Error 17\n");
		abort();
	}
	return 0.0;
}

Real PoissonPhysBound::get_dphi(int i, int j, int k) const {
	if (!reflect) {
		return 0.0;
	} else {
		if (direction == 0) {
			if (i == 1) {
				return static_cast<const Poisson*> (sibling)->get_dphi(PNX - 2, j, k);
			} else if (i == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_dphi(1, j, k);
			}
		} else if (direction == 1) {
			if (j == 1) {
				return static_cast<const Poisson*> (sibling)->get_dphi(i, PNX - 2, k);
			} else if (j == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_dphi(i, 1, k);
			}
		} else if (direction == 2) {
			if (k == 1) {
				return static_cast<const Poisson*> (sibling)->get_dphi(i, j, PNX - 2);
			} else if (k == PNX - 2) {
				return static_cast<const Poisson*> (sibling)->get_dphi(i, j, 1);
			}
		}
	}
	printf("Error 18\n");
	abort();
}
