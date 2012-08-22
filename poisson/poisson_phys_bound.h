#ifndef POISSON_PHYS_BOUND_H_
#define POISSON_PHYS_BOUND_H_

#include "../grid/grid_phys_bound.h"
#include "poisson_interface.h"

#define PNX (GNX-2*BW+2)

class Poisson;

class PoissonPhysBound: public GridPhysBound, public PoissonInterface {
private:
	bool reflect;
	Real phi[PNX][PNX];
	Real x[PNX][PNX];
	Real y[PNX][PNX];
	Real z[PNX][PNX];
	int direction;
public:
	PoissonPhysBound();
	void compute_boundaries(const Poisson*);
	virtual ~PoissonPhysBound();
	virtual Real get_phi(int, int, int) const;
	virtual Real get_dphi(int, int, int) const;
	virtual void set_sibling(OctFace, OctNode*);
};

#endif /* POISSON_PHYS_BOUND_H_ */
