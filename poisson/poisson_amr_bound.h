

#ifndef POISSON_AMR_BOUND_H_
#define POISSON_AMR_BOUND_H_

#include "grid/grid_amr_bound.h"
#include "poisson_interface.h"

class Poisson;

class PoissonAMRBound: public GridAMRBound, public PoissonInterface {
private:
	Vector<int, 3> rel_offset;
	const Poisson* sib;
	OctFace face;
public:
	PoissonAMRBound();
	virtual ~PoissonAMRBound();
	virtual Real get_phi(int, int, int) const;
	virtual Real get_dphi(int, int, int) const;
	virtual void set_sibling(OctFace, OctNode*);
};

#endif /* POISSON_AMR_BOUND_H_ */
