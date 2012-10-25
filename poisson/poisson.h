#ifndef GRIDPOISSON_H_
#define GRIDPOISSON_H_

#define PNX (GNX-2*BW+2)
#define HIORDER

#include "grid/grid_node.h"
#include "poisson_interface.h"
#include "poisson_amr_bound.h"
#include "poisson_phys_bound.h"

typedef double PReal;

class Poisson: public GridNode, public PoissonInterface {
private:

	Array3d<PReal, PNX, PNX, PNX> phi0;
	Array3d<PReal, PNX, PNX, PNX> dphi1;
	Array3d<PReal, PNX, PNX, PNX> dphi;
	Real last_resid, resid;
protected:
	Array3d<PReal, PNX, PNX, PNX> S;
	Array3d<PReal, PNX, PNX, PNX> phi;
	Array3d<PReal, PNX, PNX, PNX> fx;
	Array3d<PReal, PNX, PNX, PNX> fy;
	Array3d<PReal, PNX, PNX, PNX> fz;
	void adjust_x_force();
	void adjust_y_force();
	void adjust_z_force();
	void adjust_x_dphi();
	void adjust_y_dphi();
	void adjust_z_dphi();
	void enforce_phi_boundaries(int);
	void enforce_phi_boundaries_edge(int);
	void enforce_phi_boundaries_vertex(int);
	void enforce_dphi_boundaries(int);
	void init_vcycle();
	void vcycle_down(int);
	void vcycle_coarse_correction(int);
	Real vcycle_relax(int);
	void vcycle_retire_dphi(int);
	bool poisson_zone_is_refined(int, int, int) const;
	virtual void inject_from_parent(ChildIndex);
public:
	void phi_from_children();
	void restore_phi_n(Real dt);
	void store_phi_n12(Real dt);
	void store_phi_n(Real dt);
	void set_source();
	void store_phi();
	void restore_half_phi();
	void restore_phi();
	Real solution_error();
	Real pxc(int) const;
	Real pyc(int) const;
	Real pzc(int) const;
	Poisson();
	virtual Real compute_phi(Real, Real, Real) const = 0;
	virtual Poisson* new_octnode() const = 0;
	void compute_forces();
	virtual PoissonPhysBound* new_phys_bound() const;
	virtual PoissonAMRBound* new_amr_bound() const;
	virtual ~Poisson();
	Real error();
	virtual Real get_phi(int, int, int) const;
	virtual Real get_dphi(int, int, int) const;
	Real vcycle();
	void vcycle0();
};

#endif /* GRIDPOISSON_H_ */
