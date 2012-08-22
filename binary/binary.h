/*
 * Binary.h
 *
 *  Created on: May 8, 2012
 *      Author: dmarce1
 */

#ifndef BINARY_H_
#define BINARY_H_

#include "../poisson/poisson.h"

class Binary: public Poisson {
public:
	static Real phase;
	static Real Omega0;
	static Real M1;
	static Real R1;
	static Real n;
	static Real a;
	static Real q;
	static Real fill_factor;
	static Real Omega;
	static Real K1;
	static Real K2;
	static Real phi0_1;
	static Real phi0_2;
	static Real l1x;
	static Real scf_dt;
	static Real R2;
private:
	Real src_0[8], src_x[8], src_y[8], src_z[8];
	bool all_refined;protected:
	virtual void inject_from_parent(ChildIndex);
	Array3d<Real, GNX, GNX, GNX> lobe;
public:
	struct binary_integrals_t {
		Real m1, m2, mc, js1, js2, jorb, jc, a, q, total_energy, kinetic, I1, I2, Ic, V1, V2;
		_3Vec x, xdot;
	};
	virtual const char* output_field_names(int) const;
	virtual int nvar_output() const;
	virtual void load_output(grid_output_t* go, int, int, int) const;
	Real omega() const;
	Real Iz() const;
	Real lz() const;
	Real find_K(Real, int) const;
	Real find_phimax(Real, int) const;
	Real find_phimin(Real, int) const;
	Real next_omega(Real* f = NULL, Real* df = NULL) const;
	void mark_lobes(Real, Real, Real, Real,Real,Real, Real);
	Vector<Real, 4> find_l1(const Vector<Real, 4>& m1, const Vector<Real, 4>& m2, Vector<Real, 4>* = NULL) const;
	Real domega() const;
	Real sum_dlz() const;
	Real sum_virial(Real* norm = NULL) const;
	void M1M2data(binary_integrals_t*);
	Real sum_lz() const;
	Real pot_sum() const;
	void pot_to_grid();
	void pot_from_grid();
	void compute_physical_boundaries(Binary*);
	virtual Real compute_phi(Real, Real, Real) const;
	void compute_source_sums();
	void add_pot_et();
	void sub_pot_et();
	Binary();
	virtual ~Binary();
	virtual Binary* new_octnode() const;
	virtual void add_difs(Real, Real);
	virtual _3Vec V(int, int, int) const;
	virtual _3Vec Vfx(int, int, int) const;
	virtual _3Vec Vfy(int, int, int) const;
	virtual _3Vec Vfz(int, int, int) const;
	virtual Real max_dt() const;
	virtual void write(FILE*) const;
	virtual void read(FILE*);
};

#endif /* BINARY_H_ */
