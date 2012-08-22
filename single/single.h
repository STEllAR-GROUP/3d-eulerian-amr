/*
 * Single.h
 *
 *  Created on: May 8, 2012
 *      Author: dmarce1
 */

#ifndef SINGLE_H_
#define SINGLE_H_

#include "../poisson/poisson.h"

class Single: public Poisson {
private:
	Real src_0[8], src_x[8], src_y[8], src_z[8];
	bool all_refined;
protected:
	virtual void inject_from_parent(ChildIndex);
public:
	struct Single_integrals_t {
		Real m1, m2, mc, js1, js2, jorb, jc, a, q, virial, kinetic, I1, I2, Ic;
		_3Vec x, xdot;
	};
	void pot_to_grid();
	Real pot_sum() const;
	void pot_from_grid();
	void compute_physical_boundaries(Single*);
	virtual Real compute_phi(Real, Real, Real) const;
	void compute_source_sums();
	void add_pot_et();
	void sub_pot_et();
	Single();
	virtual ~Single();
	virtual Single* new_octnode() const;
	virtual void add_difs(Real, Real);
	virtual Real max_dt() const;
	virtual void write(FILE*) const;
	virtual void read(FILE*);
	virtual void sum_x_difs();
	virtual void sum_y_difs();
	virtual void sum_z_difs();
};

#endif /* Single_H_ */
