#ifndef EULER0_STATE_H_
#define EULER0_STATE_H_

#include "defs.h"
#include "real.h"
#include "state.h"
#include "vector.h"
#include "oct_node/oct_face.h"

#ifndef NRHO
#define NRHO 1
#endif

#define G_CON (1.0)

#define STATE_NF (6)

class State: public Vector<Real, STATE_NF> {
public:
	static const Real gamma;
	static const Real ei_floor;
	static const Real rho_floor;
	static const int d_index;
	static const int sx_index;
	static const int sy_index;
	static const int sz_index;
	static const int et_index;
	static const int tau_index;
	State();
	State(const Vector<Real, STATE_NF>&);
	virtual Real refine_value() const;
	virtual ~State() {
		return;
	}
	virtual void initialize() {
	}
	static bool low_order_variable(int i) {
		return false;
	}
	static bool smooth_variable(int i) {
		return false;
	}
	virtual Vector<Real, STATE_NF> source(const _3Vec& X) const;
	virtual Real max_abs_x_eigen(const _3Vec& X) const;
	virtual Real max_abs_y_eigen(const _3Vec& X) const;
	virtual Real max_abs_z_eigen(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
	virtual Real poisson_source() const;
	virtual void to_prim(const _3Vec& X);
	virtual void to_con(const _3Vec& X);
	static const char* field_name(int i);
	virtual void floor(const _3Vec&);
	Real rho() const;
	Real tau() const;
	Real sz() const;
	Real sx() const;
	Real sy() const;
	Real vx() const;
	Real vy() const;
	Real vz() const;
	Real ek() const;
	_3Vec V() const {
		_3Vec v;
		v[0] = vx();
		v[1] = vy();
		v[2] = vz();
		return v;
	}
	Real et() const;
	Real ei() const;
	Real cs() const;
	Real pg() const;
	Real pot() const;
	Real get_rho(int i) const {
		if (i < 0) {
			return rho();
		} else {
			return (*this)[d_index + i];
		}
	}
	void set_pot(Real);
	void set_rho(Real, int = 0);
	void set_sx(Real);
	void set_sy(Real);
	void set_sz(Real);
	void set_et(Real);
	void set_tau(Real);
	void reflect_on_z();
	virtual void enforce_outflow(const _3Vec&,const OctFace&);
};

#endif /* EULER_STATE_H_ */
