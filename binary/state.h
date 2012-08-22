#ifndef EULER_STATE_H_
#define EULER_STATE_H_
#include "defs.h"
#include "real.h"
#include "state.h"
#include "vector.h"
#include "oct_node/oct_face.h"

#define NRHO 2

#define G_CON (1.0)

#define STATE_NF (6+NRHO)

class State: public Vector<Real, STATE_NF> {
public:
#ifdef CENTER_OF_MASS_CORRECTION
	static Real fR, ftheta;
#endif
	static const Real gamma;
	static Real rho_floor;
	static const Real ei_floor;
	static const int d_index;
	static const int sr_index;
	static const int lz_index;
	static const int sz_index;
	static const int pot_index;
	static const int et_index;
	static const int tau_index;
	static bool smooth_variable(int i) {
		return (pot_index == i);
	}
	static bool low_order_variable(int i) {
		return (pot_index == i);
	}
	State();
	State(const Vector<Real, STATE_NF>&);
	virtual Vector<Real, STATE_NF> source(const _3Vec&) const;
	static Vector<Real, STATE_NF> source_x(const State& fl, const _3Vec&);
	static Vector<Real, STATE_NF> source_y(const State& fl, const _3Vec&);
	static Vector<Real, STATE_NF> source_z(const State& fl, const _3Vec&);
	virtual Vector<Real, STATE_NF> scf_source(Real, Real, int) const;
	virtual Real max_abs_x_eigen(const _3Vec& X) const;
	virtual Real max_abs_y_eigen(const _3Vec& X) const;
	virtual Real max_abs_z_eigen(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
	virtual Real poisson_source() const;
	virtual void to_prim(const _3Vec& x);
	virtual void to_con(const _3Vec& x);
	static const char* field_name(int i);
	virtual void floor(const _3Vec& X);
	Real enthalpy() const;
	static Real vx0(const _3Vec& X);
	static Real vy0(const _3Vec& X);
	Real rho() const;
	Real sr() const;
	Real st(const _3Vec&) const;
	Real lz() const;
	static Real R(const _3Vec& X) {
		return sqrt(X[0] * X[0] + X[1] * X[1]);
	}
	Real sz() const;
	Real vx(const _3Vec& X) const;
	Real vy(const _3Vec& X) const;
	Real vz(const _3Vec& X) const;
	Real lz(const _3Vec& X) const;
	_3Vec V(const _3Vec& X) const;
	Real ek(const _3Vec& X) const;
	Real ei(const _3Vec& X) const;
	Real et() const;
	Real cs(const _3Vec& X) const;
	Real pg(const _3Vec& X) const;
	Real pot() const;
	Real get_rho(int i) const;
	void set_et(Real e) {
		(*this)[et_index] = e;
	}
	void set_pot(Real);
	void set_rho(Real, int = 0);
	void set_sx(Real);
	void set_sy(Real);
	void set_sz(Real);
	void set_tau();
	Real et_inertial(const _3Vec& X) const;
	Real pot_inertial(const _3Vec& X) const;
	void reflect_on_z();
	virtual void enforce_outflow(const _3Vec& X, const OctFace&);
	Vector<Real, STATE_NF> gravity_source(const _3Vec& g, Vector<Real, STATE_NF>& D, const _3Vec&) const;
};

#endif /* EULER_STATE_H_ */
