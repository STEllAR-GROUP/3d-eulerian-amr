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
	static const Real gamma;
	static const Real ei_floor;
	static const int d_index;
	static const int sr_index;
	static const int lz_index;
	static const int sz_index;
	static const int pot_index;
	static const int et_index;
	static const int tau_index;
	static Real fR, ftheta;
	static Real rho_floor;
	static bool low_order_variable(int i);
	static bool smooth_variable(int i);
	static const char* field_name(int i);
	virtual Vector<Real, STATE_NF> source(const _3Vec&) const;
	virtual Real max_abs_x_eigen(const _3Vec& X) const;
	virtual Real max_abs_y_eigen(const _3Vec& X) const;
	virtual Real max_abs_z_eigen(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
	virtual Real poisson_source() const;
	virtual void to_prim(const _3Vec& x);
	virtual void enforce_outflow(const _3Vec& X, const OctFace&);
	virtual void to_con(const _3Vec& x);
	virtual void floor(const _3Vec& X);
	State();
	State(const Vector<Real, STATE_NF>&);
	_3Vec V(const _3Vec& X) const;
	Real et_inertial(const _3Vec& X) const;
	Real pot_inertial(const _3Vec& X) const;
	Real pot() const;
	Real rho(int i) const;
	Real enthalpy() const;
	Real rho() const;
	Real et() const;
	Real lz() const;
	Vector<Real, STATE_NF> scf_source(Real, Real, int) const;
	Vector<Real, STATE_NF> gravity_source(const _3Vec& g, Vector<Real, STATE_NF>& D, const _3Vec&) const;
	void reflect_on_z();
private:
	static Real R(const _3Vec& X);
	Real dvx(const _3Vec& X) const;
	Real dvy(const _3Vec& X) const;
	Real sr() const;
	Real st(const _3Vec&) const;
	Real sz() const;
	Real vx(const _3Vec& X) const;
	Real vy(const _3Vec& X) const;
	Real vz(const _3Vec& X) const;
	Real ek(const _3Vec& X) const;
	Real ei(const _3Vec& X) const;
	Real cs(const _3Vec& X) const;
	Real pg(const _3Vec& X) const;
};

#endif /* EULER_STATE_H_ */
