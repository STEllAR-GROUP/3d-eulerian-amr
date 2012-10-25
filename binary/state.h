#ifndef EULER_STATE_H_
#define EULER_STATE_H_
#include <stdio.h>
#include "defs.h"
#include "real.h"
#include "state.h"
#include "vector.h"
#include "oct_node/oct_face.h"
#include "helmholtz/helmholtz.h"

#define G_CON (1.0)

#define NSPECIES 3

#define STATE_NF (NSPECIES+10)

class State: public Vector<Real, STATE_NF> {
public:
	static Real abar[NSPECIES];
	eos_t compute_eos(const _3Vec& X) const;
	static const Real gamma;
	static const Real ei_floor;
	static const int d_index;
	static const int d_acc_index;
	static const int d_don_index;
	static const int sr_index;
	static const int lz_index;
	static const int sz_index;
	static const int pot_index;
	static const int et_index;
	static const int tau_index;
	static const int T_index;
	static const int specie_index;
	static Real acc_rho_cut, don_rho_cut;
	static int donor_type, accretor_type;
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
#ifdef ZTWD
	Real A() const;
	Real Z() const;
#endif
	Real et_inertial(const _3Vec& X) const;
	Real pot_inertial(const _3Vec& X) const;
	Real pot() const;
	Real rho_accretor() const;
	Real rho_donor() const;
	Real enthalpy() const;
	Real rho() const;
	Real tau() const {
		return (*this)[tau_index];
	}
	Real et() const;
	Real lz() const;
	Vector<Real, STATE_NF> scf_source(Real, Real, int) const;
	Vector<Real, STATE_NF> gravity_source(const _3Vec& g, Vector<Real, STATE_NF>& D, const _3Vec&) const;
	Vector<Real, STATE_NF> internal_energy_source(const _3Vec& X, Real divU);
	void reflect_on_z();
	void adjust_rho_fluxes();
	Real specie_rho(int) const;
	Real vx(const _3Vec& X) const;
	Real vy(const _3Vec& X) const;
	Real vz(const _3Vec& X) const;
private:
	static Real R(const _3Vec& X);
	Real dvx(const _3Vec& X) const;
	Real dvy(const _3Vec& X) const;
	Real sr() const;
	Real st(const _3Vec&) const;
	Real sz() const;
	Real ek(const _3Vec& X) const;
	Real ei(const _3Vec& X) const;
	Real cs(const _3Vec& X) const;
	Real pg(const _3Vec& X) const;
};

#endif /* EULER_STATE_H_ */
