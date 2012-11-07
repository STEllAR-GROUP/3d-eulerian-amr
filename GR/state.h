#ifndef EULER0_STATE_H_
#define EULER0_STATE_H_

#include "defs.h"
#include "real.h"
#include "state.h"
#include "vector.h"
#include "oct_node/oct_face.h"

#define STATE_NF (36)
#define PRIM_NF 5

#define POLY_K 1.0e-4

class State: public Vector<Real, STATE_NF> {
private:
#ifndef NDEBUG
	bool prims_called;
#endif
	static const int symtable[3][3];
	static const Real delta[3][3];
	static const int Theta_index;
	static const int A_index;
	static const int Z_index;
	static const int K_index;
	static const int D_index[3];
	static const int tau_index;
	static const int d_index;
	static const int s_index;
	static const int W_index;
	static const int rho_index;
	static const int v_index;
	static const int eta_index;
	static bool fluid_on;
	Vector<Real, STATE_NF> flux(int direction) const;
	Vector<Real, PRIM_NF> prims;
public:
	static const Real fgamma;
	virtual State& operator=(const State& v) {
		Vector<Real, STATE_NF>::operator=(Vector<Real, STATE_NF>(v));
		prims = v.prims;
#ifndef NDEBUG
		prims_called = true;
#endif
		return *this;
	}
	static Real d_floor;
	static void turn_fluid_on() {
		fluid_on = true;
	}
	static void turn_fluid_off() {
		fluid_on = false;
	}
	Real eta() const {
		return prims[eta_index];
	}
	Real W() const {
		Real v2 = 0.0;
		for (int i = 0; i < 3; i++) {
			v2 += v(i) * v(i);
		}
		return sqrt(1.0 - v2);
	}
	Real h() const {
		assert( prims_called );
		return 1.0 + eta() + pg() / rho();
	}
	Real Christoffel(int k, int i, int j) const {
		Real c;
		if (k == 3) {
			if (i == 3) {
				if (j == 3) {
					c = TrK() - 2.0 * Theta();
				} else {
					c = A(j);
				}
			} else {
				if (j == 3) {
					c = A(i);
				} else {
					c = -K(i, j);
				}
			}
		} else {
			if (i == 3) {
				if (j == 3) {
					c = A(k);
				} else {
					c = -K(k, j);
				}
			} else {
				if (j == 3) {
					c = -K(k, i);
				} else {
					c = D(i, j, k) + D(j, i, k) - D(k, i, j);
				}
			}
		}
		return c;
	}
	Real pg() const {
		return rho() * eta() * (fgamma - 1.0);
	}
	Real d() const {
		return max((*this)[d_index], d_floor);
	}
	Real s(int i) const {
		return (*this)[s_index + i];
	}
	void compute_prims();
	Real v(int i) const {
		return prims[v_index + i];
	}
	Real& d() {
		return (*this)[d_index];
	}
	Real& s(int i) {
		return (*this)[s_index + i];
	}
	Real Theta() const {
		return (*this)[Theta_index];
	}
	Real A(int i) const {
		return (*this)[A_index + i];
	}
	Real Z(int i) const {
		return (*this)[Z_index + i];
	}
	Real K(int i, int j) const {
		return (*this)[K_index + symtable[i][j]];
	}
	Real D(int i) const {
		return D(i, 0, 0) + D(i, 1, 1) + D(i, 2, 2);
	}
	Real E(int i) const {
		return D(0, i, 0) + D(1, i, 1) + D(2, i, 2);
	}
	Real lambda(int k, int i, int j) const {
		Real l;
		const Real zeta = 0.0;
		l = 2.0 * D(k, i, j);
		l += delta[i][k] * (A(j) + D(j) - 2.0 * E(j) - 2.0 * Z(j));
		l += delta[j][k] * (A(i) + D(i) - 2.0 * E(i) - 2.0 * Z(i));
		l -= (1.0 + zeta) * (D(i, j, k) + D(j, i, k));
		l -= (1.0 + zeta) * (-delta[i][k] * E(j) - delta[j][k] * E(i));
		l *= 0.5;
		return l;
	}
	Real D(int k, int i, int j) const {
		return (*this)[D_index[k] + symtable[i][j]];
	}
	Real& Theta() {
		return (*this)[Theta_index];
	}
	Real& A(int i) {
		return (*this)[A_index + i];
	}
	Real& Z(int i) {
		return (*this)[Z_index + i];
	}
	Real& K(int i, int j) {
		return (*this)[K_index + symtable[i][j]];
	}
	Real& D(int k, int i, int j) {
		return (*this)[D_index[k] + symtable[i][j]];
	}
	Real TrK() const {
		return K(0, 0) + K(1, 1) + K(2, 2);
	}
	Real TrS() const {
		return S(0, 0) + S(1, 1) + S(2, 2);
	}
	Real S(int i) const {
		return s(i);
	}
	Real S(int i, int j) const {
		if (i == 3) {
			if (j == 3) {
				return tau() + d();
			} else {
				return S(j);
			}
		} else {
			if (j == 3) {
				return S(i);
			} else {
				return d() * v(j) * v(i) * h() * W() + delta[i][j] * pg();
			}
		}
	}
	Real tau() const {
		return (*this)[tau_index];
	}
	Real& tau() {
		return (*this)[tau_index];
	}
	State();
	State(const Vector<Real, STATE_NF>&);
	virtual ~State() {
		return;
	}
	static bool low_order_variable(int i) {
		return false;
	}
	static bool smooth_variable(int i) {
		return false;
	}
	virtual Real poisson_source() const {
		return 0.0;
	}
	virtual Vector<Real, STATE_NF> source(const _3Vec& X) const;
	virtual Real max_abs_x_eigen(const _3Vec& X) const;
	virtual Real max_abs_y_eigen(const _3Vec& X) const;
	virtual Real max_abs_z_eigen(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
	virtual void to_prim(const _3Vec& X);
	virtual void to_con(const _3Vec& X);
	static const char* field_name(int i);
	virtual void floor(const _3Vec&);
	void reflect_on_z();
	virtual void enforce_outflow(const _3Vec&, const OctFace&);
	virtual Real rho() const;
};

#endif /* EULER_STATE_H_ */
