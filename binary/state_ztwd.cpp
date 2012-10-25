

#include "../defs.h"

#ifdef BINARY

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "oct_node/oct_face.h"
#include "state.h"
#include "helmholtz/helmholtz.h"
#include "binary.h"

#ifdef ZTWD

const Real State::gamma = EULER_GAMMA;
const Real State::ei_floor = 1.0e-20;
const int State::d_index = 0;
const int State::d_acc_index = 1;
const int State::d_don_index = 2;
const int State::sr_index = 3;
const int State::lz_index = 4;
const int State::sz_index = 5;
const int State::pot_index = 6;
const int State::et_index = 7;
const int State::tau_index = 8;
const int State::T_index = 9;
const int State::specie_index = 10;
Real State::abar[NSPECIES] = { 4.0, 14.12, 17.52 };

Real State::ftheta, State::fR;
Real State::rho_floor = 1.0e-10;

Real State::A() const {
	Real sum = 0.0;
	for (int i = 0; i < NSPECIES; i++) {
		sum += specie_rho(i) / abar[i];
	}
	sum = rho() / sum;
//	printf( "%e\n", sum );
	return sum;
}
Real State::Z() const {
	return A() / 2.0;
}

Real State::specie_rho(int i) const {
	return (*this)[specie_index + i];
}

bool State::low_order_variable(int i) {
	return (pot_index == i);
}

bool State::smooth_variable(int i) {
	return (pot_index == i);
}

const char* State::field_name(int i) {
	static char sstr[3];
	sstr[2] = '\0';
	sstr[0] = 's';
	assert(i >= 0);assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d_tot";
	case d_acc_index:
		return "d_acc";
	case d_don_index:
		return "d_don";
	case sr_index:
		return "sr";
	case lz_index:
		return "lz";
	case sz_index:
		return "sz";
	case pot_index:
		return "pot";
	case et_index:
		return "et";
	case tau_index:
		return "tau";
	case T_index:
		return "T";
	default:
		sstr[1] = '1' + i - specie_index;
		return sstr;
	}
}

Real State::R(const _3Vec& X) {
	return sqrt(X[0] * X[0] + X[1] * X[1]);
}

State::State() :
		Vector<Real, STATE_NF>() {
}

State::State(const Vector<Real, STATE_NF>& v) :
		Vector<Real, STATE_NF>(v) {
}

_3Vec State::V(const _3Vec& X) const {
	_3Vec v;
	v[0] = vx(X);
	v[1] = vy(X);
	v[2] = vz(X);
	return v;
}

Real State::cs(const _3Vec& X) const {
#ifdef ZTWD
	return compute_eos(X).cs;
#else
	return sqrt(gamma * pg(X) / rho());
#endif
}

Real State::enthalpy() const {
	Real h;
	Real x, h1, h2, A;
	if (rho_accretor() > rho_donor()) {
		A = Binary::K1;
	} else {
		A = Binary::K2;
	}
	x = pow(rho(), 1.0 / 3.0);
	h = 8.0 * A * (sqrt(x * x + 1.0) - 1.0);
	return h;
}

Real State::rho_accretor() const {
	return (*this)[d_acc_index];
}

Real State::rho_donor() const {
	return (*this)[d_don_index];
}

Real State::ei(const _3Vec& X) const {
	Real e;
	e = et() - ek(X);
	if (e < et() * 0.001) {
		return tau();
	} else {
		return e;
	}
#ifndef ZTWD
	if (e <= 0.001 * et()) {
		e = pow((*this)[tau_index], gamma);
	}
#else

#endif
	return max(ei_floor, e);
}

Real State::ek(const _3Vec& X) const {
	return 0.5 * (sr() * sr() + st(X) * st(X) + sz() * sz()) / rho();
}

Real State::et() const {
	return (*this)[et_index];
}

Real State::et_inertial(const _3Vec& X) const {
	const Real O = Binary::Omega0;
	return et() + O * (lz() - 0.5 * rho() * O * (X[0] * X[0] + X[1] * X[1]));
}

Real State::lz() const {
	return (*this)[lz_index];
}

#ifdef ZTWD
eos_t State::compute_eos(const _3Vec& X) const {
	eos_t eos;
	eos.rho = rho();
	eos.T = (*this)[T_index];
	eos.abar = A();
	eos.zbar = Z();
	helmholtz_eos(&eos);
	return eos;
}
#endif

Real State::max_abs_x_eigen(const _3Vec& X) const {
#ifndef SCF_CODE
	return fabs(vx(X)) + cs(X);
#else
	return 0.0;
#endif
}

Real State::max_abs_y_eigen(const _3Vec& X) const {
#ifndef SCF_CODE
	return fabs(vy(X)) + cs(X);
#else
	return 0.0;
#endif
}

Real State::max_abs_z_eigen(const _3Vec& X) const {
#ifndef SCF_CODE
	return fabs(vz(X)) + cs(X);
#else
	return 0.0;
#endif
}

Real State::pg(const _3Vec& X) const {
#ifdef ZTWD
	return compute_eos(X).p;
#else
	return (gamma - 1.0) * ei(X);
#endif
}

Real State::poisson_source() const {
	Real a = 4.0 * M_PI * G_CON * rho();
	return a;
}

Real State::pot() const {
	return (*this)[pot_index];
}

Real State::pot_inertial(const _3Vec& X) const {
	const Real O = Binary::Omega0;
	return pot() + 0.5 * rho() * O * O * (X[0] * X[0] + X[1] * X[1]);
}

Real State::rho() const {
	return (*this)[d_index];
}

Real State::sr() const {
	return (*this)[sr_index];
}

Real State::st(const _3Vec& X) const {
	return lz() / R(X) - rho() * R(X) * Binary::Omega;
}

Real State::sz() const {
	return (*this)[sz_index];
}

Real State::vx(const _3Vec& X) const {
	return (X[0] * sr() - X[1] * st(X)) / R(X) / rho();
}

Real State::vy(const _3Vec& X) const {
	return (X[1] * sr() + X[0] * st(X)) / R(X) / rho();
}

Real State::dvx(const _3Vec& X) const {
	return -X[1] * (Binary::Omega - Binary::Omega0);
}

Real State::dvy(const _3Vec& X) const {
	return X[0] * (Binary::Omega - Binary::Omega0);
}

Real State::vz(const _3Vec&) const {
	return sz() / rho();
}

Vector<Real, STATE_NF> State::internal_energy_source(const _3Vec& X, Real divU) {
	Vector<Real, STATE_NF> s = 0.0;
//	printf( "%e\n", divU );
	s[tau_index] = -compute_eos(X).p * divU;
	return s;
}

Vector<Real, STATE_NF> State::gravity_source(const _3Vec& g, Vector<Real, STATE_NF>& D, const _3Vec& X) const {
	State s = Vector<Real, STATE_NF>(0.0);
	const Real d = rho();
	s[sr_index] += d * (X[0] * g[0] + X[1] * g[1]) / R(X);
	s[lz_index] += d * (X[0] * g[1] - X[1] * g[0]);
	s[sz_index] += d * g[2];
	s[et_index] += D[pot_index];
	s[et_index] -= D[d_index] * pot() / d;
	s[et_index] += rho() * (dvx(X) * g[0] + dvy(X) * g[1]);
	return s;
}

Vector<Real, STATE_NF> State::scf_source(Real x, Real y, int lobe) const {
	Vector<Real, STATE_NF> s = 0.0;
	Real phi, rho_np1, h;
	const Real O2 = Binary::Omega * Binary::Omega;
	phi = pot() / rho(); // - 0.5 * (x * x + y * y) * O2;
	s[d_index] += -rho();
	s[d_acc_index] += -rho_accretor();
	s[d_don_index] += -rho_donor();
	for (int i = 0; i < NSPECIES; i++) {
		s[specie_index + i] += -specie_rho(i);
	}
	double minH;
	if (Binary::K1 > 0.0) {
		double dphi, cut;
		int type = 0;
		if (lobe == 2 && phi < Binary::phi0_2) {
			dphi = Binary::phi0_2 - phi;
			h = dphi / 8.0 / Binary::K2;
			rho_np1 = pow(h * h + 2.0 * h, 1.5);
			s[d_index] += rho_np1;
			s[d_don_index] += rho_np1;
			type = Binary::donor_type;
			cut = Binary::donor_cut;
		} else if (lobe == 1 && phi < Binary::phi0_1) {
			dphi = Binary::phi0_1 - phi;
			h = dphi / 8.0 / Binary::K1;
			rho_np1 = pow(h * h + 2.0 * h, 1.5);
			s[d_index] += rho_np1;
			s[d_acc_index] += rho_np1;
			type = Binary::accretor_type;
			cut = Binary::accretor_cut;
		}
		switch (type) {
		case 0:
			break;
		case 1:
			s[specie_index + 0] += rho_np1;
			break;
		case 2:
			if (rho() < cut) {
				s[specie_index + 0] += rho_np1;
			} else {
				s[specie_index + 1] += rho_np1;
			}
			break;
		case 3:
			s[specie_index + 1] += rho_np1;
			break;
		case 4:
			s[specie_index + 2] += rho_np1;
			break;

		}
	}
	return s;
}

Vector<Real, STATE_NF> State::source(const _3Vec& X) const {
	State s = Vector<Real, STATE_NF>(0.0);
#ifndef SCF_CODE
	s[sr_index] += (pg(X) + pow(lz() / R(X), 2) / rho()) / R(X);
	const Real theta = atan2(X[1], X[0]);
	s[sr_index] += rho() * fR * cos(theta - ftheta);
	s[lz_index] -= rho() * fR * sin(theta - ftheta) * R(X);
	s[et_index] += sr() * fR * cos(theta - ftheta);
	s[et_index] -= st(X) * fR * sin(theta - ftheta) * R(X);
#endif
	return s;
}

Vector<Real, STATE_NF> State::x_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
#ifndef SCF_CODE
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = vx(X) * (*this)[i];
	}
	flux[sr_index] += X[0] * pg(X) / R(X);
	flux[lz_index] -= X[1] * pg(X);
	flux[et_index] += (vx(X) + dvx(X)) * pg(X);
#endif
	return flux;
}

Vector<Real, STATE_NF> State::y_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
#ifndef SCF_CODE
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = vy(X) * (*this)[i];
	}
	flux[sr_index] += X[1] * pg(X) / R(X);
	flux[lz_index] += X[0] * pg(X);
	flux[et_index] += (vy(X) + dvy(X)) * pg(X);
#endif
	return flux;
}

Vector<Real, STATE_NF> State::z_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
#ifndef SCF_CODE
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = vz(X) * (*this)[i];
	}
	flux[sz_index] += pg(X);
	flux[et_index] += vz(X) * pg(X);
#endif
	return flux;
}

void State::enforce_outflow(const _3Vec& X, const OctFace& f) {
	(*this)[et_index] = ei(X);
	Real v0;
	switch (f) {
	case XU:
		if (vx(X) > 0.0) {
			v0 = vy(X);
			(*this)[lz_index] = (X[0] * (v0 + X[0] * Binary::Omega)) * rho();
			(*this)[sr_index] = X[1] * v0 * rho() / R(X);
		}
		break;
	case XL:
		if (vx(X) < 0.0) {
			v0 = vy(X);
			(*this)[lz_index] = (X[0] * (v0 + X[0] * Binary::Omega)) * rho();
			(*this)[sr_index] = X[1] * v0 * rho() / R(X);
		}
		break;
	case YU:
		if (vy(X) > 0.0) {
			v0 = vx(X);
			(*this)[lz_index] = (-X[1] * (v0 - X[1] * Binary::Omega)) * rho();
			(*this)[sr_index] = X[0] * v0 * rho() / R(X);
		}
		break;
	case YL:
		if (vy(X) < 0.0) {
			v0 = vx(X);
			(*this)[lz_index] = (-X[1] * (v0 - X[1] * Binary::Omega)) * rho();
			(*this)[sr_index] = X[0] * v0 * rho() / R(X);
		}
		break;
	case ZU:
		if (sz() > 0.0) {
			(*this)[sz_index] = 0.0;
		}
		break;
	case ZL:
		if (sz() < 0.0) {
			(*this)[sz_index] = 0.0;
		}
		break;
	}
	(*this)[et_index] += 0.5 * (sr() * sr() + st(X) * st(X) + sz() * sz()) / rho();
}

void State::adjust_rho_fluxes() {
	Real sum;
	/*	if (sgn(rho()) != sgn(rho_accretor())) {
	 (*this)[d_acc_index] = 0.0;
	 }
	 if (sgn(rho()) != sgn(rho_donor())) {
	 (*this)[d_don_index] = 0.0;
	 }
	 sum = rho_accretor() + rho_donor();
	 if (sum == 0.0) {
	 (*this)[d_acc_index] = (*this)[d_don_index] = 0.5 * rho();
	 } else {
	 (*this)[d_acc_index] *= rho() / sum;
	 (*this)[d_don_index] *= rho() / sum;
	 }
	 sum = 0.0;
	 if (rho() * specie_rho(0) <= 0.0 && rho() * specie_rho(1) <= 0.0 && rho() * specie_rho(2) <= 0.0) {
	 printf("%e %e %e %e\n", rho(), specie_rho(0), specie_rho(1), specie_rho(2));
	 }
	 for (int i = 0; i < NSPECIES; i++) {
	 if (rho() * specie_rho(i) < 0.0) {
	 (*this)[specie_index + i] = 0.0;
	 }
	 sum += specie_rho(i);
	 }
	 if (sum == 0.0) {
	 printf("!!!!!!!!!!\n");
	 for (int i = 0; i < NSPECIES; i++) {
	 (*this)[specie_index + i] = rho() / Real(NSPECIES);
	 }
	 } else {
	 for (int i = 0; i < NSPECIES; i++) {
	 (*this)[specie_index + i] *= rho() / sum;
	 }
	 }*/
}

void State::floor(const _3Vec& X) {
	if (rho() < rho_floor) {
		(*this)[d_index] = rho_floor;
		(*this)[d_acc_index] = rho_floor * 0.5;
		(*this)[d_don_index] = rho_floor * 0.5;
		(*this)[specie_index + 0] = rho_floor;
		for (int i = 1; i < NSPECIES; i++) {
			(*this)[specie_index + i] = 0.0;
		}
	}
	/*	Real sum;
	 for (int i = 0; i < NSPECIES; i++) {
	 if (specie_rho(i) < 0.0) {
	 (*this)[specie_index + i] = 0.0;
	 }
	 sum += specie_rho(i);
	 }
	 if (sum != 0.0) {
	 for (int i = 0; i < NSPECIES; i++) {
	 (*this)[specie_index + i] *= rho() / sum;
	 }
	 }*/
#ifndef SCF_CODE
	Real u, d;
	d = rho();
	u = ei(X);
	eos_t eos;
	Real dT;
	eos.abar = A();
	eos.zbar = Z();
	eos.rho = d;
	eos.e = u / d;
	eos.T = (*this)[T_index];
	helmholtz_compute_T(&eos);
	(*this)[T_index] = eos.T;
	if (u > 0.1 * et()) {
		(*this)[tau_index] = u;
	}
#else
	Real e = 0.0;
	eos_t eos;
	eos.abar = A();
	eos.zbar = Z();
	eos.T = (1.0e+5 * specie_rho(0) + 2.0e+3 * (rho() - specie_rho(0))) / rho();
	eos.rho = rho();
	helmholtz_eos(&eos);
	e = eos.rho * eos.e;
	(*this)[et_index] = e;
	(*this)[tau_index] = e;
	(*this)[T_index] = eos.T;
	(*this)[sr_index] = 0.0;
	(*this)[lz_index] = Binary::Omega * R(X) * R(X) * rho();
	(*this)[sz_index] = 0.0;
#endif
}

void State::reflect_on_z() {
	(*this)[sz_index] = -(*this)[sz_index];
}

void State::to_con(const _3Vec& x) {
	(*this)[sz_index] *= rho();
	(*this)[pot_index] *= rho();
	(*this)[sr_index] *= rho();
	(*this)[lz_index] += Binary::Omega * R(x);
	(*this)[lz_index] *= rho() * R(x);
	(*this)[et_index] += ek(x);
#ifndef ZTWD
	(*this)[tau_index] = pow((*this)[tau_index], 1.0 / gamma);
#else
	double sum;
	sum = 0.0;
	for (int i = 0; i < NSPECIES; i++) {
		sum += (*this)[specie_index + i];
	}
	for (int i = 0; i < NSPECIES; i++) {
		(*this)[specie_index + i] /= sum;
	}
	sum = 0.0;
	sum += (*this)[d_don_index];
	sum += (*this)[d_acc_index];
	(*this)[d_don_index] /= sum;
	(*this)[d_acc_index] /= sum;
	for (int i = 0; i < NSPECIES; i++) {
		(*this)[specie_index + i] *= rho();
	}
	(*this)[d_don_index] *= rho();
	(*this)[d_acc_index] *= rho();
#endif
}

void State::to_prim(const _3Vec& x) {
#ifndef ZTWD
	(*this)[tau_index] = pow((*this)[tau_index], gamma);
#else
	for (int i = 0; i < NSPECIES; i++) {
		(*this)[specie_index + i] /= rho();
	}
	(*this)[d_don_index] /= rho();
	(*this)[d_acc_index] /= rho();

#endif
	(*this)[et_index] -= ek(x);
	(*this)[sr_index] /= rho();
	(*this)[lz_index] /= rho() * R(x);
	(*this)[lz_index] -= Binary::Omega * R(x);
	(*this)[sz_index] /= rho();
	(*this)[pot_index] /= rho();
}




#endif


#else



#endif
