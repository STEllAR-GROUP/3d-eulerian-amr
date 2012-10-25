#include "../defs.h"

#ifdef BINARY

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "oct_node/oct_face.h"
#include "state.h"
#include "binary.h"

#ifndef ZTWD

const Real State::gamma = EULER_GAMMA;
const Real State::ei_floor = 1.0e-20;
const int State::d_index = 0;
const int State::sr_index = 0 + NRHO;
const int State::lz_index = 1 + NRHO;
const int State::sz_index = 2 + NRHO;
const int State::pot_index = 3 + NRHO;
const int State::et_index = 4 + NRHO;
const int State::tau_index = 5 + NRHO;
Real State::ftheta, State::fR;
Real State::rho_floor = 1.0e-10;

void State::initialize(int acc, int don) {
	switch (acc) {
	case 1:
		abar[0] = 4.0;
		abar[2] = 4.0;
		break;
	case 2:
		abar[0] = 4.0;
		abar[2] = 14.118;
		break;
	case 3:
		abar[0] = 14.118;
		abar[2] = 14.118;
		break;
	case 4:
		abar[0] = 17.518;
		abar[2] = 17.518;
		break;
	}
	switch (don) {
	case 1:
		abar[1] = 4.0;
		abar[3] = 4.0;
		break;
	case 2:
		abar[1] = 4.0;
		abar[3] = 14.118;
		break;
	case 3:
		abar[1] = 14.118;
		abar[3] = 14.118;
		break;
	case 4:
		abar[1] = 17.518;
		abar[3] = 17.518;
		break;
	}
	for (int i = 0; i < 4; i++) {
		zbar[i] = abar[i] / 2.0;
	}
}

bool State::low_order_variable(int i) {
	return (pot_index == i);
}

bool State::smooth_variable(int i) {
	return (pot_index == i);
}

const char* State::field_name(int i) {
	static char dstr[3];
	dstr[2] = '\0';
	dstr[0] = 'd';
	assert(i >= 0);assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d1";
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
	default:
		dstr[1] = '1' + i;
		return dstr;
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
	return sqrt(gamma * pg(X) / rho());
}

Real State::enthalpy() const {
	Real h;
	h = Binary::K1 * (Binary::n + 1.0) * pow(rho(0), 1.0 / Binary::n);
	h += Binary::K2 * (Binary::n + 1.0) * pow(rho(1), 1.0 / Binary::n);
	return h;
}

Real State::rho(int i) const {
	if (i < 0) {
		return rho();
	} else {
		return max(rho_floor / Real(NRHO), (*this)[d_index + i]);
	}
}

Real State::ei(const _3Vec& X) const {
	Real e;
	e = et() - ek(X);
	if (e <= 0.001 * et()) {
		e = pow((*this)[tau_index], gamma);
	}
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
	return (gamma - 1.0) * ei(X);
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
	Real rho_tot = 0.0;
	for (int i = 0; i < NRHO; i++) {
		rho_tot += rho(i);
	}
	return rho_tot;
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

Vector<Real, STATE_NF> State::gravity_source(const _3Vec& g, Vector<Real, STATE_NF>& D, const _3Vec& X) const {
	State s = Vector<Real, STATE_NF>(0.0);
	const Real d = rho();
	s[sr_index] += d * (X[0] * g[0] + X[1] * g[1]) / R(X);
	s[lz_index] += d * (X[0] * g[1] - X[1] * g[0]);
	s[sz_index] += d * g[2];
	s[et_index] += D[pot_index];
	for (int i = 0; i < NRHO; i++) {
		s[et_index] -= D[d_index + i] * pot() / d;
	}
	s[et_index] += rho() * (dvx(X) * g[0] + dvy(X) * g[1]);
	return s;
}

Vector<Real, STATE_NF> State::scf_source(Real x, Real y, int lobe) const {
	Vector<Real, STATE_NF> s = 0.0;
	Real phi, rho_np1, h;
	const Real O2 = Binary::Omega * Binary::Omega;
	phi = pot() / rho(); // - 0.5 * (x * x + y * y) * O2;
	for (int i = 0; i < NRHO; i++) {
		s[d_index + i] += (-rho(i));
	}
	if (lobe == 2) {
		if (phi < Binary::phi0_2) {
			rho_np1 = pow((Binary::phi0_2 - phi) / (Binary::K2 * (Binary::n + 1.0)), Binary::n);
			s[d_index + 1] += rho_np1;
		}
	} else if (lobe == 1) {
		if (phi < Binary::phi0_1) {
			rho_np1 = pow((Binary::phi0_1 - phi) / (Binary::K1 * (Binary::n + 1.0)), Binary::n);
			s[d_index + 0] += rho_np1;
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
	switch (f) {
	case XU:
		if (vx(X) > 0.0) {
			(*this)[lz_index] = X[0] * vy(X) * rho();
			(*this)[sr_index] = X[1] * vy(X) * rho() / R(X);
		}
		break;
	case XL:
		if (vx(X) < 0.0) {
			(*this)[lz_index] = X[0] * vy(X) * rho();
			(*this)[sr_index] = X[1] * vy(X) * rho() / R(X);
		}
		break;
	case YU:
		if (vy(X) > 0.0) {
			(*this)[lz_index] = -X[1] * vx(X) * rho();
			(*this)[sr_index] = X[0] * vx(X) * rho() / R(X);
		}
		break;
	case YL:
		if (vy(X) < 0.0) {
			(*this)[lz_index] = -X[1] * vx(X) * rho();
			(*this)[sr_index] = X[0] * vx(X) * rho() / R(X);
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

void State::floor(const _3Vec& X) {
#ifdef SCF_CODE
	for (int i = 0; i < NRHO; i++) {
		(*this)[d_index + i] = max((*this)[d_index + i], rho_floor / Real(NRHO));
	}

	Real e = 0.0;
	e += Binary::K1 * pow((*this)[d_index + 0], 1.0 + 1.0 / Binary::n) / (gamma - 1.0);
	e += Binary::K2 * pow((*this)[d_index + 1], 1.0 + 1.0 / Binary::n) / (gamma - 1.0);
	(*this)[et_index] = e;
	(*this)[tau_index] = pow(e, 1.0 / gamma);
	(*this)[sr_index] = 0.0;
	(*this)[lz_index] = Binary::Omega * R(X) * R(X) * rho();
	(*this)[sz_index] = 0.0;
#else
	for (int i = 0; i < NRHO; i++) {
		(*this)[d_index + i] = max((*this)[d_index + i], rho_floor / Real(NRHO));
	}
	Real e;
	e = et() - ek(X);
	if (e > 0.1 * et()) {
		(*this)[tau_index] = pow(e, 1.0 / gamma);
	}
	(*this)[tau_index] = max(pow(ei_floor, 1.0 / gamma), (*this)[tau_index]);
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
}

void State::to_prim(const _3Vec& x) {
	(*this)[et_index] -= ek(x);
	(*this)[sr_index] /= rho();
	(*this)[lz_index] /= rho() * R(x);
	(*this)[lz_index] -= Binary::Omega * R(x);
	(*this)[sz_index] /= rho();
	(*this)[pot_index] /= rho();
}

#endif
#endif
