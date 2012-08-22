#include "../defs.h"
#ifdef SINGLE

//#define POLYTROPIC

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "oct_node/oct_face.h"
#include "state.h"
const Real State::gamma = 5.0 / 3.0;
const Real State::rho_floor = 1.0e-20;
const Real State::ei_floor = 1.0e-20;
const int State::d_index = 0;
const int State::sx_index = 1;
const int State::sy_index = 2;
const int State::sz_index = 3;
const int State::et_index = 4;
const int State::tau_index = 5;
const int State::pot_index = 6;

Vector<Real, STATE_NF> State::scalar_x_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
#ifdef USE_LZ
	c[sx_index] = X[0] / sqrt(X[0] * X[0] + X[1] * X[1]);
#endif
	return c;
}

Vector<Real, STATE_NF> State::scalar_y_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
#ifdef USE_LZ
	c[sx_index] = X[1] / sqrt(X[0] * X[0] + X[1] * X[1]);
#endif
	return c;
}

Vector<Real, STATE_NF> State::scalar_z_coeff(const _3Vec& X) {
	Vector<Real, STATE_NF> c = 0.0;
	return c;
}

Vector<Real, STATE_NF> State::scalar_x_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
#ifdef USE_LZ
	flux[sx_index] = pg(X);
#endif
	return flux;
}

Vector<Real, STATE_NF> State::scalar_y_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
#ifdef USE_LZ
	flux[sx_index] = pg(X);
#endif
	return flux;
}

Vector<Real, STATE_NF> State::scalar_z_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	return flux;
}

Vector<Real, STATE_NF> State::source(const _3Vec& X) const {
#ifdef USE_LZ
	State s = Vector<Real, STATE_NF> (0.0);
	const Real R = sqrt(X[0] * X[0] + X[1] * X[1]);
	const Real theta = atan2(X[1], X[0]);
	Real d, R2;
	s[sx_index] = pow(sy(), 2) * pow(R, -3) / rho();
	//	s[sx_index] += fR * cos(theta - ftheta) * rho();
	//	s[sy_index] -= fR * sin(theta - ftheta) * rho() * R;
	return s;
#else
	State s = Vector<Real, STATE_NF> (0.0);
	return s;
#endif
}

const char* State::field_name(int i) {
	static char dstr[3];
	dstr[2] = '\0';
	dstr[0] = 'd';
	assert(i >= 0);assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d";
	case sz_index:
		return "sz";
	case et_index:
		return "et";
#ifdef USE_LZ
	case sx_index:
		return "sr";
	case sy_index:
		return "lz";
#else
		case sx_index:
		return "sx";
		case sy_index:
		return "sy";
#endif
	case pot_index:
		return "pot";
	}
	return "tau";
}

void State::enforce_outflow(const OctFace& f, const _3Vec& X) {
	switch (f) {
#ifdef USE_LZ
	case XU:
		if (vx(X) > 0.0) {
			(*this)[sy_index] = X[0] * vy(X) * rho();
			(*this)[sx_index] = X[1] * vy(X) * rho() / sqrt(X[0] * X[0] + X[1] * X[1]);
		}
		break;
	case XL:
		if (vx(X) < 0.0) {
			(*this)[sy_index] = X[0] * vy(X) * rho();
			(*this)[sx_index] = X[1] * vy(X) * rho() / sqrt(X[0] * X[0] + X[1] * X[1]);
		}
		break;
	case YU:
		if (vy(X) > 0.0) {
			(*this)[sy_index] = -X[1] * vx(X) * rho();
			(*this)[sx_index] = X[0] * vx(X) * rho() / sqrt(X[0] * X[0] + X[1] * X[1]);
		}
		break;
	case YL:
		if (vy(X) < 0.0) {
			(*this)[sy_index] = -X[1] * vx(X) * rho();
			(*this)[sx_index] = X[0] * vx(X) * rho() / sqrt(X[0] * X[0] + X[1] * X[1]);
		}
		break;
#else
		case XU:
		if (vx() > 0.0) {
			set_et(et() - 0.5 * sx() * sx() / rho());
			set_sx(0.0);
		}
		break;
		case XL:
		if (vx() < 0.0) {
			set_et(et() - 0.5 * sx() * sx() / rho());
			set_sx(0.0);
		}
		break;
		case YU:
		if (vy() > 0.0) {
			set_et(et() - 0.5 * sy() * sy() / rho());
			set_sy(0.0);
		}
		break;
		case YL:
		if (vy() < 0.0) {
			set_et(et() - 0.5 * sy() * sy() / rho());
			set_sy(0.0);
		}
		break;
#endif
	case ZU:
		if (sz() > 0.0) {
			set_et(et() - 0.5 * sz() * sz() / rho());
			set_sz(0.0);
		}
		break;
	case ZL:
		if (sz() < 0.0) {
			set_et(et() - 0.5 * sz() * sz() / rho());
			set_sz(0.0);
		}
		break;
	}
}

void State::floor(const _3Vec& X) {
	(*this)[d_index] = max((*this)[d_index], rho_floor);
#ifndef POLYTROPIC
	const Real ei0 = et() - ek(X);
	if (ei0 > 0.1 * et()) {
		(*this)[tau_index] = pow(max(ei0, ei_floor), 1.0 / gamma);
	}
#endif
}

void State::reflect_on_z() {
	set_sz(-sz());
}

Real State::poisson_source() const {
	Real a = 4.0 * M_PI * G_CON * rho();
	return a;
}

Real State::refine_value() const {
	return rho();
}

State::State() :
	Vector<Real, STATE_NF> () {
}

State::State(const Vector<Real, STATE_NF>& v) :
	Vector<Real, STATE_NF> (v) {
	return;
}

Real State::rho() const {
	return (*this)[d_index];
}

Real State::sx() const {
	return (*this)[sx_index];
}

Real State::sy() const {
	return (*this)[sy_index];
}

Real State::sz() const {
	return (*this)[sz_index];
}

Real State::et() const {
	return (*this)[et_index];
}

Real State::vx(const _3Vec& X) const {
#ifndef USE_LZ
	return sx() / rho();
#else
	const Real Rinv = 1.0 / sqrt(X[0] * X[0] + X[1] * X[1]);
	return (-X[1] * sy() * Rinv + X[0] * sx()) * Rinv / rho();
#endif
}

Real State::vy(const _3Vec& X) const {
#ifndef USE_LZ
	return sy() / rho();
#else
	const Real Rinv = 1.0 / sqrt(X[0] * X[0] + X[1] * X[1]);
	return (X[0] * sy() * Rinv + X[1] * sx()) * Rinv / rho();
#endif
}

Real State::vz(const _3Vec& X) const {
	return sz() / rho();
}

Real State::ek(const _3Vec& X) const {
#ifdef USE_LZ
	return 0.5 * (sx() * sx() + sy() * sy() / (X[0] * X[0] + X[1] * X[1]) + sz() * sz()) / rho();
#else
	return 0.5 * (sx() * sx() + sy() * sy()+ sz() * sz()) / rho();
#endif
}

Real State::ei(const _3Vec& X) const {
#ifdef POLYTROPIC
	Real ei0;
#else POLYTROPIC
	Real ei0 = et() - ek(X);
	if (ei0 <= 0.001 * et()) {
#endif
		ei0 = pow((*this)[tau_index], gamma);
#ifndef POLYTROPIC
	}
#endif
	return max(ei0, ei_floor);
}

Real State::pg(const _3Vec& X) const {
	return (gamma - 1.0) * ei(X);
}

Real State::cs(const _3Vec& X) const {
	assert(rho() > 0.0);assert(pg(X) >= 0.0);
	return sqrt(gamma * pg(X) / rho());
}

Real State::max_abs_x_eigen(const _3Vec& X, const _3Vec& V) const {
	return fabs(vx(X) - V[0]) + cs(X);
}

Real State::max_abs_y_eigen(const _3Vec& X, const _3Vec& V) const {
	return fabs(vy(X) - V[1]) + cs(X);
}

Real State::max_abs_z_eigen(const _3Vec& X, const _3Vec& V) const {
	return fabs(vz(X) - V[2]) + cs(X);
}

Vector<Real, STATE_NF> State::vector_x_flux(const _3Vec& X, const _3Vec& V) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, v0, p;
	v = vx(X);
	v0 = v - V[0];
	p = pg(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v0;
	}
	flux[et_index] += v * p;
#ifndef USE_LZ
	flux[sx_index] += p;
#else
	flux[sy_index] -= X[1] * p;
#endif
	return flux;
}

Vector<Real, STATE_NF> State::vector_y_flux(const _3Vec& X, const _3Vec& V) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, v0, p;
	v = vy(X);
	v0 = v - V[1];
	p = pg(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v0;
	}
	flux[et_index] += v * p;
#ifndef USE_LZ
	flux[sy_index] += p;
#else
	flux[sy_index] += X[0] * p;
#endif
	return flux;
}

Vector<Real, STATE_NF> State::vector_z_flux(const _3Vec& X, const _3Vec& V) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, v0, p;
	v = vz(X);
	v0 = v - V[2];
	p = pg(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v0;
	}
	flux[sz_index] += p;
	flux[et_index] += v * p;
	return flux;
}

void State::set_rho(Real a, int frac) {
	(*this)[d_index + frac] = a;
}

void State::set_sx(Real a) {
	(*this)[sx_index] = a;
}

void State::set_sy(Real a) {
	(*this)[sy_index] = a;
}

void State::set_sz(Real a) {
	(*this)[sz_index] = a;
}

void State::set_et(Real a) {
	(*this)[et_index] = a;
}

void State::set_tau(Real a) {
	(*this)[tau_index] = a;
}

void State::to_prim(const _3Vec& X) {
	(*this)[et_index] -= ek(X);
	(*this)[sx_index] /= rho();
	(*this)[sy_index] /= rho();
	(*this)[sz_index] /= rho();
}

void State::to_con(const _3Vec& X) {
	(*this)[sx_index] *= rho();
	(*this)[sy_index] *= rho();
	(*this)[sz_index] *= rho();
	(*this)[et_index] += ek(X);
}

#endif
