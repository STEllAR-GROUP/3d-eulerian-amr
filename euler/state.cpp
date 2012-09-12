#include "../defs.h"
#ifdef EULER

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

Vector<Real, STATE_NF> State::source(const _3Vec& X) const {
	State s = Vector<Real, STATE_NF> (0.0);
	return s;
}

const char* State::field_name(int i) {
	static char dstr[3];
	dstr[2] = '\0';
	dstr[0] = 'd';
	assert(i >= 0);
	assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d";
	case sz_index:
		return "sz";
	case et_index:
		return "et";
	case sx_index:
		return "sx";
	case sy_index:
		return "sy";
	}
	return "tau";
}

void State::enforce_outflow(const _3Vec& X, const OctFace& f) {
	switch (f) {
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
	const Real ei0 = et() - ek();
	if (ei0 > 0.1 * et()) {
		(*this)[tau_index] = pow(max(ei0, ei_floor), 1.0 / gamma);
	}
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

Real State::vx() const {
	return sx() / rho();
}

Real State::vy() const {
	return sy() / rho();
}

Real State::vz() const {
	return sz() / rho();
}

Real State::ek() const {
	return 0.5 * (sx() * sx() + sy() * sy() + sz() * sz()) / rho();
}

Real State::ei() const {
	return max(et() - ek(), ei_floor);
}

Real State::pg() const {
	const Real ei0 = et() - ek();
	if (ei0 < 0.001 * et()) {
		return (gamma - 1.0) * pow((*this)[tau_index], gamma);
	} else {
		return (gamma - 1.0) * ei0;
	}
}

Real State::cs() const {
	assert(rho() > 0.0);
	assert(pg() >= 0.0);
	return sqrt(gamma * pg() / rho());
}

Real State::max_abs_x_eigen(const _3Vec& X) const {
	return fabs(vx()) + cs();
}

Real State::max_abs_y_eigen(const _3Vec& X) const {
	return fabs(vy()) + cs();
}

Real State::max_abs_z_eigen(const _3Vec& X) const {
	return fabs(vz()) + cs();
}

Vector<Real, STATE_NF> State::x_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p;
	v = vx();
	p = pg();
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
	}
	flux[sx_index] += p;
	flux[et_index] += v * p;
	return flux;
}

Vector<Real, STATE_NF> State::y_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p;
	v = vy();
	p = pg();
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
	}
	flux[sy_index] += p;
	flux[et_index] += v * p;
	return flux;
}

Vector<Real, STATE_NF> State::z_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p;
	v = vz();
	p = pg();
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
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
	(*this)[et_index] -= ek();
	(*this)[sx_index] /= rho();
	(*this)[sy_index] /= rho();
	(*this)[sz_index] /= rho();
}

void State::to_con(const _3Vec& X) {
	(*this)[sx_index] *= rho();
	(*this)[sy_index] *= rho();
	(*this)[sz_index] *= rho();
	(*this)[et_index] += ek();
}

#endif
