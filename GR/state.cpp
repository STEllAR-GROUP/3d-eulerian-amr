#include "../defs.h"
#ifdef GR

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "oct_node/oct_face.h"
#include "state.h"

Real State::d_floor = 1.0e-20;
const Real State::fgamma = 5.0 / 3.0;
const int State::symtable[3][3] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };
const Real State::delta[3][3] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
const int State::d_index = 0;
const int State::tau_index = 1;
const int State::s_index = 2;
const int State::Theta_index = 5;
const int State::A_index = 6;
const int State::Z_index = 9;
const int State::K_index = 12;
const int State::D_index[3] = { 18, 24, 30 };
const int State::rho_index = 0;
const int State::v_index = 2;
const int State::eta_index = 1;

bool State::fluid_on = false;

Vector<Real, STATE_NF> State::source(const _3Vec& X) const {
	State src = Vector<Real, STATE_NF>(0.0);
	src.Theta() += -(tau() + d());
	for (int i = 0; i < 3; i++) {
		src.Z(i) += -S(i);
		for (int j = i; j < 3; j++) {
			src.K(i, j) += -S(i, j) + delta[i][j] * (TrS() - (tau() + d())) * 0.5;
		}
	}
	if (fluid_on) {
		for (int k = 0; k < 3; k++) {
			for (int i = 0; i < 4; i++) {
				for (int m = 0; m < 4; m++) {
					src.s(k) += -Christoffel(k, i, m) * S(i, m);
				}
			}
		}
		for (int i = 0; i < 4; i++) {
			for (int m = 0; m < 4; m++) {
				src.tau() += -Christoffel(3, i, m) * S(i, m);
			}
		}
	}
	return src;
}

Vector<Real, STATE_NF> State::flux(int k) const {
	State flux = Vector<Real, STATE_NF>(0.0);
	if (fluid_on) {
		flux.d() += d() * v(k);
		flux.tau() += s(k) - flux.d();
		for (int i = 0; i < 3; i++) {
			flux.s(i) += s(i) * v(k) + delta[i][k] * pg();
		}
	}
	flux.Theta() += D(k) - E(k) - Z(k);
	flux.A(k) += TrK() - 2.0 * Theta();
	for (int i = 0; i < 3; i++) {
		flux.Z(i) += delta[i][k] * (TrK() - Theta()) - K(k, i);
		for (int j = i; j < 3; j++) {
			flux.D(k, i, j) += K(i, j);
			flux.K(i, j) += lambda(k, i, j);
		}
	}
	return flux;
}

const char* State::field_name(int i) {
	const char* names[STATE_NF] = { "d", "tau", "sx", "sy", "sz", "Theta", "A_1", "A_2", "A_3", "Z_1", "Z_2", "Z_3",
			"K_11", "K_12", "K_13", "K_22", "K_23", "K_33", "D1_11", "D1_12", "D1_13", "D1_22", "D1_23", "D1_33",
			"D2_11", "D2_12", "D2_13", "D2_22", "D2_23", "D2_33", "D3_11", "D3_12", "D3_13", "D3_22", "D3_23", "D3_33" };
	return names[i];
}

void State::enforce_outflow(const _3Vec& X, const OctFace& f) {
}

void State::compute_prims() {
	Real pbar, p, f, err, v[3], W, eta, rho, dfdpbar, v2, cs2, h;
	pbar = 0.0;
	do {
		for (int i = 0; i < 3; i++) {
			v2 += pow(s(i) / (tau() + d() + pbar), 2.0);
		}
		if (v2 > 1.0) {
			printf("Unphysical momentum\n");
			abort();
		}
		W = 1.0 / sqrt(1.0 - v2);
		eta = max((tau() + d() * (1.0 - W) + pbar * (1.0 - W * W)) / (d() * W), 0.0);
		rho = d() / W;
		p = (fgamma - 1.0) * rho * eta;
		f = p - pbar;
		cs2 = fgamma * (fgamma - 1.0) * eta / (fgamma * eta + 1.0);
		dfdpbar = v2 * cs2 - 1.0;
		pbar -= f / dfdpbar;
		err = fabs(f);
	} while (err > 1.0e-10);
	tau() = max(max(tau(), sqrt(s(0) * s(0) + s(1) * s(1) + s(2) * s(2)) - d() - p), 0.0);
	h = 1.0 + p * W / d() + eta;
	prims[rho_index] = d() / W;
	assert( prims[rho_index]!= 0.0);
	prims[eta_index] = eta;
	for (int i = 0; i < 3; i++) {
		prims[v_index + i] = s(i) / rho / W / W / h;
	}
#ifndef NDEBUG
	prims_called = true;
#endif
}

void State::floor(const _3Vec& X) {
	d() = max(d(), d_floor);
	compute_prims();
}

State::State() :
		Vector<Real, STATE_NF>() {
	prims = 0.0;
#ifndef NDEBUG
	prims_called = false;
#endif
}

State::State(const Vector<Real, STATE_NF>& v) :
		Vector<Real, STATE_NF>(v) {
	return;
}

Real State::rho() const {
	return prims[rho_index];
}

Real State::max_abs_x_eigen(const _3Vec& X) const {
	return 1.0;
}

Real State::max_abs_y_eigen(const _3Vec& X) const {
	return 1.0;
}

Real State::max_abs_z_eigen(const _3Vec& X) const {
	return 1.0;
}

void State::reflect_on_z() {
}

Vector<Real, STATE_NF> State::x_flux(const _3Vec& X) const {
	return flux(0);
}

Vector<Real, STATE_NF> State::y_flux(const _3Vec& X) const {
	return flux(1);
}

Vector<Real, STATE_NF> State::z_flux(const _3Vec& X) const {
	return flux(2);
}

void State::to_prim(const _3Vec& X) {
//	compute_prims();
	for (int i = 0; i < PRIM_NF; i++) {
		(*this)[i] = prims[i];
	}
//	}
}

void State::to_con(const _3Vec& X) {

//	if (fluid_on) {
	for (int i = 0; i < PRIM_NF; i++) {
		prims[i] = (*this)[i];
	}
#ifndef NDEBUG
	prims_called = true;
#endif
	for (int i = 0; i < 3; i++) {
		s(i) = rho() * h() * W() * W() * v(i);
	}
	d() = W() * rho();
	tau() = rho() * h() * W() * W() - d() - pg();
//	}
}

#endif
