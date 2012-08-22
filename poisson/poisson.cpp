#include "poisson.h"
#include "indexer_2d_by2.h"
#include "indexer_2d.h"
#include <omp.h>

#define RELAX_RESID 0.1
#define VDOWN_RESID 1.0
#define VUP_RESID 0.01

const int interp_c[64] = { 0, 0, 0, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, -1, 11, 11, -1, -1, 11,
		11, -1, 0, -1, -1, 0, 0, -1, -1, 0, -1, 11, 11, -1, -1, 11, 11, -1, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0,
		-1, -1, 0, 0, 0, 0, 0 };
const Real interp_c0 = 64.0;

void Poisson::vcycle_coarse_correction(int level) {
	if (get_level() == level) {
		if (get_level() > 0) {
			ChildIndex c = get_index();
			const Poisson* p = static_cast<const Poisson*> (get_parent());
			const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
			int i, j, k, i0, j0, k0;
			Real a, dx, dy, dz;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(a,i,j,k,i0,j0,k0,dx,dy,dz)
			for (int index = 0; index <= indexer.max_index(); index++) {
				k = indexer.y(index);
				j = indexer.x(index);
				k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
				j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
				for (i = 1, i0 = 1 + c.get_x() * (PNX / 2 - 1); i < PNX - 1; i += 2, i0++) {
					a = p->dphi(i0, j0, k0) - dphi.oct_avg(i, j, k);
					dphi(i, j, k) += a;
					dphi(i + 1, j, k) += a;
					dphi(i, j + 1, k) += a;
					dphi(i + 1, j + 1, k) += a;
					dphi(i, j, k + 1) += a;
					dphi(i + 1, j, k + 1) += a;
					dphi(i, j + 1, k + 1) += a;
					dphi(i + 1, j + 1, k + 1) += a;
				}
			}
		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->vcycle_coarse_correction(level);
			}
		}
	}
}

void Poisson::vcycle0() {
	const int maxlevel = this->max_level();
	int i;
	phi_from_children();
	Real resid0, resid;
	for (int l = 0; l <= maxlevel; l++) {
		enforce_phi_boundaries(l);
	}
	init_vcycle();
	for (int l = maxlevel; l >= 0; l--) {
		enforce_dphi_boundaries(l);
		vcycle_down(l);
		enforce_dphi_boundaries(l);
		resid0 = vcycle_relax(l);
		do {
			enforce_dphi_boundaries(l);
			resid = vcycle_relax(l);
		} while (resid / resid0 > VDOWN_RESID);
	}
	for (int l = 0; l <= maxlevel; l++) {
		enforce_dphi_boundaries(l);
		vcycle_coarse_correction(l);
		enforce_dphi_boundaries(l);
		resid0 = vcycle_relax(l);
		do {
			enforce_dphi_boundaries(l);
			resid = vcycle_relax(l);
		} while (resid / resid0 > VUP_RESID);
		vcycle_retire_dphi(l);
		enforce_phi_boundaries(l);
	}
}

Real Poisson::vcycle_relax(int level) {
	Array3d<Real, PNX, PNX, PNX> tmp0;
	Real resid = 0.0;
	if (get_level() == level) {
		bool first_pass = true;
		Real resid0;
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int i, j, k;
		do {
			resid = 0.0;
#pragma omp parallel for schedule(OMP_SCHEDULE) collapse(2)
			for (int k = 1; k < PNX - 1; k++) {
				for (int j = 1; j < PNX - 1; j++) {
					for (int i = 1; i < PNX - 1; i++) {
						tmp0(i, j, k) = dphi1(i, j, k) + dphi.divergence(i, j, k) * (1.0 / 6.0);
					}
				}
			}
#pragma omp parallel for schedule(OMP_SCHEDULE) collapse(2) reduction(+:resid)
			for (int k = 1; k < PNX - 1; k++) {
				for (int j = 1; j < PNX - 1; j++) {
					for (int i = 1; i < PNX - 1; i++) {
						dphi(i, j, k) += tmp0(i, j, k);
						resid += tmp0(i, j, k) * tmp0(i, j, k);
					}
				}
			}
			if (first_pass) {
				resid0 = resid;
				first_pass = false;
			}
		} while (resid / resid0 > RELAX_RESID);
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				resid += static_cast<Poisson*> (get_child(i))->vcycle_relax(level);
			}
		}
	}
	return resid;
}

void Poisson::inject_from_parent(ChildIndex c) {
	GridNode::inject_from_parent(c);
	const Poisson* p = static_cast<const Poisson*> (get_parent());
	const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
	Real u;
	int k, j, k0, j0, i, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(u,k,j,i,k0,j0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
		j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
		for (i = 1, i0 = 1 + c.get_x() * (PNX / 2 - 1); i < PNX - 1; i += 2, i0++) {
			u = p->phi(i0, j0, k0);
			phi(i + 0, j + 0, k + 0) = u;
			phi(i + 1, j + 0, k + 0) = u;
			phi(i + 0, j + 1, k + 0) = u;
			phi(i + 1, j + 1, k + 0) = u;
			phi(i + 0, j + 0, k + 1) = u;
			phi(i + 1, j + 0, k + 1) = u;
			phi(i + 0, j + 1, k + 1) = u;
			phi(i + 1, j + 1, k + 1) = u;
			phi0(i + 0, j + 0, k + 0) = u;
			phi0(i + 1, j + 0, k + 0) = u;
			phi0(i + 0, j + 1, k + 0) = u;
			phi0(i + 1, j + 1, k + 0) = u;
			phi0(i + 0, j + 0, k + 1) = u;
			phi0(i + 1, j + 0, k + 1) = u;
			phi0(i + 0, j + 1, k + 1) = u;
			phi0(i + 1, j + 1, k + 1) = u;
		}
	}
}

void Poisson::store_phi() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			phi0(i, j, k) = phi(i, j, k);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->store_phi();
		}
	}
}

void Poisson::restore_phi() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			phi(i, j, k) = phi0(i, j, k);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->restore_phi();
		}
	}
}

void Poisson::restore_half_phi() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i;
	Real r;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,r)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			r = phi(i, j, k);
			phi(i, j, k) = (phi(i, j, k) + phi0(i, j, k)) * 0.5;
			phi0(i, j, k) = r;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->restore_half_phi();
		}
	}
}

void Poisson::compute_forces() {
	const Indexer2d indexer(1, PNX - 1, 1, PNX - 1);
	const Real h = get_dx();
	int j, k, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX; i++) {
			fx(i, j, k) = -(phi(i, j, k) - phi(i - 1, j, k)) / h;
			fy(i, j, k) = -(phi(i, j, k) - phi(i, j - 1, k)) / h;
			fz(i, j, k) = -(phi(i, j, k) - phi(i, j, k - 1)) / h;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->compute_forces();
		}
	}
	if (get_level() == 0) {
		adjust_x_force();
		adjust_y_force();
		adjust_z_force();
	}
}

void Poisson::adjust_x_dphi() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(0, cj, ck);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const Poisson*> (get_sibling(XL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_x(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(XU)->get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v, d;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,d,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->phi(i, j, k) - child->phi(i - 1, j, k);
						v += child->phi(i, j + 1, k) - child->phi(i - 1, j + 1, k);
						v += child->phi(i, j, k + 1) - child->phi(i - 1, j, k + 1);
						v += child->phi(i, j + 1, k + 1) - child->phi(i - 1, j + 1, k + 1);
						v *= 0.5;
						d = (phi(i0, j0, k0) - phi(i0 - 1, j0, k0) - v) * (1.0 / 6.0);
						if (!poisson_zone_is_refined(i0, j0, k0)) {
							dphi(i0, j0, k0) += d;
						}
						if (!poisson_zone_is_refined(i0 - 1, j0, k0)) {
							dphi(i0 - 1, j0, k0) -= d;
						}

					}
				}
			}
		}
	}
}

void Poisson::adjust_y_dphi() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, 0, ck);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const Poisson*> (get_sibling(YL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_y(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(YU)->get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v, d;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,d,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->phi(j, i, k) - child->phi(j, i - 1, k);
						v += child->phi(j + 1, i, k) - child->phi(j + 1, i - 1, k);
						v += child->phi(j, i, k + 1) - child->phi(j, i - 1, k + 1);
						v += child->phi(j + 1, i, k + 1) - child->phi(j + 1, i - 1, k + 1);
						v *= 0.5;
						d = (phi(j0, i0, k0) - phi(j0, i0 - 1, k0) - v) * (1.0 / 6.0);
						if (!poisson_zone_is_refined(j0, i0, k0)) {
							dphi(j0, i0, k0) += d;
						}
						if (!poisson_zone_is_refined(j0, i0 - 1, k0)) {
							dphi(j0, i0 - 1, k0) -= d;
						}
					}
				}
			}
		}
	}
}

void Poisson::adjust_z_dphi() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, ck, 0);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const Poisson*> (get_sibling(ZL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_z(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(ZU)->get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v, d;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,d,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->phi(j, k, i) - child->phi(j, k, i - 1);
						v += child->phi(j + 1, k, i) - child->phi(j + 1, k, i - 1);
						v += child->phi(j, k + 1, i) - child->phi(j, k + 1, i - 1);
						v += child->phi(j + 1, k + 1, i) - child->phi(j + 1, k + 1, i - 1);
						v *= 0.5;
						d = (phi(j0, k0, i0) - phi(j0, k0, i0 - 1) - v) * (1.0 / 6.0);
						if (!poisson_zone_is_refined(j0, k0, i0)) {
							dphi(j0, k0, i0) += d;
						}
						if (!poisson_zone_is_refined(j0, k0, i0 - 1)) {
							dphi(j0, k0, i0 - 1) -= d;
						}
					}
				}
			}
		}
	}
}

void Poisson::adjust_x_force() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->adjust_x_force();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(0, cj, ck);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const Poisson*> (get_sibling(XL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_x(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(XU)->get_child(ci));
					ci.set_x(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->fx(i, j, k);
						v += child->fx(i, j + 1, k);
						v += child->fx(i, j, k + 1);
						v += child->fx(i, j + 1, k + 1);
						fx(i0, j0, k0) = v * 0.25;
					}
				}
			}
		}
	}
}

void Poisson::adjust_y_force() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->adjust_y_force();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, 0, ck);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const Poisson*> (get_sibling(YL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_y(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(YU)->get_child(ci));
					ci.set_y(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->fy(j, i, k);
						v += child->fy(j + 1, i, k);
						v += child->fy(j, i, k + 1);
						v += child->fy(j + 1, i, k + 1);
						fy(j0, i0, k0) = v * 0.25;
					}
				}
			}
		}
	}
}

void Poisson::adjust_z_force() {
	ChildIndex ci;
	const Poisson* child_l;
	const Poisson* child_r;
	const Poisson* child;
	int i, i0;

	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->adjust_z_force();
		}
	}

	for (int l = 0; l < 3; l++) {
		for (int ck = 0; ck < 2; ck++) {
			for (int cj = 0; cj < 2; cj++) {
				ci.set_index(cj, ck, 0);
				if (l == 0) {
					child_r = static_cast<const Poisson*> (get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const Poisson*> (get_sibling(ZL)->get_child(ci));
					i0 = 1;
				} else if (l == 1) {
					child_l = static_cast<const Poisson*> (get_child(ci));
					ci.set_z(1);
					child_r = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX / 2;
				} else if (l == 2) {
					child_r = static_cast<const Poisson*> (get_sibling(ZU)->get_child(ci));
					ci.set_z(1);
					child_l = static_cast<const Poisson*> (get_child(ci));
					i0 = PNX - 1;
				}
				if ((child_r == NULL) && (child_l != NULL)) {
					child = child_l;
					i = PNX - 1;
				} else if ((child_l == NULL) && (child_r != NULL)) {
					child = child_r;
					i = 1;
				} else {
					child = NULL;
				}
				if (child != NULL) {
					const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
					Real v;
					int j, k, j0, k0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,j,k,j0,k0)
					for (int index = 0; index <= indexer.max_index(); index++) {
						k = indexer.y(index);
						j = indexer.x(index);
						k0 = (1 + k) / 2 + ck * (PNX / 2 - 1);
						j0 = (1 + j) / 2 + cj * (PNX / 2 - 1);
						v = +child->fz(j, k, i);
						v += child->fz(j + 1, k, i);
						v += child->fz(j, k + 1, i);
						v += child->fz(j + 1, k + 1, i);
						fz(j0, k0, i0) = v * 0.25;
					}
				}
			}
		}
	}
}

Real Poisson::solution_error() {
	Real sum;
	const Real h3 = pow(get_dx(), 3);
	const Real h2 = pow(get_dx(), 2);
	const Real h1 = get_dx();
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	sum = 0.0;
	int i, j, k;
	Real df;
	const Real h1inv = 1.0 / h1;
#pragma omp parallel for schedule(OMP_SCHEDULE) reduction(+:sum) private(i,j,k,df)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			if (!poisson_zone_is_refined(i, j, k)) {
				df = -(fx(i + 1, j, k) - fx(i, j, k)) * h1inv;
				df -= (fy(i, j + 1, k) - fy(i, j, k)) * h1inv;
				df -= (fz(i, j, k + 1) - fz(i, j, k)) * h1inv;
				dphi(i, j, k) = (df - S(i, j, k));
				sum += fabs(dphi(i, j, k)) * h3;
			} else {
				dphi(i, j, k) = 0.0;
			}
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			sum += static_cast<Poisson*> (get_child(i))->solution_error();
		}
	}
	return sum;
}

Poisson::Poisson() :
	GridNode(), PoissonInterface() {
	const Indexer2d indexer(0, PNX - 1, 0, PNX - 1);
	int i, j, k;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i,j,k)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 0; i < PNX; i++) {
			phi(i, j, k) = 0.0;
			dphi1(i, j, k) = 0.0;
			dphi(i, j, k) = 0.0;
			phi0(i, j, k) = 0.0;
		}
	}
}

Poisson::~Poisson() {
}

Real Poisson::get_phi(int i, int j, int k) const {
	return phi(i, j, k);
}

Real Poisson::get_dphi(int i, int j, int k) const {
	return dphi(i, j, k);
}

PoissonPhysBound* Poisson::new_phys_bound() const {
	return new PoissonPhysBound;
}

PoissonAMRBound* Poisson::new_amr_bound() const {
	return new PoissonAMRBound;
}

void Poisson::phi_from_children() {

	for (int i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->phi_from_children();
		}
	}

	ChildIndex c;
	Poisson* child;
	for (c = 0; c < OCT_NCHILD; c++) {
		child = static_cast<Poisson*> (get_child(c));
		const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
		int k, j, k0, j0, i, i0;
		int n, o, p, ii;
		Real sum;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,k0,j0,i,i0,ii,sum,n,o,p)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
			j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
			if (child != NULL) {
				for (i = 1; i < PNX - 1; i += 2) {
					i0 = (1 + i) / 2 + c.get_x() * (PNX / 2 - 1);
					sum = 0.0;
					ii = 0;
					for (n = -1; n < 3; n++) {
						for (o = -1; o < 3; o++) {
							for (p = -1; p < 3; p++) {
								sum += interp_c[ii] * child->phi(i + p, j + o, k + n);
								ii++;
							}
						}
					}
					phi(i0, j0, k0) = sum / interp_c0;
				}
			}
		}
	}
}

void Poisson::enforce_phi_boundaries(int l) {
	if (get_level() == l) {
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int i, j, k;
		const PoissonInterface* sib;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i,j,k,sib)
		for (int index = 0; index <= indexer.max_index(); index++) {
			j = indexer.x(index);
			k = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XL));
			phi(0, j, k) = sib->get_phi(PNX - 2, j, k);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XU));
			phi(PNX - 1, j, k) = sib->get_phi(1, j, k);
			i = indexer.x(index);
			k = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YL));
			phi(i, 0, k) = sib->get_phi(i, PNX - 2, k);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YU));
			phi(i, PNX - 1, k) = sib->get_phi(i, 1, k);
			i = indexer.x(index);
			j = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZL));
			phi(i, j, 0) = sib->get_phi(i, j, PNX - 2);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZU));
			phi(i, j, PNX - 1) = sib->get_phi(i, j, 1);
		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->enforce_phi_boundaries(l);
			}
		}
	}
	if (get_level() == 0) {
		enforce_phi_boundaries_edge(l);
		enforce_phi_boundaries_vertex(l);
	}
}

void Poisson::enforce_phi_boundaries_edge(int l) {
	if (get_level() == l) {
		const PoissonInterface* sib;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(sib)
		for (int k = 1; k < PNX - 1; k++) {
			if (is_phys_bound(XL) || is_phys_bound(XU)) {
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(XL));
				phi(0, 0, k) = sib->get_phi(PNX - 2, 0, k);
				phi(0, PNX - 1, k) = sib->get_phi(PNX - 2, PNX - 1, k);
				phi(0, k, 0) = sib->get_phi(PNX - 2, k, 0);
				phi(0, k, PNX - 1) = sib->get_phi(PNX - 2, k, PNX - 1);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(XU));
				phi(PNX - 1, 0, k) = sib->get_phi(1, 0, k);
				phi(PNX - 1, PNX - 1, k) = sib->get_phi(1, PNX - 1, k);
				phi(PNX - 1, k, 0) = sib->get_phi(1, k, 0);
				phi(PNX - 1, k, PNX - 1) = sib->get_phi(1, k, PNX - 1);
			} else {
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(YL));
				phi(0, 0, k) = sib->get_phi(0, PNX - 2, k);
				phi(PNX - 1, 0, k) = sib->get_phi(PNX - 1, PNX - 2, k);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(YU));
				phi(0, PNX - 1, k) = sib->get_phi(0, 1, k);
				phi(PNX - 1, PNX - 1, k) = sib->get_phi(PNX - 1, 1, k);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZL));
				phi(0, k, 0) = sib->get_phi(0, k, PNX - 2);
				phi(PNX - 1, k, 0) = sib->get_phi(PNX - 1, k, PNX - 2);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZU));
				phi(0, k, PNX - 1) = sib->get_phi(0, k, 1);
				phi(PNX - 1, k, PNX - 1) = sib->get_phi(PNX - 1, k, 1);
			}
			if (is_phys_bound(YL) || is_phys_bound(YU)) {
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(YL));
				phi(k, 0, 0) = sib->get_phi(k, PNX - 2, 0);
				phi(k, 0, PNX - 1) = sib->get_phi(k, PNX - 2, PNX - 1);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(YU));
				phi(k, PNX - 1, 0) = sib->get_phi(k, 1, 0);
				phi(k, PNX - 1, PNX - 1) = sib->get_phi(k, 1, PNX - 1);
			} else {
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZL));
				phi(k, 0, 0) = sib->get_phi(k, 0, PNX - 2);
				phi(k, PNX - 1, 0) = sib->get_phi(k, PNX - 1, PNX - 2);
				sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZU));
				phi(k, 0, PNX - 1) = sib->get_phi(k, 0, 1);
				phi(k, PNX - 1, PNX - 1) = sib->get_phi(k, PNX - 1, 1);
			}

		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->enforce_phi_boundaries_edge(l);
			}
		}
	}
}

void Poisson::enforce_phi_boundaries_vertex(int l) {
	if (get_level() == l) {
		const PoissonInterface* sib;
		if (is_phys_bound(XL) || is_phys_bound(XU)) {
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XL));
			phi(0, 0, 0) = sib->get_phi(PNX - 2, 0, 0);
			phi(0, PNX - 1, 0) = sib->get_phi(PNX - 2, PNX - 1, 0);
			phi(0, 0, PNX - 1) = sib->get_phi(PNX - 2, 0, PNX - 1);
			phi(0, PNX - 1, PNX - 1) = sib->get_phi(PNX - 2, PNX - 1, PNX - 1);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XU));
			phi(PNX - 1, 0, 0) = sib->get_phi(1, 0, 0);
			phi(PNX - 1, PNX - 1, 0) = sib->get_phi(1, PNX - 1, 0);
			phi(PNX - 1, 0, PNX - 1) = sib->get_phi(1, 0, PNX - 1);
			phi(PNX - 1, PNX - 1, PNX - 1) = sib->get_phi(1, PNX - 1, PNX - 1);
		} else if (is_phys_bound(YL) || is_phys_bound(YU)) {
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YL));
			phi(0, 0, 0) = sib->get_phi(0, PNX - 2, 0);
			phi(PNX - 1, 0, 0) = sib->get_phi(PNX - 1, PNX - 2, 0);
			phi(0, 0, PNX - 1) = sib->get_phi(0, PNX - 2, PNX - 1);
			phi(PNX - 1, 0, PNX - 1) = sib->get_phi(PNX - 1, PNX - 2, PNX - 1);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YU));
			phi(0, PNX - 1, 0) = sib->get_phi(0, 1, 0);
			phi(PNX - 1, PNX - 1, 0) = sib->get_phi(PNX - 1, 1, 0);
			phi(0, PNX - 1, PNX - 1) = sib->get_phi(0, 1, PNX - 1);
			phi(PNX - 1, PNX - 1, PNX - 1) = sib->get_phi(PNX - 1, 1, PNX - 1);
		} else {
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZL));
			phi(0, 0, 0) = sib->get_phi(0, 0, PNX - 2);
			phi(0, PNX - 1, 0) = sib->get_phi(0, PNX - 1, PNX - 2);
			phi(PNX - 1, 0, 0) = sib->get_phi(PNX - 1, 0, PNX - 2);
			phi(PNX - 1, PNX - 1, 0) = sib->get_phi(PNX - 1, PNX - 1, PNX - 2);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZU));
			phi(0, 0, PNX - 1) = sib->get_phi(0, 0, 1);
			phi(0, PNX - 1, PNX - 1) = sib->get_phi(0, PNX - 1, 1);
			phi(PNX - 1, 0, PNX - 1) = sib->get_phi(PNX - 1, 0, 1);
			phi(PNX - 1, PNX - 1, PNX - 1) = sib->get_phi(PNX - 1, PNX - 1, 1);
		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->enforce_phi_boundaries_vertex(l);
			}
		}
	}
}

void Poisson::enforce_dphi_boundaries(int l) {
	if (get_level() == l) {
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int i, j, k;
		const PoissonInterface* sib;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(i,j,k,sib)
		for (int index = 0; index <= indexer.max_index(); index++) {
			j = indexer.x(index);
			k = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XL));
			assert(sib!=NULL);
			dphi(0, j, k) = sib->get_dphi(PNX - 2, j, k);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(XU));
			assert(sib!=NULL);
			dphi(PNX - 1, j, k) = sib->get_dphi(1, j, k);
			i = indexer.x(index);
			k = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YL));
			assert(sib!=NULL);
			dphi(i, 0, k) = sib->get_dphi(i, PNX - 2, k);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(YU));
			assert(sib!=NULL);
			dphi(i, PNX - 1, k) = sib->get_dphi(i, 1, k);
			i = indexer.x(index);
			j = indexer.y(index);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZL));
			assert(sib!=NULL);
			dphi(i, j, 0) = sib->get_dphi(i, j, PNX - 2);
			sib = dynamic_cast<const PoissonInterface*> (get_sibling(ZU));
			assert(sib!=NULL);
			dphi(i, j, PNX - 1) = sib->get_dphi(i, j, 1);
		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->enforce_dphi_boundaries(l);
			}
		}
	}
}

void Poisson::init_vcycle() {
	const Indexer2d indexer(0, PNX - 1, 0, PNX - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < PNX; i++) {
			dphi(i, j, k) = 0.0;
			dphi1(i, j, k) = 0.0;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->init_vcycle();
		}
	}
}

void Poisson::set_source() {
	const Indexer2d indexer(0, PNX - 1, 0, PNX - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < PNX; i++) {
			S(i, j, k) = (*this)(i + BW - 1, j + BW - 1, k + BW - 1).poisson_source();
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Poisson*> (get_child(i))->set_source();
		}
	}
}

bool Poisson::poisson_zone_is_refined(int i, int j, int k) const {
	ChildIndex c;
	i = (2 * i) / PNX;
	j = (2 * j) / PNX;
	k = (2 * k) / PNX;
	c.set_index(i, j, k);
	return (get_child(c) != NULL);
}

Real Poisson::pxc(int i) const {
	return xc(i + BW - 1);
}

Real Poisson::pyc(int i) const {
	return yc(i + BW - 1);
}

Real Poisson::pzc(int i) const {
	return zc(i + BW - 1);
}

void Poisson::vcycle_down(int level) {
	if (get_level() == level) {
		ChildIndex c;
		Poisson* child;
		const Real h2 = get_dx() * get_dx();
		const Indexer2d_by2 indexer(1, PNX - 2, 1, PNX - 2);
		int k, j, i, k0, j0, i0, i0min, i0max;
		for (c = 0; c < OCT_NCHILD; c++) {
			child = static_cast<Poisson*> (get_child(c));
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,k0,j0,i0,i0min,i0max)
			for (int index = 0; index <= indexer.max_index(); index++) {
				k = indexer.y(index);
				j = indexer.x(index);
				k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
				j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
				if (child == NULL) {
					i0min = 1 + c.get_x() * (PNX / 2 - 1);
					i0max = PNX / 2 + c.get_x() * (PNX / 2 - 1) - 1;
					for (i0 = i0min; i0 <= i0max; i0++) {
						dphi(i0, j0, k0) = (phi.divergence(i0, j0, k0) - h2 * S(i0, j0, k0)) * (1.0 / 6.0);
					}
				} else {
					for (i = 1; i < PNX - 1; i += 2) {
						i0 = (1 + i) / 2 + c.get_x() * (PNX / 2 - 1);
						dphi1(i0, j0, k0) = child->dphi1.oct_avg(i, j, k);
						dphi(i0, j0, k0) = child->dphi.oct_avg(i, j, k);
					}
				}
			}
		}
		adjust_x_dphi();
		adjust_y_dphi();
		adjust_z_dphi();
		for (c = 0; c < OCT_NCHILD; c++) {
			child = static_cast<Poisson*> (get_child(c));
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,k0,j0,i0,i0min,i0max)
			for (int index = 0; index <= indexer.max_index(); index++) {
				k = indexer.y(index);
				j = indexer.x(index);
				k0 = (1 + k) / 2 + c.get_z() * (PNX / 2 - 1);
				j0 = (1 + j) / 2 + c.get_y() * (PNX / 2 - 1);
				if (child == NULL) {
					i0min = 1 + c.get_x() * (PNX / 2 - 1);
					i0max = PNX / 2 + c.get_x() * (PNX / 2 - 1) - 1;
					for (i0 = i0min; i0 <= i0max; i0++) {
						dphi1(i0, j0, k0) = dphi(i0, j0, k0);
					}
				}
			}
		}
	} else {
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->vcycle_down(level);
			}
		}
	}
}

void Poisson::vcycle_retire_dphi(int level) {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int i, j, k;
	if (get_level() == level) {
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = 1; i < PNX - 1; i++) {
				phi(i, j, k) += dphi(i, j, k);
			}
		}
	} else {
		for (i = 0; i < OCT_NCHILD; i++) {
			if (this->get_child(i) != NULL) {
				static_cast<Poisson*> (get_child(i))->vcycle_retire_dphi(level);
			}
		}
	}
}

Real Poisson::vcycle() {
	set_source();
	vcycle0();
	compute_forces();
	return solution_error();
}
