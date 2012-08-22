#include "../defs.h"

#ifdef SINGLE
#define ZSTAR

#include "single.h"
#include "indexer_2d.h"

void Single::write(FILE* fp) const {
	if (get_level() == 0) {
		Real h;
		bool dummy;
		fwrite(&dummy, sizeof(bool), 1, fp);
		fwrite(&h, sizeof(Real), 1, fp);
	}
	GridNode::write(fp);
}

void Single::read(FILE* fp) {
	if (get_level() == 0) {
		Real h;
		bool dummy;
		fread(&dummy, sizeof(bool), 1, fp);
		fread(&h, sizeof(Real), 1, fp);
	}
	GridNode::read(fp);
}

void Single::pot_to_grid() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i, j0, k0, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,j0,k0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		j0 = j - 1 + BW;
		k0 = k - 1 + BW;
		for (i = 1, i0 = BW; i < PNX - 1; i++, i0++) {
			(*this)(i0, j0, k0)[State::pot_index] = (phi(i, j, k) * (*this)(i0, j0, k0).rho());
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->pot_to_grid();
		}
	}
}

void Single::add_difs(Real a, Real b) {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	Vector<Real, STATE_NF>* D;
	int j, k, i, i0, k0, j0;
	State s;
	Real dsum, fx0, fy0, fz0, x, y;
	Real angular;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,D,s,dsum,i0,j0,k0,fx0,fy0,fz0,x,y,angular)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			i0 = i + BW - 1;
			j0 = j + BW - 1;
			k0 = k + BW - 1;
			D = &differential(i0, j0, k0);
			s = (*this)(i0, j0, k0);
			fx0 = 0.5 * s.rho() * (fx(i, j, k) + fx(i + 1, j, k));
			fy0 = 0.5 * s.rho() * (fy(i, j, k) + fy(i, j + 1, k));
			fz0 = 0.5 * s.rho() * (fz(i, j, k) + fz(i, j, k + 1));
			(*D)[State::sz_index] += fz0;
#ifndef USE_LZ
			(*D)[State::sx_index] += fx0;
			(*D)[State::sy_index] += fy0;
#else
			(*D)[State::sy_index] += pxc(i) * fy0 - pyc(j) * fx0;
			(*D)[State::sx_index] += (pxc(i) * fx0 + pyc(j) * fy0) / sqrt(pxc(i) * pxc(i) + pyc(j) * pyc(j));

#endif
#ifndef ZSTAR
			(*D)[State::et_index] += (*D)[State::pot_index];
			dsum = 0.0;
			for (int l = 0; l < NRHO; l++) {
				dsum += (*D)[State::d_index + l];
			}
			(*D)[State::et_index] -= (s.pot() / s.rho()) * dsum;
			(*D)[State::pot_index] = 0.0;
#endif
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->add_difs(a, b);
		}
	}
	if (get_level() == 0) {
		GridNode::add_difs(a, b);
	}
}

Real Single::pot_sum() const {
	Real h3 = pow(get_dx(), 3);
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	Real a0, ds;
	int k, j, i;
	a0 = 0.0;
#pragma omp parallel for schedule(OMP_SCHEDULE) reduction(+:a0) private(ds,k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 1; i < PNX - 1; i++) {
			if (!poisson_zone_is_refined(i, j, k)) {
				ds = 0.5 * phi(i, j, k) * (*this)(i + BW - 1, j + BW - 1, k + BW - 1).rho() * h3;
				a0 += ds;
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			a0 += static_cast<const Single*> (get_child(i))->pot_sum();
		}
	}
	return a0;
}

void Single::inject_from_parent(ChildIndex c) {
	Poisson::inject_from_parent(c);
	pot_from_grid();
}

void Single::pot_from_grid() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i, j0, k0, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,j0,k0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		j0 = j - 1 + BW;
		k0 = k - 1 + BW;
		for (i = 1, i0 = BW; i < PNX - 1; i++, i0++) {
			phi(i, j, k) = (*this)(i0, j0, k0).pot() / (*this)(i0, j0, k0).rho();
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->pot_from_grid();
		}
	}
	if (get_level() == 0) {
		const int l1 = max_level();
		for (int l = 0; l <= l1; l++) {
			enforce_phi_boundaries(l);
		}
	}
}

void Single::add_pot_et() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int k, j, i;
	Real e;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,e)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 1; i < PNX - 1; i++) {
			const int i0 = i + BW - 1;
			const int j0 = j + BW - 1;
			const int k0 = k + BW - 1;
			e = (*this)(i0, j0, k0).et();
			e += 0.5 * (*this)(i0, j0, k0).pot();
			(*this)(i0, j0, k0).set_et(e);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->add_pot_et();
		}
	}
}

Single* Single::new_octnode() const {
	return new Single;
}

void Single::sub_pot_et() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int k, j, i;
	Real e;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,e)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 1; i < PNX - 1; i++) {
			const int i0 = i + BW - 1;
			const int j0 = j + BW - 1;
			const int k0 = k + BW - 1;
			e = (*this)(i0, j0, k0).et();
			e -= 0.5 * (*this)(i0, j0, k0).pot();
			(*this)(i0, j0, k0).set_et(e);
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->sub_pot_et();
		}
	}
}

void Single::compute_physical_boundaries(Single* root) {
	for (int f = 0; f < OCT_NSIB; f++) {
		if (is_phys_bound(f)) {
			static_cast<PoissonPhysBound*> (get_sibling(f))->compute_boundaries(root);
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->compute_physical_boundaries(root);
		}
	}

}

Real Single::compute_phi(Real x, Real y, Real z) const {
	Real sum, dx, dy, dz, alpha;
	alpha = 1.0;
#ifdef Z_REFLECT
	alpha *= 2.0;
#endif
	sum = 0.0;
	for (int i = 0; i < 8; i++) {
		dx = x - src_x[i];
		dy = y - src_y[i];
		dz = z - src_z[i];
		if (!all_refined) {
			sum -= alpha * src_0[i] / sqrt(dx * dx + dy * dy + dz * dz);
		} else {
			sum = 0.0;
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			sum += static_cast<const Single*> (get_child(i))->compute_phi(x, y, z);
		}
	}
	return sum;
}

void Single::compute_source_sums() {
	if (get_level() == 0) {
		set_source();
	}
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	const Real h3 = get_dx() * get_dx() * get_dx();
	Real s0[8], sx[8], sy[8], sz[8];
	int k, j, i;
	for (int i = 0; i < 8; i++) {
		s0[i] = sx[i] = sy[i] = sz[i] = 0.0;
	}
	all_refined = true;
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 1; i < PNX - 1; i++) {
			if (!poisson_zone_is_refined(i, j, k)) {
				all_refined = false;
				const int ci = i / (PNX / 2) + ((j / (PNX / 2)) << 1) + ((k / (PNX / 2)) << 2);
				s0[ci] += S(i, j, k) * h3;
				sx[ci] += S(i, j, k) * pxc(i) * h3;
				sy[ci] += S(i, j, k) * pyc(j) * h3;
				sz[ci] += S(i, j, k) * pzc(k) * h3;
			}
		}
	}
	for (int i = 0; i < 8; i++) {
		if (!all_refined && s0[i] != 0.0) {
			src_x[i] = sx[i] / s0[i];
			src_y[i] = sy[i] / s0[i];
			src_z[i] = sz[i] / s0[i];
			src_0[i] = s0[i] / (4.0 * M_PI);
		} else {
			src_0[i] = 0.0;
			src_x[i] = src_y[i] = src_z[i] = 1.0;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->compute_source_sums();
		}
	}
}

Real Single::max_dt() const {
	return GridNode::max_dt();
}

Single::Single() :
	Poisson() {
	// TODO Auto-generated constructor stub

}

Single::~Single() {
	// TODO Auto-generated destructor stub
}

void Single::sum_x_difs(bool scalar) {
	Grid::sum_x_difs(scalar);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->sum_x_difs(scalar);
		}
	}
#ifdef ZSTAR
	if (!scalar) {
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = 1; i < PNX - 1; i++) {
				const int i0 = i + BW - 1;
				const int j0 = j + BW - 1;
				const int k0 = k + BW - 1;
				D(i0, j0, k0)[State::et_index] += 0.5 * fx(i + 0, j, k) * (F(i0 + 0, j0, k0)[State::d_index]);
				D(i0, j0, k0)[State::et_index] += 0.5 * fx(i + 1, j, k) * (F(i0 + 1, j0, k0)[State::d_index]);
			}
		}
	}
#endif
}

void Single::sum_y_difs(bool scalar) {
	Grid::sum_y_difs(scalar);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->sum_y_difs(scalar);
		}
	}
#ifdef ZSTAR
	if (!scalar) {
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = 1; i < PNX - 1; i++) {
				const int i0 = i + BW - 1;
				const int j0 = j + BW - 1;
				const int k0 = k + BW - 1;
				D(i0, j0, k0)[State::et_index] += 0.5 * fy(i, j + 0, k) * (F(i0, j0 + 0, k0)[State::d_index]);
				D(i0, j0, k0)[State::et_index] += 0.5 * fy(i, j + 1, k) * (F(i0, j0 + 1, k0)[State::d_index]);
			}
		}
	}
#endif
}

void Single::sum_z_difs(bool scalar) {
	Grid::sum_z_difs(scalar);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Single*> (get_child(i))->sum_z_difs(scalar);
		}
	}
#ifdef ZSTAR
	if (!scalar) {
		const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
		int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			k = indexer.y(index);
			j = indexer.x(index);
			for (i = 1; i < PNX - 1; i++) {
				const int i0 = i + BW - 1;
				const int j0 = j + BW - 1;
				const int k0 = k + BW - 1;
				D(i0, j0, k0)[State::et_index] += 0.5 * fz(i, j, k + 0) * (F(i0, j0, k0 + 0)[State::d_index]);
				D(i0, j0, k0)[State::et_index] += 0.5 * fz(i, j, k + 1) * (F(i0, j0, k0 + 1)[State::d_index]);
			}
		}
	}
#endif
}

#endif
