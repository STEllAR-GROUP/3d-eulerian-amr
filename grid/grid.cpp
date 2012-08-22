#include "grid_output.h"
#include "grid.h"
#include "indexer_2d.h"
#include <stdlib.h>
#include <stdio.h>
#include <silo.h>
#include <omp.h>

void Grid::set_origin(const _3Vec& v) {
	origin = v;
}

_3Vec Grid::X(int i, int j, int k) const {
	_3Vec x;
	x[0] = xc(i);
	x[1] = yc(j);
	x[2] = zc(k);
	return x;

}

_3Vec Grid::Xfx(int i, int j, int k) const {
	_3Vec x;
	x[0] = xf(i);
	x[1] = yc(j);
	x[2] = zc(k);
	return x;

}

_3Vec Grid::Xfy(int i, int j, int k) const {
	_3Vec x;
	x[0] = xc(i);
	x[1] = yf(j);
	x[2] = zc(k);
	return x;

}

_3Vec Grid::Xfz(int i, int j, int k) const {
	_3Vec x;
	x[0] = xc(i);
	x[1] = yc(j);
	x[2] = zf(k);
	return x;

}

_3Vec Grid::V(int i, int j, int k) const {
	_3Vec v = 0.0;
	return v;
}

_3Vec Grid::Vfx(int i, int j, int k) const {
	_3Vec v = 0.0;
	return v;
}

_3Vec Grid::Vfy(int i, int j, int k) const {
	_3Vec v = 0.0;
	return v;
}

_3Vec Grid::Vfz(int i, int j, int k) const {
	_3Vec v = 0.0;
	return v;
}

Reconstruct Grid::reconstruct;
Real Grid::h0 = (2.0 * GRID_DIM / Real(GNX - 2 * BW));
Real Grid::v0x = 0.0;
Real Grid::v0y = 0.0;
Real Grid::v0z = 0.0;

Vector<Real, STATE_NF>& Grid::differential(int i, int j, int k) {
	return D(i, j, k);
}

Vector<Real, STATE_NF> Grid::differential(int i, int j, int k) const {
	return D(i, j, k);
}

Vector<Real, STATE_NF> Grid::getD(int j, int k, int l) const {
	return D(j, k, l);

}

void Grid::incD(const Vector<Real, STATE_NF>& v, int i, int j, int k) {
	D(i, j, k) += v;
}

Real Grid::rc(int x, int y, int z) const {
	return sqrt(xc(x) * xc(x) + yc(y) * yc(y) + zc(z) * zc(z));
}

void Grid::add_to_dif(const Vector<Real, STATE_NF>& v, int x, int y, int z) {
	D(x, y, z) += v;
}

void Grid::set_flux(const Vector<Real, STATE_NF>& f, int a, int b, int c) {
	F(a, b, c) = f;
}

void Grid::write(FILE* fp) const {
	fwrite(&h, sizeof(Real), 1, fp);
	fwrite(&h0, sizeof(Real), 1, fp);
	fwrite(&time, sizeof(Real), 1, fp);
	for (int l = 0; l < 3; l++) {
		fwrite(&(offset[l]), sizeof(int), 1, fp);
		fwrite(&(origin[l]), sizeof(Real), 1, fp);
	}
	for (int l = 0; l < STATE_NF; l++) {
		fwrite(&(FO[l]), sizeof(Real), 1, fp);
	}
	for (int i = 0; i < GNX; i++) {
		for (int j = 0; j < GNX; j++) {
			for (int k = 0; k < GNX; k++) {
				for (int l = 0; l < STATE_NF; l++) {
					fwrite(&((*this)(i, j, k)[l]), sizeof(Real), 1, fp);
				}
			}
		}
	}
}

void Grid::read(FILE* fp) {
	fread(&h, sizeof(Real), 1, fp);
	fread(&h0, sizeof(Real), 1, fp);
	fread(&time, sizeof(Real), 1, fp);
	for (int l = 0; l < 3; l++) {
		fread(&(offset[l]), sizeof(int), 1, fp);
		fread(&(origin[l]), sizeof(Real), 1, fp);
	}
	for (int l = 0; l < STATE_NF; l++) {
		fread(&(FO[l]), sizeof(Real), 1, fp);
	}
	for (int i = 0; i < GNX; i++) {
		for (int j = 0; j < GNX; j++) {
			for (int k = 0; k < GNX; k++) {
				for (int l = 0; l < STATE_NF; l++) {
					fread(&((*this)(i, j, k)[l]), sizeof(Real), 1, fp);
				}
			}
		}
	}
}

const Vector<Real, STATE_NF> Grid::get_flux(int a, int b, int c) const {
	assert(a >= bw );assert(b >= bw );assert(c >= bw );assert( a <= GNX-bw );assert( b <= GNX-bw );assert( c <= GNX-bw );
	return F(a, b, c);
}

State Grid::operator()(int i, int j, int k) const {
	return U(i, j, k);
}

State& Grid::operator()(int i, int j, int k) {
	return U(i, j, k);
}

Real Grid::get_time() const {
	return time;
}

Real Grid::get_dx() const {
	return h;
}

const Vector<int, 3>& Grid::get_offset() const {
	return offset;
}

void Grid::store() {
	FO0 = FO;
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			U0(i, j, k) = U(i, j, k);
		}
	}
}

void Grid::restore() {
	FO = FO0;
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			U(i, j, k) = U0(i, j, k);
		}
	}
}

void Grid::clear_difs() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	DFO = Vector<Real, STATE_NF> (0.0);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			D(i, j, k) = Vector<Real, STATE_NF> (0.0);
		}
	}
}

void Grid::compute_x_flux() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State q0[GNX], ql[GNX], qr[GNX];
	int k, j, i;
	_3Vec x, v;
	Real a;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(q0,ql,qr,a,k,j,i,x,v)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			q0[i] = U(i, j, k);
			q0[i].to_prim(this->X(i, j, k));
		}
		reconstruct(q0, ql, qr);
		for (i = bw; i < GNX - bw + 1; i++) {
			v = this->Vfx(i, j, k);
			x = this->Xfx(i, j, k);
			ql[i].to_con(x);
			qr[i].to_con(x);
			a = max(ql[i].max_abs_x_eigen(x), qr[i].max_abs_x_eigen(x));
			F(i, j, k) = ((ql[i].x_flux(x) + qr[i].x_flux(x)) - (qr[i] - ql[i]) * a) * 0.5;

		}
	}
}

void Grid::compute_y_flux() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State q0[GNX], ql[GNX], qr[GNX];
	int k, j, i;
	_3Vec x, v;
	Real a;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(q0,ql,qr,a,k,j,i,x,v)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		i = indexer.x(index);
		for (j = 0; j < GNX; j++) {
			q0[j] = U(i, j, k);
			q0[j].to_prim(this->X(i, j, k));
		}
		reconstruct(q0, ql, qr);
		for (j = bw; j < GNX - bw + 1; j++) {
			v = this->Vfy(i, j, k);
			x = this->Xfy(i, j, k);
			ql[j].to_con(x);
			qr[j].to_con(x);
			a = max(ql[j].max_abs_y_eigen(x), qr[j].max_abs_y_eigen(x));
			F(i, j, k) = ((ql[j].y_flux(x) + qr[j].y_flux(x)) - (qr[j] - ql[j]) * a) * 0.5;

		}
	}
}

void Grid::compute_z_flux() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State q0[GNX], ql[GNX], qr[GNX];
	int k, j, i;
	Real a;
	_3Vec x, v;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(q0,ql,qr,a,k,j,i,x,v)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.y(index);
		i = indexer.x(index);
		for (k = 0; k < GNX; k++) {
			q0[k] = U(i, j, k);
			q0[k].to_prim(this->X(i, j, k));
		}
		reconstruct(q0, ql, qr);
		for (k = bw; k < GNX - bw + 1; k++) {
			v = this->Vfz(i, j, k);
			x = this->Xfz(i, j, k);
			ql[k].to_con(x);
			qr[k].to_con(x);
			a = max(ql[k].max_abs_z_eigen(x), qr[k].max_abs_z_eigen(x));
			F(i, j, k) = ((ql[k].z_flux(x) + qr[k].z_flux(x)) - (qr[k] - ql[k]) * a) * 0.5;
		}
	}
}

void Grid::sum_x_difs() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	const Real hinv = 1.0 / h;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			D(i, j, k) -= (F(i + 1, j, k) - F(i, j, k)) * hinv;
		}
	}
}

void Grid::sum_y_difs() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	const Real hinv = 1.0 / h;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			D(i, j, k) -= (F(i, j + 1, k) - F(i, j, k)) * hinv;
		}
	}
}

void Grid::sum_z_difs() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	const Real hinv = 1.0 / h;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			D(i, j, k) -= (F(i, j, k + 1) - F(i, j, k)) * hinv;
		}
	}
}

State Grid::flow_off() const {
	return FO;
}

void Grid::sum_x_flow_off() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State tmp;
	int k, j;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(tmp,k,j)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		tmp = (F(GNX - bw, j, k) - F(bw, j, k)) * h * h;
#pragma omp critical
		DFO += tmp;
	}
}

void Grid::sum_y_flow_off() {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State tmp;
	int k, j;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(tmp,k,j)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		tmp = (F(j, GNX - bw, k) - F(j, bw, k)) * h * h;
#pragma omp critical
		DFO += tmp;
	}
}

void Grid::sum_z_flow_off() {
#ifdef Z_REFLECT
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State tmp;
	int k, j;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(tmp,k,j)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		tmp = (F(j, k, GNX - bw)) * h * h;// * 2.0;
#pragma omp critical
		DFO += tmp;
	}
#else
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	State tmp;
	int k, j;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(tmp,k,j)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		tmp = (F(j, k, GNX - bw) - F(j, k, bw)) * h * h;
#pragma omp critical
		DFO += tmp;
	}
#endif
}

void Grid::add_difs(Real dt, Real beta) {
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			D(i, j, k) += U(i, j, k).source(this->X(i, j, k));
			U(i, j, k) = (U(i, j, k) + D(i, j, k) * dt) * beta + U0(i, j, k) * (1.0 - beta);
			U(i, j, k).floor(this->X(i, j, k));
		}
	}
	FO = (FO + DFO * dt) * beta + FO0 * (1.0 - beta);
}

Real Grid::max_dt() const {
	Real dtinv = 0.0;
	Real this_dt;
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	_3Vec x, v;
	State q0[GNX], ql[GNX], qr[GNX];
	int k, j, i;
	Real this_dtinv;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(q0,ql,qr,k,j,i,this_dtinv,x,v)
	for (int index = 0; index <= indexer.max_index(); index++) {
		this_dtinv = 0.0;
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = 0; i < GNX; i++) {
			q0[i] = U(i, j, k);
			q0[i].to_prim(this->X(i, j, k));
		}
		reconstruct(q0, ql, qr);
		for (i = bw; i < GNX - bw + 1; i++) {
			x = this->Xfx(i, j, k);
			v = this->Vfx(i, j, k);
			ql[i].to_con(x);
			qr[i].to_con(x);
			this_dtinv = max(this_dtinv, ql[i].max_abs_x_eigen(x), qr[i].max_abs_x_eigen(x));
		}
		for (i = 0; i < GNX; i++) {
			q0[i] = U(j, i, k);
			q0[i].to_prim(this->X(j, i, k));
		}
		reconstruct(q0, ql, qr);
		for (i = bw; i < GNX - bw + 1; i++) {
			v = this->Vfy(j, i, k);
			x = this->Xfy(j, i, k);
			ql[i].to_con(x);
			qr[i].to_con(x);
			this_dtinv = max(this_dtinv, ql[i].max_abs_y_eigen(x), qr[i].max_abs_y_eigen(x));
		}
		for (i = 0; i < GNX; i++) {
			q0[i] = U(j, k, i);
			q0[i].to_prim(this->X(j, k, i));
		}
		reconstruct(q0, ql, qr);
		for (i = bw; i < GNX - bw + 1; i++) {
			v = this->Vfz(j, k, i);
			x = this->Xfz(j, k, i);
			ql[i].to_con(x);
			qr[i].to_con(x);
			this_dtinv = max(this_dtinv, ql[i].max_abs_z_eigen(x), qr[i].max_abs_z_eigen(x));
		}
#pragma omp critical
		dtinv = max(this_dtinv, dtinv);
	}
	dtinv /= h;
	if (dtinv == 0.0) {
		dtinv = 1.0E-10;
	}
	this_dt = 1.0 / dtinv;
	this_dt *= GRID_CFL_FACTOR;

	return this_dt;

}

Grid::Grid() :
	GridInterface() {
	offset = 0;
	h = h0;
	time = 0.0;
	FO0 = Vector<Real, STATE_NF> (0.0);
	FO = Vector<Real, STATE_NF> (0.0);
	DFO = Vector<Real, STATE_NF> (0.0);
	origin = 0.0;
}

Grid::~Grid() {
	return;
}

Real Grid::xc(int i) const {
	return xf(i) + 0.5 * h;
}

Real Grid::yc(int j) const {
	return yf(j) + 0.5 * h;
}

Real Grid::zc(int k) const {
	return zf(k) + 0.5 * h;
}

Real Grid::xf(int i) const {
	return Real(offset[0] + i) * h - GRID_DIM - bw * h0 - origin[0];
}

Real Grid::yf(int i) const {
	return Real(offset[1] + i) * h - GRID_DIM - bw * h0 - origin[1];
}

Real Grid::zf(int i) const {
#ifndef Z_REFLECT
	return Real(offset[2] + i) * h - GRID_DIM - bw * h0;
#else
	return Real(offset[2] + i) * h - bw * h0 - origin[2];
#endif
}

void Grid::set_time(Real a) {
	time = a;
}

void Grid::set_offset(const Vector<int, 3>& a) {
	offset = a;
}

void Grid::set_dx(Real a) {
	h = a;

}
