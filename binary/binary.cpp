#include "../defs.h"

#ifdef BINARY

#include "binary.h"
#include "indexer_2d.h"

int Binary::accretor_type = 1, Binary::donor_type = 1;
Real Binary::Ug, Binary::Us, Binary::UK, Binary::Ucm;
Real Binary::AccretorMass = 1.0;
Real Binary::M1 = 2.0e-3;
Real Binary::fill_factor = 0.99;
Real Binary::n = 1.5;
Real Binary::R1 = 2.5e-2;
Real Binary::R2 = 2.5e-2;
Real Binary::a = 1.0;
Real Binary::q = 0.20;
Real Binary::K1 = -1.0;
Real Binary::K2 = -1.0;
Real Binary::phi0_1 = 0.0;
Real Binary::phi0_2 = 0.0;
Real Binary::scf_dt = 0.2;
Real Binary::phase = 0.0;
Real Binary::Omega0 = 0.0;
Real Binary::l1x = 0.0, Binary::accretor_cut, Binary::donor_cut;
Real Binary::Omega = sqrt((M1 * (q + 1.0) / (a * a * a)));

Vector<Real, 4> Binary::mass_sum(int frac) const {
	const Real h3 = pow(get_dx(), 3);
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	Vector<Real, 4> a;
	int k, j, i;
	Real a0, ax, ay, az, d;
	a0 = ax = ay = az = 0.0;
#pragma omp parallel for schedule(OMP_SCHEDULE) reduction(+:a0,ax,ay,az) private(k,j,i,d)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw; i++) {
			if (!zone_is_refined(i, j, k)) {
				d = frac == 0 ? (*this)(i, j, k).rho_accretor() : (*this)(i, j, k).rho_donor();
				a0 += d * h3;
				ax += d * h3 * xc(i);
				ay += d * h3 * yc(j);
				az += d * h3 * zc(k);
			}
		}
	}
	a[0] = a0;
	a[1] = ax;
	a[2] = ay;
#ifdef Z_REFLECT
	a[3] = 0.0;
#else
	a[3] = az;
#endif
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			a += static_cast<const Binary*>(get_child(i))->mass_sum(frac);
		}
	}
	if (get_level() == 0) {
		a[1] /= a[0];
		a[2] /= a[0];
		a[3] /= a[0];
	}
	return a;
}

void Binary::write(FILE* fp) const {
	if (get_level() == 0) {
		Real h;
		bool dummy;
		fwrite(&dummy, sizeof(bool), 1, fp);
		fwrite(&l1x, sizeof(Real), 1, fp);
		fwrite(&h, sizeof(Real), 1, fp);
		fwrite(&M1, sizeof(Real), 1, fp);
		fwrite(&fill_factor, sizeof(Real), 1, fp);
		fwrite(&n, sizeof(Real), 1, fp);
		fwrite(&R1, sizeof(Real), 1, fp);
		fwrite(&q, sizeof(Real), 1, fp);
		fwrite(&a, sizeof(Real), 1, fp);
		fwrite(&K1, sizeof(Real), 1, fp);
		fwrite(&K2, sizeof(Real), 1, fp);
		fwrite(&phi0_1, sizeof(Real), 1, fp);
		fwrite(&phi0_2, sizeof(Real), 1, fp);
		fwrite(&Omega, sizeof(Real), 1, fp);
		fwrite(&phase, sizeof(Real), 1, fp);
		fwrite(&Omega0, sizeof(Real), 1, fp);
		fwrite(&Binary::accretor_type, sizeof(int), 1, fp);
		fwrite(&Binary::donor_type, sizeof(int), 1, fp);
	}
	GridNode::write(fp);
}

int Binary::nvar_output() const {
	return STATE_NF + 2;
}

void Binary::output(grid_output_t* ptr) const {
	for (int k = 0; k < GNX + 1; k++) {
		for (int j = 0; j < GNX + 1; j++) {
			for (int i = 0; i < GNX + 1; i++) {
				Real x, y, z;
				Real p;
#ifdef SCF_CODE
				p = 0.0;
#else
				p = Binary::phase;
#endif
				z = zf(k);
				x = cos(p) * xf(i) - sin(p) * yf(j);
				y = cos(p) * yf(j) + sin(p) * xf(i);
				*(ptr->x) = x;
				*(ptr->y) = y;
				*(ptr->z) = z;
				ptr->x++;
				ptr->y++;
				ptr->z++;
			}
		}
	}
	for (int k = bw; k < GNX - bw; k++) {
		for (int j = bw; j < GNX - bw; j++) {
			for (int i = bw; i < GNX - bw; i++) {
				if (zone_is_refined(i, j, k)) {
					continue;
				}
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 0)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 0) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 1) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 1)) + ptr->pi;
				(ptr->nodelist)[ptr->ni++] = (i + 0) + (GNX + 1) * ((j + 1) + (GNX + 1) * (k + 1)) + ptr->pi;
				this->load_output(ptr, i, j, k);
				ptr->ei++;
			}
		}
	}
	ptr->pi += GRID_NNODES;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<const Binary*>(get_child(i))->output(ptr);
		}
	}
}

void Binary::load_output(grid_output_t* go, int i, int j, int k) const {
	for (int l = 0; l < STATE_NF; l++) {
		go->ele[l][go->ei] = (*this)(i, j, k)[l];
	}
	go->ele[STATE_NF + 0][go->ei] = phi(i - BW + 1, j - BW + 1, k - BW + 1);
	go->ele[STATE_NF + 1][go->ei] = lobe(i, j, k);
}

const char* Binary::output_field_names(int i) const {
	if (i < STATE_NF) {
		return State::field_name(i);
	} else if (i == STATE_NF) {
		return "phi";
	} else if (i == STATE_NF + 1) {
		return "lobe";
	} else if (i == STATE_NF + 2) {
		return "dlz";
	} else {
		static char str[256];
		sprintf(str, "D_%s", State::field_name(i - STATE_NF));
		return str;
	}

}

void Binary::read(FILE* fp) {
	if (get_level() == 0) {
		Real h;
		bool dummy;
		fread(&dummy, sizeof(bool), 1, fp);
		fread(&l1x, sizeof(Real), 1, fp);
		fread(&h, sizeof(Real), 1, fp);
		fread(&M1, sizeof(Real), 1, fp);
		fread(&fill_factor, sizeof(Real), 1, fp);
		fread(&n, sizeof(Real), 1, fp);
		fread(&R1, sizeof(Real), 1, fp);
		fread(&q, sizeof(Real), 1, fp);
		fread(&a, sizeof(Real), 1, fp);
		fread(&K1, sizeof(Real), 1, fp);
		fread(&K2, sizeof(Real), 1, fp);
		fread(&phi0_1, sizeof(Real), 1, fp);
		fread(&phi0_2, sizeof(Real), 1, fp);
		fread(&Omega, sizeof(Real), 1, fp);
		fread(&phase, sizeof(Real), 1, fp);
		fread(&Omega0, sizeof(Real), 1, fp);
#ifdef ZTWD
		fread(&Binary::accretor_type, sizeof(int), 1, fp);
		fread(&Binary::donor_type, sizeof(int), 1, fp);
		//	State::initialize(Binary::accretor_type, Binary::donor_type);
#endif
		//	printf("K1,K2 %f %f\n", K1, K2);
	}
	GridNode::read(fp);
}

void Binary::add_difs(Real a, Real b) {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	Vector<Real, STATE_NF>* D;
	int j, k, i, i0, k0, j0;
	State s;
	Real fx0, fy0, fz0, du;
	const Real O = Binary::Omega0;
	_3Vec gforce;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,D,s,i0,j0,k0,gforce,du)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		for (i = 1; i < PNX - 1; i++) {
			i0 = i + BW - 1;
			j0 = j + BW - 1;
			k0 = k + BW - 1;
			D = &differential(i0, j0, k0);
			s = (*this)(i0, j0, k0);
#ifndef SCF_CODE
			gforce[0] = 0.5 * (fx(i, j, k) + fx(i + 1, j, k));
			gforce[1] = 0.5 * (fy(i, j, k) + fy(i, j + 1, k));
			gforce[2] = 0.5 * (fz(i, j, k) + fz(i, j, k + 1));
			(*D) += s.gravity_source(gforce, *D, X(i0, j0, k0));
			du = 0.0;
			du += Uxf(i0 + 1, j0, k0).vx(Xfx(i0 + 1, j0, k0));
			du += Uyf(i0, j0 + 1, k0).vy(Xfy(i0, j0 + 1, k0));
			du += Uzf(i0, j0, k0 + 1).vz(Xfz(i0, j0, k0 + 1));
			du -= Uxf(i0, j0, k0).vx(Xfx(i0, j0, k0));
			du -= Uyf(i0, j0, k0).vy(Xfy(i0, j0, k0));
			du -= Uzf(i0, j0, k0).vz(Xfz(i0, j0, k0));
			du /= get_dx();
			(*D) += s.internal_energy_source(X(i, j, k), du);
#else
			(*D) += s.scf_source(pxc(i), pyc(j), lobe(i0, j0, k0));
#endif
			(*D)[State::pot_index] = 0.0;
#ifdef ZTWD
			(*D)[State::T_index] = 0.0;
#endif
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->add_difs(a, b);
		}
	}
	if (get_level() == 0) {
		GridNode::add_difs(a, b);
	}
}

Real Binary::next_omega(Real* f, Real *df) const {
	int i0, j0, k0;
	Real df0 = 0.0;
	Real f0 = 0.0;
	Real dx = get_dx();
	Real d1, d2;
	const Grid& g = *this;
	if (f == NULL) {
		f = &f0;
		df = &df0;
	}
	const Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		i0 = i + BW - 1;
		for (int j = 1; j < PNX - 1; j++) {
			j0 = j + BW - 1;
			for (int k = 1; k < PNX - 1; k++) {
				k0 = k + BW - 1;
				if (!zone_is_refined(i0, j0, k0) && xc(i) < 0.0) {
					d1 = -pxc(i) * Binary::Omega * Binary::Omega;
					d2 = d1;
					d2 += (g(i0 + 1, j0, k0).enthalpy() - g(i0 - 1, j0, k0).enthalpy()) / dx / 2.0;
					d2 += (phi(i + 1, j, k) - phi(i - 1, j, k)) / dx / 2.0;
					*df += 2.0 * d1 / Binary::Omega * dx3 * g(i0, j0, k0).rho();
					*f += d2 * dx3 * g(i0, j0, k0).rho();
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<const Binary*>(get_child(i))->next_omega(f, df);
		}
	}
	return fabs(Binary::Omega - *f / *df);
}

Real Binary::sum_dlz() const {
	int i0, j0, k0;
	Real s = 0.0;
	Real dx = get_dx();
	Real rho, x, y;
	const Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		i0 = i + BW - 1;
		for (int j = 1; j < PNX - 1; j++) {
			j0 = j + BW - 1;
#pragma omp parallel for private(k0,rho,x,y) reduction(+:s)
			for (int k = 1; k < PNX - 1; k++) {
				k0 = k + BW - 1;
				if (!zone_is_refined(i0, j0, k0)) {
					x = pxc(i);
					y = pyc(j);
					rho = (*this)(i0, j0, k0).rho();
					s += (x * (fy(i, j, k) + fy(i, j + 1, k))) * rho * dx3 * 0.5;
					s -= (y * (fx(i, j, k) + fx(i + 1, j, k))) * rho * dx3 * 0.5;
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			s += static_cast<const Binary*>(get_child(i))->sum_dlz();
		}
	}
	return s;
}

Real Binary::domega() const {
	int i0, j0, k0;
	Real domega = 0.0;
	Real dx = get_dx();
	Real dsum;
	int l;
	const Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		i0 = i + BW - 1;
		for (int j = 1; j < PNX - 1; j++) {
			j0 = j + BW - 1;
#pragma omp parallel for private(k0,dsum,l) reduction(+:domega)
			for (int k = 1; k < PNX - 1; k++) {
				k0 = k + BW - 1;
				if (!zone_is_refined(i0, j0, k0)) {
					dsum = (*this)(i0, j0, k0).rho() - (*this).U0(i0, j0, k0).rho();
					if (pxc(i) > 0.0) {
						domega += pyc(j) * dsum * dx3;
					} else {
						domega -= pyc(j) * dsum * dx3;
					}
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			domega += static_cast<const Binary*>(get_child(i))->domega();
		}
	}
	if (get_level() == 0) {
		Vector<Real, 4> m1, m2;
		m1 = mass_sum(0);
		m2 = mass_sum(1);
		domega /= (m1[1] - m2[1]);
	}
	return domega;
}

Real Binary::sum_lz() const {
	int i0, j0, k0;
	Real s = 0.0;
	Real dx = get_dx();
	Real x, y;
	const Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		i0 = i + BW - 1;
		for (int j = 1; j < PNX - 1; j++) {
			j0 = j + BW - 1;
#pragma omp parallel for private(k0,x,y) reduction(+:s)
			for (int k = 1; k < PNX - 1; k++) {
				k0 = k + BW - 1;
				if (!zone_is_refined(i0, j0, k0)) {
					s += (*this)(i0, j0, k0).lz() * dx3;
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			s += static_cast<const Binary*>(get_child(i))->sum_lz();
		}
	}
	return s;
}

Vector<Real, 4> Binary::find_l1(const Vector<Real, 4>& m1, const Vector<Real, 4>& m2, Vector<Real, 4>* v) const {
	Vector<Real, 4> v0 = -1.0e99;
	if (v == NULL) {
		v = &v0;
	}
	int i0, j0, k0;
	Real x, y, O2, dx;
	Real fxm, fxp, fyp, fym;
	O2 = Binary::Omega * Binary::Omega;
	dx = get_dx();
	Real phixp, phixm, phiyp, phiym, phi0;
	for (int i = 1; i < PNX - 1; i++) {
		i0 = i + BW - 1;
		for (int j = 1; j < PNX - 1; j++) {
			j0 = j + BW - 1;
			for (int k = 1; k < PNX - 1; k++) {
				x = pxc(i);
				y = pyc(j);
				if (fabs(pzc(k)) > this->get_dx()) {
					continue;
				}
				Real dx0 = m2[1] - m1[1];
				Real dx1 = x - m1[1];
				Real dx2 = x - m2[1];
				Real dy0 = m2[2] - m1[2];
				Real dy1 = y - m1[2];
				Real dy2 = y - m2[2];
				Real a2 = dx0 * dx0 + dy0 * dy0;
				Real b2 = dx1 * dx1 + dy1 * dy1;
				Real c2 = dx2 * dx2 + dy2 * dy2;
				Real z2 = 0.5 * ((b2 + c2) - (0.5 / a2) * pow(b2 - c2, 2) - 0.5 * a2);
				if (z2 > (3.0 * get_dx() * get_dx() / 4.0) || z2 < 0.0) {
					continue;
				}
				if (dx1 * dx2 + dy1 * dy2 > 0.0) {
					continue;
				}
				k0 = k + BW - 1;
				phi0 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;
				if (!this->poisson_zone_is_refined(i, j, k)) {
					if (phi0 > (*v)[0]) {
						(*v)[0] = phi0;
						(*v)[1] = pxc(i);
						(*v)[2] = pyc(j);
						(*v)[3] = pzc(k);
					}

				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			*v = static_cast<const Binary*>(get_child(i))->find_l1(m1, m2, v);
		}
	}
	return *v;
}

void Binary::mark_lobes(Real l1, Real l1x, Real l1y, Real c1x, Real c2x, Real c1y, Real c2y) {
	if (get_level() == 0) {
		find_phimin(&c1x, &c1y, 1);
		find_phimin(&c2x, &c2y, 2);
		//	printf("%e %e %e %e\n", c1x, c1y, c2x, c2y);
	}
	Real x, y, z, phi0, nx, ny, nz, phip, phim, dphi, frac;
	int i0, j0, k0;
	Real O2 = Binary::Omega * Binary::Omega;
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				i0 = i + BW - 1;
				j0 = j + BW - 1;
				k0 = k + BW - 1;
				x = pxc(i);
				y = pyc(j);
				z = pzc(k);
				phi0 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;

				lobe(i0, j0, k0) = 0;
				if (phi0 < l1) {
					nx = 0.5 * (fx(i, j, k) + fx(i + 1, j, k)) + x * O2;
					ny = 0.5 * (fy(i, j, k) + fy(i, j + 1, k)) + y * O2;
					nz = 0.5 * (fz(i, j, k) + fz(i, j, k + 1));
					if ((x - l1x) * (c1x - c2x) + (y - l1y) * (c1y - c2y) > 0.0) {
						if (nx * (x - c1x) + ny * (y - c1y) + nz * z < 0.0) {
							lobe(i0, j0, k0) = 1;
						}
					} else {
						if (nx * (x - c2x) + ny * (y - c2y) + nz * z < 0.0) {
							lobe(i0, j0, k0) = 2;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->mark_lobes(l1, l1x, l1y, c1x, c2x, c1y, c2y);
		}
	}
}

Real Binary::find_K(Real phi0, int ln) const {
	Real m2;
	if (ln == 2) {
		m2 = Binary::M1 * Binary::q;
	} else {
		m2 = Binary::M1;
	}
#ifdef Z_REFLECT
	m2 /= 2.0;
#endif

	Real x, y, nx, ny, nz, phi1;
	Real O2 = Binary::Omega * Binary::Omega;
	Real K = 0.0;
	int i0, j0, k0;
	Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				i0 = i + BW - 1;
				j0 = j + BW - 1;
				k0 = k + BW - 1;
				if (lobe(i0, j0, k0) == ln && !poisson_zone_is_refined(i, j, k)) {
					x = pxc(i);
					y = pyc(j);
					phi1 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;
					if (phi1 <= phi0) {
						K += pow(phi0 - phi1, Binary::n);
					}
				}
			}
		}
	}
	K *= dx3;
	K /= pow(Binary::n + 1.0, Binary::n) * m2;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			K += static_cast<const Binary*>(get_child(i))->find_K(phi0, ln);
		}
	}
	if (get_level() == 0) {
		K = pow(K, 1.0 / Binary::n);
	}
//	printf("%e\n", K);
	return K;
}

Real Binary::find_A(Real phi0, int ln, Real* dmda) const {
	Real m2;
	Real dmda0, df;
	if (get_level() == 0) {
		dmda = &dmda0;
		dmda0 = 0.0;
	}
	if (ln == 2) {
		m2 = Binary::M1 * Binary::q;
	} else {
		m2 = Binary::M1;
	}
#ifdef Z_REFLECT
	m2 /= 2.0;
#endif
	df = 0.0;
	Real x, y, nx, ny, nz, phi1;
	Real O2 = Binary::Omega * Binary::Omega;
	Real K = 0.0, H, A, y0;
	int i0, j0, k0;
	Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				i0 = i + BW - 1;
				j0 = j + BW - 1;
				k0 = k + BW - 1;
				if (lobe(i0, j0, k0) == ln && !poisson_zone_is_refined(i, j, k)) {
					x = pxc(i);
					y = pyc(j);
					phi1 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;
					if (ln == 1) {
						H = Binary::phi0_1 - phi1;
						A = Binary::K1;
						y0 = H / A / 8.0;
						if (y0 > 0.0) {
							df += -3.0 * y0 * (y0 + 1.0) / A * sqrt(y0 * y0 + 2.0 * y0);
							K += pow(y0 * y0 + 2.0 * y0, 1.5);
						}
					} else {
						H = Binary::phi0_2 - phi1;
						A = Binary::K2;
						y0 = H / A / 8.0;
						if (y0 > 0.0) {
							//printf( "%e\n", df);
							df += -3.0 * y0 * (y0 + 1.0) / A * sqrt(y0 * y0 + 2.0 * y0);
							K += pow(y0 * y0 + 2.0 * y0, 1.5);
						}
					}
				}
			}
		}
	}
	K *= dx3;
	df *= dx3;
	*dmda += df;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			K += static_cast<const Binary*>(get_child(i))->find_A(phi0, ln, dmda);
		}
	}
	if (get_level() == 0) {
		Real dm = K - m2;
		K = dm / dmda0;
		//	printf("%e %e\n", dm, dmda0);
////		abort();
	}
	return K;
}

Real Binary::find_phimax(Real rho_cut, int ln) const {
	Real phi_max = -1.0e+99;
	Real x, y, nx, ny, nz, phi1;
	Real O2 = Binary::Omega * Binary::Omega;
	Real K = 0.0;
	int i0, j0, k0;
	Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				i0 = i + BW - 1;
				j0 = j + BW - 1;
				k0 = k + BW - 1;
				if (lobe(i0, j0, k0) == ln && !poisson_zone_is_refined(i, j, k)
						&& (*this)(i0, j0, k0).rho() > rho_cut) {
					x = pxc(i);
					y = pyc(j);
					phi1 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;
					phi_max = max(phi_max, phi1);
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			phi_max = max(phi_max, static_cast<const Binary*>(get_child(i))->find_phimax(rho_cut, ln));
		}
	}
	return phi_max;
}

Real Binary::find_phimin(Real* x0, Real* y0, int ln, Real phi_min) const {
	Real x00, y00;
	if (x0 == NULL) {
		x0 = &x00;
		x00 = 0.0;
	}
	if (y0 == NULL) {
		y0 = &y00;
		y00 = 0.0;
	}
	Real phipx, phimx, phipy, phimy;
	Real x, y, nx, ny, nz, phi1;
	Real phi_px, phi_mx, phi_py, phi_my, phi_pz, phi_mz;
	Real O2 = Binary::Omega * Binary::Omega;
	Real K = 0.0;
	int i0, j0, k0;
	Real dx3 = get_dx() * get_dx() * get_dx();
	for (int i = 1; i < PNX - 1; i++) {
		for (int j = 1; j < PNX - 1; j++) {
			for (int k = 1; k < PNX - 1; k++) {
				i0 = i + BW - 1;
				j0 = j + BW - 1;
				k0 = k + BW - 1;
				x = pxc(i);
				y = pyc(j);
				if (((ln == 1 && x > 0.0) || (ln == 2 && x < 0.0)) && !poisson_zone_is_refined(i, j, k)) {
					phi1 = phi(i, j, k) - 0.5 * (x * x + y * y) * O2;
					phi_px = phi(i + 1, j, k) - 0.5 * (pow(pxc(i + 1), 2) + pow(pyc(j), 2)) * O2;
					phi_mx = phi(i - 1, j, k) - 0.5 * (pow(pxc(i - 1), 2) + pow(pyc(j), 2)) * O2;
					phi_py = phi(i, j + 1, k) - 0.5 * (pow(pxc(i), 2) + pow(pyc(j + 1), 2)) * O2;
					phi_my = phi(i, j - 1, k) - 0.5 * (pow(pxc(i), 2) + pow(pyc(j - 1), 2)) * O2;
					phi_pz = phi(i, j, k + 1) - 0.5 * (pow(pxc(i), 2) + pow(pyc(j), 2)) * O2;
					phi_mz = phi(i, j, k - 1) - 0.5 * (pow(pxc(i), 2) + pow(pyc(j), 2)) * O2;
					if (phi_min > phi1 && phi1 <= phi_px && phi1 <= phi_py && phi1 <= phi_pz && phi1 <= phi_mx
							&& phi1 <= phi_my && phi1 <= phi_mz) {
						phi_min = phi1;
						*x0 = x;
						*y0 = y;
					}
				}
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			phi_min = static_cast<const Binary*>(get_child(i))->find_phimin(x0, y0, ln, phi_min);
		}
	}
	return phi_min;
}

Real Binary::pot_sum() const {
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
			a0 += static_cast<const Binary*>(get_child(i))->pot_sum();
		}
	}
	return a0;
}

void Binary::M1M2data(binary_integrals_t* b) {
	if (get_level() == 0) {
		Vector<Real, 4> mf1 = mass_sum(0);
		Vector<Real, 4> mf2 = mass_sum(1);
		Vector<Real, 4> L1 = find_l1(mf1, mf2);
		mark_lobes(L1[0], L1[1], L1[2], mf1[1], mf2[1], mf1[2], mf2[2]);
		b->jc = b->mc = b->m1 = b->m2 = b->js1 = b->js2 = 0.0;
		b->x1 = 0.0;
		b->xdot1 = 0.0;
		b->x2 = 0.0;
		b->xdot2 = 0.0;
		b->total_energy = 0.0;
		b->kinetic = 0.0;
		b->I1 = b->I2 = b->Ic = 0.0;
		b->V1 = b->V2 = 0.0;
	}
	Real h3 = pow(get_dx(), 3);
	const Indexer2d indexer(BW, GNX - BW - 1, BW, GNX - BW - 1);
	int k, j, i;
	State s;
	_3Vec x;
	binary_integrals_t b0;
	b0.mc = b0.jc = b0.m1 = b0.m2 = b0.js1 = b0.js2 = b0.total_energy = b0.kinetic = 0.0;
	b0.x1 = 0.0;
	b0.xdot1 = 0.0;
	b0.x2 = 0.0;
	b0.xdot2 = 0.0;
	b0.I1 = b0.I2 = b0.Ic = b0.V1 = b0.V2 = 0.0;
	int l;
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = BW; i < GNX - BW; i++) {
			if (!zone_is_refined(i, j, k)) {
				s = (*this)(i, j, k);
				x = X(i, j, k);
				if (lobe(i, j, k) == 1) {
					b0.V1 += h3;
					b0.m1 += s.rho() * h3;
					b0.x1 += x * s.rho() * h3;
					b0.xdot1 += s.V(x) * s.rho() * h3;
					b0.xdot1[0] -= Binary::Omega * x[1] * s.rho() * h3;
					b0.xdot1[1] += Binary::Omega * x[0] * s.rho() * h3;
					b0.js1 += s.lz() * h3;
					b0.I1 += (x[0] * x[0] + x[1] * x[1]) * s.rho() * h3;
				} else if (lobe(i, j, k) == 2) {
					b0.V2 += h3;
					b0.m2 += s.rho() * h3;
					b0.x2 += x * s.rho() * h3;
					b0.xdot2 += s.V(x) * s.rho() * h3;
					b0.xdot2[0] -= Binary::Omega * x[1] * s.rho() * h3;
					b0.xdot2[1] += Binary::Omega * x[0] * s.rho() * h3;
					b0.js2 += s.lz() * h3;
					b0.I2 += (x[0] * x[0] + x[1] * x[1]) * s.rho() * h3;
				}
				b0.mc += s.rho() * h3;
				b0.jc += s.lz() * h3;
				b0.total_energy += s.et_inertial(x) * h3 + s.pot_inertial(x) * h3 * 0.5;
				b0.kinetic += s.et_inertial(x) * h3;
				b0.Ic += (x[0] * x[0] + x[1] * x[1]) * s.rho() * h3;
			}
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->M1M2data(b);
		}
	}
	b->Ic += b0.Ic;
	b->I1 += b0.I1;
	b->I2 += b0.I2;
	b->V1 += b0.V1;
	b->V2 += b0.V2;
	b->jc += b0.jc;
	b->mc += b0.mc;
	b->js1 += b0.js1;
	b->js2 += b0.js2;
	b->m1 += b0.m1;
	b->m2 += b0.m2;
	b->x1 += b0.x1;
	b->xdot1 += b0.xdot1;
	b->x2 += b0.x2;
	b->xdot2 += b0.xdot2;
	b->kinetic += b0.kinetic;
	b->total_energy += b0.total_energy;
	if (get_level() == 0) {
#ifdef Z_REFLECT
		b->x1[2] = b->xdot1[2] = 0.0;
		b->x2[2] = b->xdot2[2] = 0.0;
#endif
		Real j1, j2;
		b->q = b->m2 / b->m1;
		b->x1 /= b->m1;
		b->x2 /= b->m2;
		b->xdot1 /= b->m1;
		b->xdot2 /= b->m2;
		b->a = sqrt((b->x1 - b->x2).dot(b->x1 - b->x2));
		j1 = b->m1 * (b->x1[0] * b->xdot1[1] - b->x1[1] * b->xdot1[0]);
		j2 = b->m2 * (b->x2[0] * b->xdot2[1] - b->x2[1] * b->xdot2[0]);
		b->jorb = j1 + j2;
		b->js1 -= j1;
		b->js2 -= j2;
		b->I1 -= b->m1 * b->x1.dot(b->x1);
		b->I2 -= b->m2 * b->x2.dot(b->x2);
	}
}

void Binary::pot_to_grid() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i, j0, k0, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,j0,k0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		j0 = j - 1 + BW;
		k0 = k - 1 + BW;
		for (i = 1, i0 = BW; i < PNX - 1; i++, i0++) {
			const Real d = (*this)(i0, j0, k0).rho();
			const Real R2 = pxc(i) * pxc(i) + pyc(j) * pyc(j);
#ifdef SCF_CODE
			Real pot = d * (phi(i, j, k) - 0.5 * R2 * Binary::Omega * Binary::Omega);
#else
			Real pot = d * (phi(i, j, k) - 0.5 * R2 * Binary::Omega0 * Binary::Omega0);
#endif
			(*this)(i0, j0, k0)[State::pot_index] = pot;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->pot_to_grid();
		}
	}
}

void Binary::inject_from_parent(ChildIndex c) {
	Poisson::inject_from_parent(c);
	pot_from_grid();
}

void Binary::pot_from_grid() {
	const Indexer2d indexer(1, PNX - 2, 1, PNX - 2);
	int j, k, i, j0, k0, i0;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(j,k,i,j0,k0,i0)
	for (int index = 0; index <= indexer.max_index(); index++) {
		j = indexer.x(index);
		k = indexer.y(index);
		j0 = j - 1 + BW;
		k0 = k - 1 + BW;
		for (i = 1, i0 = BW; i < PNX - 1; i++, i0++) {
			const Real d = (*this)(i0, j0, k0).rho();
			const Real R2 = pxc(i) * pxc(i) + pyc(j) * pyc(j);
			Real ph = ((*this)(i0, j0, k0).pot() / d + 0.5 * R2 * Binary::Omega0 * Binary::Omega0);
			phi(i, j, k) = ph;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->pot_from_grid();
		}
	}
	if (get_level() == 0) {
		const int l1 = max_level();
		for (int l = 0; l <= l1; l++) {
			enforce_phi_boundaries(l);
		}
	}
}

void Binary::add_pot_et() {
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
			e += (*this)(i0, j0, k0).pot() - 0.5 * phi(i, j, k) * (*this)(i0, j0, k0).rho();
			(*this)(i0, j0, k0)[State::et_index] = e;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->add_pot_et();
		}
	}
}

Binary* Binary::new_octnode() const {
	return new Binary;
}

void Binary::sub_pot_et() {
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
			e -= (*this)(i0, j0, k0).pot() - 0.5 * phi(i, j, k) * (*this)(i0, j0, k0).rho();
			(*this)(i0, j0, k0)[State::et_index] = e;
		}
	}
	for (i = 0; i < OCT_NCHILD; i++) {
		if (get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->sub_pot_et();
		}
	}
}

void Binary::compute_physical_boundaries(Binary* root) {
	for (int f = 0; f < OCT_NSIB; f++) {
		if (is_phys_bound(f)) {
			static_cast<PoissonPhysBound*>(get_sibling(f))->compute_boundaries(root);
		}
	}
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->compute_physical_boundaries(root);
		}
	}

}

Real Binary::compute_phi(Real x, Real y, Real z) const {
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
			sum += static_cast<const Binary*>(get_child(i))->compute_phi(x, y, z);
		}
	}
	return sum;
}

void Binary::compute_source_sums() {
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
			static_cast<Binary*>(get_child(i))->compute_source_sums();
		}
	}
}

Real Binary::max_dt() const {
#ifdef SCF_CODE
	return Binary::scf_dt;
#else
	return GridNode::max_dt();
#endif
}

Binary::Binary() :
		Poisson() {
	// TODO Auto-generated constructor stub

}

Binary::~Binary() {
	// TODO Auto-generated destructor stub
}

_3Vec Binary::V(int i, int j, int k) const {
	_3Vec v;
	v[0] = -yc(j) * Binary::Omega;
	v[1] = +xc(i) * Binary::Omega;
	v[2] = 0.0;
	return v;
}

_3Vec Binary::Vfx(int i, int j, int k) const {
	_3Vec v;
	v[0] = -yc(j) * Binary::Omega;
	v[1] = +xf(i) * Binary::Omega;
	v[2] = 0.0;
	return v;
}

_3Vec Binary::Vfy(int i, int j, int k) const {
	_3Vec v;
	v[0] = -yf(j) * Binary::Omega;
	v[1] = +xc(i) * Binary::Omega;
	v[2] = 0.0;
	return v;
}

_3Vec Binary::Vfz(int i, int j, int k) const {
	return V(i, j, k);
}

void Binary::add_point(const Vector<Real, 3>& loc, const State& st) {
	Vector<Real, 3> tmp;
	Vector<int, 3> index, ci;
	tmp = loc - this->X(0, 0, 0);
	for (int i = 0; i < 3; i++) {
		index[i] = int(tmp[i] / get_dx() + 0.5);
	}
	tmp = -tmp;
	for (int i = 0; i < 3; i++) {
		tmp[i] += index[i] * get_dx();
	}
	Real d = sqrt(tmp.dot(tmp));
	if (d < get_dx() * 1.0e-10) {
//printf("(%e %e %e) ", loc[0], loc[1], loc[2]);
//	printf("(%e %e %e %e)\n", tmp[0] / get_dx(), tmp[1] / get_dx(), tmp[2] / get_dx(), get_dx());
		(*this)(index) = st;
	} else {
		if (this->get_level() >= 5) {
			printf("(%e %e %e %e)\n", tmp[0] / get_dx(), tmp[1] / get_dx(), tmp[2] / get_dx(), get_dx());
			abort();
		}
		ChildIndex c;
		ci = index / (GNX / 2);
		c.set_index(ci[0], ci[1], ci[2]);
		if (this->get_child(int(c)) == NULL) {
			this->create_child(c);
		}
		static_cast<Binary*>(this->get_child(c))->add_point(loc, st);
	}
}

void Binary::input() {
	DBfile* db = DBOpen( "X.silo", DB_PDB, DB_READ );
	DBucdmesh* mesh;
	DBucdvar* vars[STATE_NF];
	Vector<Real,
3>max_x	, min_x, origin;
	mesh = DBGetUcdmesh(db, "mesh");
	for (int i = 0; i < 3; i++) {
		min_x[i] = reinterpret_cast<Real*>(mesh->min_extents)[i];
		max_x[i] = reinterpret_cast<Real*>(mesh->max_extents)[i];
	}
	origin = max_x + min_x;
	origin *= -0.5;
#ifdef Z_REFLECT
	origin[2] = 0.0;
#endif
	this->set_origin(origin);
	for (int i = 0; i < STATE_NF; i++) {
		vars[i] = DBGetUcdvar(db, output_field_names(i));
	}
	printf("%e %e %e\n", origin[0], origin[1], origin[2]);
	OReal* coords[3];
	Real x0, x1, y0, y1;
	coords[0] = reinterpret_cast<OReal*>(mesh->coords[0]);
	coords[1] = reinterpret_cast<OReal*>(mesh->coords[1]);
	coords[2] = reinterpret_cast<OReal*>(mesh->coords[2]);
	x0 = coords[0][mesh->zones->nodelist[0]];
	x1 = coords[0][mesh->zones->nodelist[1]];
	y0 = coords[1][mesh->zones->nodelist[0]];
	y1 = coords[1][mesh->zones->nodelist[1]];
	Binary::phase = atan2(y1 - y0, x1 - x0);
	if (Binary::phase < 0.0) {
		Binary::phase += 2.0 * M_PI;
	}
//	printf( "%e %e %e %e %e \n",Binary::phase, x0, x1, y0, y1 );
	for (int i = 0; i < mesh->nnodes; i++) {
		Real x, y, p;
		p = Binary::phase;
		x = coords[0][i];
		y = coords[1][i];
		coords[0][i] = cos(p) * x + sin(p) * y;
		coords[1][i] = cos(p) * y - sin(p) * x;
	}
//	printf( "%i %i %i\n", mesh->nnodes, mesh->zones->nzones, mesh->zones->lnodelist);
	for (int i = 0; i < mesh->zones->lnodelist / 8; i++) {
		Vector<Real, 3> x = 0.0;
		for (int j = 0; j < 8; j++) {
			for (int k = 0; k < 3; k++) {
				x[k] += coords[k][mesh->zones->nodelist[8 * i + j]];
			}
		}
		x *= 0.125;
//printf( "%e %e %e\n", x[0 ]- this->get_origin()[0], x[1], x[2]);
//abort();
		State st;
		for (int j = 0; j < STATE_NF; j++) {
			st[j] = ((Real*) vars[j]->vals[0])[i];
		}
		this->add_point(x, st);
	}
	Binary::Omega = 2.38840683;
	Binary::Omega0 = 2.39110996;
	this->set_time(mesh->time * 2.0 * M_PI / Binary::Omega0);
	this->set_time_normal(2.0 * M_PI / Binary::Omega0);
	printf("time = %e\n", get_time());
	for (int i = 0; i < STATE_NF; i++) {
		DBFreeUcdvar(vars[i]);
	}
	DBFreeUcdmesh(mesh);
	DBClose(db);
	printf("Done reading Silo\n");
}

void Binary::compute_x_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->compute_x_flux();
		}
	}
	Grid::compute_x_flux();
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	State u;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,u)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw + 1; i++) {
			u = F(i, j, k);
			u.adjust_rho_fluxes();
			F(i, j, k) = u;
		}
	}
}

void Binary::compute_y_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->compute_y_flux();
		}
	}
	Grid::compute_y_flux();
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	State u;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,u)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw + 1; i++) {
			u = F(j, i, k);
			u.adjust_rho_fluxes();
			F(j, i, k) = u;
		}
	}
}

void Binary::compute_z_flux() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (this->get_child(i) != NULL) {
			static_cast<Binary*>(get_child(i))->compute_z_flux();
		}
	}
	Grid::compute_z_flux();
	const Indexer2d indexer(bw, GNX - bw - 1, bw, GNX - bw - 1);
	int k, j, i;
	State u;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(k,j,i,u)
	for (int index = 0; index <= indexer.max_index(); index++) {
		k = indexer.y(index);
		j = indexer.x(index);
		for (i = bw; i < GNX - bw + 1; i++) {
			u = F(j, k, i);
			u.adjust_rho_fluxes();
			F(j, k, i) = u;
		}
	}
}

#endif
