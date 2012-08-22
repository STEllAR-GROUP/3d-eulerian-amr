#include "poisson.h"
#include "poisson_amr_bound.h"
/*
 -6.2779017857162E-004	6.2779017857142E-003	1.0463169642857E-003
 6.2779017857148E-003	-6.2779017857143E-002	-1.0463169642857E-002
 1.0463169642854E-003	-1.0463169642858E-002	-1.7438616071427E-003

 7.3242187500005E-003	-7.3242187500001E-002	-1.2207031250000E-002
 -7.3242187500001E-002	7.3242187500000E-001	1.2207031250000E-001
 -1.2207031250000E-002	1.2207031250000E-001	2.0345052083334E-002

 0.0000000000000E+000	0.0	0.0000000000000E+000
 0.0000000000000E+000	2.3809523809524E-001	0.0000000000000E+000
 0.0	0.0	0.0


 */

//#define V2
#define V1

#ifdef V1
const double interp[3][3][3] = { { {-8.239746093750E-004, +8.239746093750E-003, +1.373291015625E-003}, {
			+8.239746093750E-003, -8.239746093750E-002, -1.373291015625E-002}, {+1.373291015625E-003, -1.373291015625E-002,
			-2.288818359375E-003}}, { {+8.239746093750E-003, -8.239746093750E-002, -1.373291015625E-002}, {
			-8.239746093750E-002, +8.239746093750E-001, +1.373291015625E-001}, {-1.373291015625E-002, +1.373291015625E-001,
			+2.288818359375E-002}}, { {+1.373291015625E-003, -1.373291015625E-002, -2.288818359375E-003}, {
			-1.373291015625E-002, +1.373291015625E-001, +2.288818359375E-002}, {-2.288818359375E-003, +2.288818359375E-002,
			+3.814697265625E-003}}};
#endif

PoissonAMRBound::PoissonAMRBound() :
	GridAMRBound(), PoissonInterface() {
}

PoissonAMRBound::~PoissonAMRBound() {
}

Real PoissonAMRBound::get_phi(int i, int j, int k) const {
	const Poisson* parent = static_cast<const Poisson*> (get_amr_parent());
	assert( parent != NULL );
	int i0, j0, k0;
	int i1, j1, k1;
	int l, n, m;
	int i2, j2, k2, dir;
	Real sum;
	Real phi_f;
	i0 = (rel_offset[0] + i);
	j0 = (rel_offset[1] + j);
	k0 = (rel_offset[2] + k);
	i1 = i0 % 2;
	j1 = j0 % 2;
	k1 = k0 % 2;
	i0 /= 2;
	j0 /= 2;
	k0 /= 2;
	sum = 0.0;
	if (i0 > PNX - 2) {
		parent = static_cast<const Poisson*> (parent->get_sibling(XU));
		if( !parent->is_real()) {
			return 0.0;
		}
		i0 -= PNX - 2;
	} else if (i0 < 1) {
		parent = static_cast<const Poisson*> (parent->get_sibling(XL));
		if( !parent->is_real()) {
			return 0.0;
		}
		i0 += PNX - 2;
	}
	if (j0 > PNX - 2) {
		parent = static_cast<const Poisson*> (parent->get_sibling(YU));
		if( !parent->is_real()) {
			return 0.0;
		}
		j0 -= PNX - 2;
	} else if (j0 < 1) {
		parent = static_cast<const Poisson*> (parent->get_sibling(YL));
		if( !parent->is_real()) {
			return 0.0;
		}
		j0 += PNX - 2;
	}
	if (k0 > PNX - 2) {
		parent = static_cast<const Poisson*> (parent->get_sibling(ZU));
		if( !parent->is_real()) {
			return 0.0;
		}
		k0 -= PNX - 2;
	} else if (k0 < 1) {
		parent = static_cast<const Poisson*> (parent->get_sibling(ZL));
		if( !parent->is_real()) {
			return 0.0;
		}
		k0 += PNX - 2;
	}
	assert( parent->is_real());
	for (l = 0; l < 3; l++) {
		for (n = 0; n < 3; n++) {
			for (m = 0; m < 3; m++) {
				i2 = i0 + (l - 1) * (2 * i1 - 1);
				j2 = j0 + (n - 1) * (2 * j1 - 1);
				k2 = k0 + (m - 1) * (2 * k1 - 1);
				sum += interp[l][n][m] * parent->get_phi(i2, j2, k2);
			}
		}
	}
	return sum;
}

Real PoissonAMRBound::get_dphi(int i, int j, int k) const {
	const int i0 = (rel_offset[0] + i) / 2;
	const int j0 = (rel_offset[1] + j) / 2;
	const int k0 = (rel_offset[2] + k) / 2;
	return (static_cast<const Poisson*> (get_amr_parent()))->get_dphi(i0, j0, k0);
}

void PoissonAMRBound::set_sibling(OctFace f, OctNode* sib0) {
	GridAMRBound::set_sibling(f, sib0);
	Vector<int, 3> v(0);
	sib = static_cast<const Poisson*> (sib0);
	const GridNode* parent = static_cast<const GridNode*> (sib->get_parent());
	switch (f) {
	case XU:
		v[0] = -1;
		break;
	case XL:
		v[0] = 1;
		break;
	case YU:
		v[1] = -1;
		break;
	case YL:
		v[1] = 1;
		break;
	case ZU:
		v[2] = -1;
		break;
	case ZL:
		v[2] = 1;
		break;
	}
	Vector<int, 3> co, po;
	co = static_cast<const GridNode*> (sib)->get_offset();
	po = parent->get_offset();
	rel_offset = co + v * (PNX - 2);
	rel_offset -= po * 2;
	rel_offset -= BW - 1;
	face = f;
}
