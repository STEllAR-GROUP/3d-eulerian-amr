#include "grid_amr_bound.h"
#include "grid_node.h"
#include <stdlib.h>

GridAMRBound::GridAMRBound() :
	GridInterface(), OctNode() {
}

Vector<int, 3> GridAMRBound::get_rel_offset() const {
	return rel_offset;
}

bool GridAMRBound::is_real() const {
	return false;
}

const Grid* GridAMRBound::get_amr_parent() const {
	return parent;
}

OctNode* GridAMRBound::get_sibling(int) const {
	printf("Error 4\n");
	abort();
	return NULL;
}

State GridAMRBound::operator()(int i, int j, int k) const {
	assert( i >= bw );
	assert( j >= bw );
	assert( k >= bw );
	assert( i < GNX-bw );
	assert( j < GNX-bw );
	assert( k < GNX-bw );
	int i1, j1, k1;
	int i0, j0, k0;
	Vector<Real, STATE_NF> m;
	State u;
	i1 = (rel_offset[0] + i);
	j1 = (rel_offset[1] + j);
	k1 = (rel_offset[2] + k);
	i0 = i1 % 2;
	j0 = j1 % 2;
	k0 = k1 % 2;
	i1 /= 2;
	j1 /= 2;
	k1 /= 2;
	u = (*parent)(i1, j1, k1);
	switch (face) {
	case XL:
	case XU:
		m = minmod((*parent)(i1 + 1, j1, k1) - u, u - (*parent)(i1 - 1, j1, k1));
		if (i0 == 1) {
			m = -m;
		}
		break;
	case YL:
	case YU:
		m = minmod((*parent)(i1, j1 + 1, k1) - u, u - (*parent)(i1, j1 - 1, k1));
		if (j0 == 1) {
			m = -m;
		}
		break;
	case ZL:
	case ZU:
		m = minmod((*parent)(i1, j1, k1 + 1) - u, u - (*parent)(i1, j1, k1 - 1));
		if (k0 == 1) {
			m = -m;
		}
		break;
	}
	u -= m * 0.25;
/*	assert( u.rho() > 0.0 );
	if (u.rho() < 0.0) {
		printf("%e %e\n", u.rho(), m[0]);
		abort();
	}*/
	return u;
}

void GridAMRBound::set_sibling(OctFace f, OctNode* sib) {
	OctNode::set_sibling(f, sib);
	Vector<int, 3> v(0);
	parent = static_cast<const GridNode*> (sib->get_parent());
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
	rel_offset = (static_cast<const GridNode*> (sib))->get_offset() + v * (GNX - 2 * bw);
	rel_offset -= parent->get_offset() * 2;
	face = f;
}

GridAMRBound::~GridAMRBound() {
}
