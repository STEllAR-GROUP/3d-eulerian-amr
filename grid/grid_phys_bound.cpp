#include "grid_phys_bound.h"
#include "grid_node.h"
#include "indexer_2d.h"
#include <stdlib.h>

GridPhysBound::GridPhysBound() :
	GridInterface(), OctNode() {
	reflect = false;
}

OctNode* GridPhysBound::get_sibling(int) const {
	printf("Error 3\n");
	abort();
	return NULL;
}

void GridPhysBound::create_child(const ChildIndex&) {
	printf("Error 0\n");
	abort();
}

void GridPhysBound::destroy_child(const ChildIndex&) {
	printf("Error 1\n");
	abort();
}

State GridPhysBound::operator()(int a, int b, int c) const {
	const Vector<int, 3> v = map(a, b, c);
	State s = sibling->operator()(v[0], v[1], v[2]);
	if (!reflect) {
		_3Vec X;
		if (face == XU) {
			X = sibling->Xfx(v[0], v[1], v[2]);
		} else if (face == XL) {
			X = sibling->Xfx(v[0] + 1, v[1], v[2]);
		} else if (face == YU) {
			X = sibling->Xfy(v[0], v[1], v[2]);
		} else if (face == YL) {
			X = sibling->Xfy(v[0], v[1] + 1, v[2]);
		} else if (face == ZU) {
			X = sibling->Xfz(v[0], v[1], v[2]);
		} else {
			X = sibling->Xfz(v[0], v[1], v[2] + 1);
		}
		s.enforce_outflow(X,face);
	} else {
		if (direction == 3) {
			s.reflect_on_z();
		} else {
			assert(false);
		}
	}
	return s;
}

bool GridPhysBound::is_real() const {
	return false;
}

void GridPhysBound::set_sibling(OctFace f, OctNode* sib) {
	face = f;
	OctNode::set_sibling(f, sib);
	const Indexer2d indexer(0, GNX - 1, 0, GNX - 1);
	sibling = static_cast<const GridNode*> (sib);
	assert( sib->get_sibling(invert(f)) == this );
	Vector<int, 3> v;
	int i;
	switch (f) {
	case XU:
#ifdef X_REFLECT
		reflect = true;
#endif
		direction = 1;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = bw;
			v[1] = indexer.x(index);
			v[2] = indexer.y(index);
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[0] = GNX - i - 1;
				}
				map(i, v[1], v[2]) = v;
			}
		}
		break;
	case XL:
		direction = 1;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = GNX - bw - 1;
			v[1] = indexer.x(index);
			v[2] = indexer.y(index);
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[0] = GNX - i - 1;
				}
				map(i, v[1], v[2]) = v;
			}
		}
		break;
	case YU:
#ifdef Y_REFLECT
		reflect = true;
#endif
		direction = 2;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = indexer.x(index);
			v[1] = bw;
			v[2] = indexer.y(index);
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[1] = GNX - i - 1;
				}
				map(v[0], i, v[2]) = v;
			}
		}
		break;
	case YL:
		direction = 2;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = indexer.x(index);
			v[1] = GNX - bw - 1;
			v[2] = indexer.y(index);
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[1] = GNX - i - 1;
				}
				map(v[0], i, v[2]) = v;
			}
		}
		break;
	case ZU:
#ifdef Z_REFLECT
		reflect = true;
#endif
		direction = 3;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = indexer.x(index);
			v[1] = indexer.y(index);
			v[2] = bw;
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[2] = GNX - i - 1;
				}
				map(v[0], v[1], i) = v;
			}
		}
		break;
	case ZL:
		direction = 3;
#pragma omp parallel for schedule(OMP_SCHEDULE) private(v,i)
		for (int index = 0; index <= indexer.max_index(); index++) {
			v[0] = indexer.x(index);
			v[1] = indexer.y(index);
			v[2] = GNX - bw - 1;
			for (i = 0; i < GNX; i++) {
				if (reflect) {
					v[2] = GNX - i - 1;
				}
				map(v[0], v[1], i) = v;
			}
		}
		break;
	}
}

GridPhysBound::~GridPhysBound() {
}
