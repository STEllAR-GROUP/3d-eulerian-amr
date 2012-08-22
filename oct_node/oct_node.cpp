#include "oct_node.h"
#include <stdlib.h>
#include <stdio.h>


ChildIndex OctNode::my_child_index() const {
	assert( get_parent() != NULL );
	ChildIndex c;
	c.set_x(get_location()[0] % 2);
	c.set_y(get_location()[1] % 2);
	c.set_z(get_location()[2] % 2);
	return c;
}

const OctNode* OctNode::get_root() const {
	int level = get_level();
	const OctNode* ptr = this;
	while( level != 0 ) {
		ptr = ptr->get_parent();
		level = ptr->get_level();
	}
	return ptr;
}

bool OctNode::is_finest() const {
	bool rc = true;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (children[i] != NULL) {
			rc = false;
			break;
		}
	}
	return rc;
}

OctNode::OctNode() {
	parent = NULL;
	level = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		children[i] = NULL;
	}
	for (int i = 0; i < OCT_NSIB; i++) {
		siblings[i] = NULL;
	}
	location = 0;
}


int OctNode::get_node_cnt() const {
	int cnt = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (children[i] != NULL) {
			cnt += children[i]->get_node_cnt();
		}
	}
	return cnt + 1;
}

int OctNode::nchildren() const {
	int cnt = 0;
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (children[i] != NULL) {
			cnt++;
		}
	}
	return cnt;
}

int OctNode::max_level() const {
	int maxl;
	int n;
	if (nchildren() == 0) {
		maxl = get_level();
	} else {
		maxl = 0;
		for (int i = 0; i < OCT_NCHILD; i++) {
			if (children[i] != NULL) {
				n = children[i]->max_level();
				maxl = n > maxl ? n : maxl;
			}
		}
	}
	return maxl;
}

bool OctNode::is_phys_bound(OctFace f) const {
	bool rc;
	switch (f) {
	case XL:
		rc = location[0] == 0;
		break;
	case XU:
		rc = location[0] == ((1 << level) - 1);
		break;
	case YL:
		rc = location[1] == 0;
		break;
	case YU:
		rc = location[1] == ((1 << level) - 1);
		break;
	case ZL:
		rc = location[2] == 0;
		break;
	case ZU:
		rc = location[2] == ((1 << level) - 1);
		break;
	}
	return rc;
}

void OctNode::set_sibling(OctFace f, OctNode* s) {
	if (siblings[f] != NULL) {
		if (!siblings[f]->is_real()) {
			delete siblings[f];
		}
	}
	siblings[f] = s;
}

int OctNode::get_level() const {
	return level;
}

OctNode* OctNode::new_octnode() const {
	return new OctNode;
}

bool OctNode::is_real() const {
	return true;
}

void OctNode::set_parent(OctNode* a) {
	parent = a;
	level = parent->get_level() + 1;
}

OctNode* OctNode::get_child(const ChildIndex& a) {
	return children[a];
}

const OctNode* OctNode::get_child(const ChildIndex& a) const {
	return children[a];
}

ChildIndex OctNode::get_index() const {
	ChildIndex c;
	c.set_index(location[0] % 2, location[1] % 2, location[2] % 2);
	return c;
}

OctNode* OctNode::get_sibling(int a) const {
	return siblings[a];
}

OctNode* OctNode::get_parent() const {
	return parent;
}

const Vector<int, 3>& OctNode::get_location() const {
	return location;
}

void OctNode::output_tree() const {
	printf("%i %i %i %i\n", level, location[0], location[1], location[2]);
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (children[i] != NULL) {
			children[i]->output_tree();
		}
	}
}

void OctNode::create_child(const ChildIndex& c) {
	assert( children[c] == NULL );
	ChildIndex s;
	children[c] = this->new_octnode();
	children[c]->set_parent(this);
	children[c]->location = location * 2 + c.vector();
	OctFace f1;
	OctFace f2;
	s = c;
	if (c.get_x() == 0) {
		s.set_x(1);
		f1 = XL;
		f2 = XU;
	} else {
		s.set_x(0);
		f2 = XL;
		f1 = XU;
	}
	if (children[s] != NULL) {
		children[c]->set_sibling(f2, children[s]);
		children[s]->set_sibling(f1, children[c]);
	}
	if (siblings[f1] != NULL) {
		children[c]->set_sibling(f1, siblings[f1]->children[s]);
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, children[c]);
		}
	}
	s = c;
	if (c.get_y() == 0) {
		s.set_y(1);
		f1 = YL;
		f2 = YU;
	} else {
		s.set_y(0);
		f2 = YL;
		f1 = YU;
	}
	if (children[s] != NULL) {
		children[c]->set_sibling(f2, children[s]);
		children[s]->set_sibling(f1, children[c]);
	}
	if (siblings[f1] != NULL) {
		children[c]->set_sibling(f1, siblings[f1]->children[s]);
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, children[c]);
		}
	}
	s = c;
	if (c.get_z() == 0) {
		s.set_z(1);
		f1 = ZL;
		f2 = ZU;
	} else {
		s.set_z(0);
		f2 = ZL;
		f1 = ZU;
	}
	if (children[s] != NULL) {
		children[c]->set_sibling(f2, children[s]);
		children[s]->set_sibling(f1, children[c]);
	}
	if (siblings[f1] != NULL) {
		children[c]->set_sibling(f1, siblings[f1]->children[s]);
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, children[c]);
		}
	}
}

void OctNode::destroy_child(const ChildIndex& c) {
	assert( children[c] != NULL );
	ChildIndex s;
	OctFace f1;
	OctFace f2;
	s = c;
	if (c.get_x() == 0) {
		s.set_x(1);
		f1 = XL;
		f2 = XU;
	} else {
		s.set_x(0);
		f2 = XL;
		f1 = XU;
	}
	if (children[s] != NULL) {
		children[s]->set_sibling(f1, NULL);
	}
	if (siblings[f1] != NULL) {
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, NULL);
		}
	}
	s = c;
	if (c.get_y() == 0) {
		s.set_y(1);
		f1 = YL;
		f2 = YU;
	} else {
		s.set_y(0);
		f2 = YL;
		f1 = YU;
	}
	if (children[s] != NULL) {
		children[s]->set_sibling(f1, NULL);
	}
	if (siblings[f1] != NULL) {
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, NULL);
		}
	}
	s = c;
	if (c.get_z() == 0) {
		s.set_z(1);
		f1 = ZL;
		f2 = ZU;
	} else {
		s.set_z(0);
		f2 = ZL;
		f1 = ZU;
	}
	if (children[s] != NULL) {
		children[s]->set_sibling(f1, NULL);
	}
	if (siblings[f1] != NULL) {
		if (siblings[f1]->children[s] != NULL) {
			(siblings[f1]->children[s])->set_sibling(f2, NULL);
		}
	}
	delete children[c];
	children[c] = NULL;
}

OctNode::~OctNode() {
	for (int i = 0; i < OCT_NCHILD; i++) {
		if (children[i] != NULL) {
			destroy_child(i);
		}
	}
}

