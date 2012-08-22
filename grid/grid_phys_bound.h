#ifndef GRID_PHYS_BOUND_H_
#define GRID_PHYS_BOUND_H_

#include "grid_node.h"
#include "oct_node/oct_node.h"

class GridNode;

class GridPhysBound: public GridInterface, public OctNode {
private:
	Array3d<Vector<int, 3> , GNX, GNX, GNX> map;
	bool reflect;
	int direction;
	OctFace face;
protected:
	const GridNode* sibling;
public:
	GridPhysBound();
	virtual void create_child(const ChildIndex&);
	virtual void destroy_child(const ChildIndex&);
	virtual bool is_real() const;
	virtual void set_sibling(OctFace, OctNode*);
	virtual State operator()(int, int, int) const;
	virtual OctNode* get_sibling(int) const;
	virtual ~GridPhysBound();
};

#endif /* GRID_PHYS_BOUND_H_ */
