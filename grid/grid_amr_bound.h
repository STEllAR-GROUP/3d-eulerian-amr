#ifndef GRID_AMR_BOUND_H_
#define GRID_AMR_BOUND_H_

#include "grid.h"
#include "oct_node/oct_node.h"

class GridAMRBound: public GridInterface, public OctNode {
protected:
	const Grid* parent;
	Vector<int, 3> rel_offset;
	OctFace face;
protected:
	const Grid* get_amr_parent() const;
	Vector<int,3> get_rel_offset() const;
public:
	GridAMRBound();
	virtual bool is_real() const;
	virtual void set_sibling(OctFace, OctNode*);
	virtual State operator()(int, int, int) const;
	virtual OctNode* get_sibling(int) const;
	virtual ~GridAMRBound();
};

#endif 
