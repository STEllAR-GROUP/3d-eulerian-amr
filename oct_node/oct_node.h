#ifndef OCOctNode_NODE_H_
#define OCOctNode_NODE_H_

#include "oct_face.h"
#include "child_index.h"

class OctRegistry;

class OctNode {
private:
	OctNode* children[OCT_NCHILD];
	OctNode* siblings[OCT_NSIB];
	OctNode* parent;
	int level;
	Vector<int,3> location;
protected:
	void set_parent(OctNode*);
public:
	const Vector<int,3>& get_location() const;
	bool is_finest() const;
	ChildIndex my_child_index() const;
	void output_tree() const;
	OctNode* get_child(const ChildIndex&);
	const OctNode* get_child(const ChildIndex&) const;
	virtual bool is_real() const;
	int get_node_cnt() const;
	int get_level() const;
	virtual void destroy_child(const ChildIndex&);
	virtual void create_child(const ChildIndex&);
	virtual bool is_phys_bound(OctFace) const;
	OctNode();
	virtual ~OctNode();
	virtual OctNode* new_octnode() const;
	virtual void set_sibling(OctFace, OctNode*);
	int max_level() const;
	int nchildren() const;
	virtual OctNode* get_sibling(int) const;
	OctNode* get_parent() const;
	const OctNode* get_root() const;
	ChildIndex get_index() const;
};

#endif
