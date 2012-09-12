#ifndef GRID_NODE_H_
#define GRID_NODE_H_

#include "grid.h"
#include "grid_output.h"
#include "grid_phys_bound.h"
#include "grid_amr_bound.h"
#include "oct_node/oct_node.h"
#include <silo.h>

class GridPhysBound;

class GridNode: public Grid, public OctNode {
private:
	Array3d<State, GNX, GNX, GNX> E0;
	Array3d<Vector<Real, STATE_NF> , GNX, GNX, GNX> E;
	bool refine[OCT_NCHILD];
	int age;
	static Real refine_alpha;
	static Real t_norm;
protected:
	void clear_refine_flags();
	void propagate_refine_flags_up();
	void mirror_refine(const GridNode* mirror,int dir);
	void set_refine_flags();
	void enforce_proper_nesting();
	bool use_refine_flags();
	virtual void inject_from_parent(ChildIndex);
	virtual int nvar_output() const;
	virtual void load_output(grid_output_t* go, int, int, int) const;
	virtual const char* output_field_names(int) const;
public:
	static void set_time_normal( Real );
	static int max_refine_level;
	static Real min_refine_density;
	virtual void write(FILE*) const;
	virtual void read(FILE*);
	virtual void output(grid_output_t*) const;
	void output(const char*, int) const;
	void set_origin(const _3Vec& );
	void debug() const;
	void error_from_parent();
	Vector<State,4> state_sum() const;
	State state_max() const;
	State state_min() const;
	virtual void init();
	void adjust_x_flux();
	void adjust_y_flux();
	void adjust_z_flux();
	void enforce_boundaries();
	void enforce_edge_boundaries();
	bool check_for_refine(bool clear=true);
	bool zone_is_refined(int, int, int) const;
	virtual void create_child(const ChildIndex&);
	virtual void destroy_child(const ChildIndex&);
	GridNode();
	virtual ~GridNode();
	virtual GridNode* new_octnode() const;
	virtual GridAMRBound* new_amr_bound() const;
	virtual GridPhysBound* new_phys_bound() const;
	virtual void add_difs(Real, Real);
	virtual void clear_difs();
	virtual void compute_x_flux();
	virtual void compute_y_flux();
	virtual void compute_z_flux();
	virtual int fine_point_count() const;
	virtual int ngrids(int) const;
	virtual int ngrids() const;
	virtual int refine_flag_cnt() const;
	virtual void store();
	virtual void restore();
	virtual void sum_x_difs();
	virtual void sum_y_difs();
	virtual void sum_z_difs();
	virtual void set_time(Real);
	virtual void set_dx(Real);
	virtual void inject_from_children();
	virtual int get_zone_cnt() const;
	virtual Real max_dt() const;
	virtual void exec_function(void(*f)(GridNode*), int = 0);
};

#endif /* GRID_NODE_H_ */
