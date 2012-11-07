#ifndef GRID_H_
#define GRID_H_

#include "array3d.h"
#include "grid_interface.h"
#include "reconstruct.h"
#include "defs.h"
#include "oct_node/oct_face.h"

#define GRID_NZONES ((GNX-2*BW)*(GNX-2*BW)*(GNX-2*BW))
#define GRID_NNODES ((GNX+1)*(GNX+1)*(GNX+1))

class Grid: public GridInterface {
private:
	static Real h0;
	static Real v0x, v0y, v0z;
	static Reconstruct reconstruct;
	Real h;
	Real time;
	Vector<int, 3> offset;
	State FO0;
	State FO;
	State DFO;
	_3Vec origin;
protected:
	Array3d<State, GNX, GNX, GNX> U;
	Real Vxx(int, int, int) const;
	Real Vyy(int, int, int) const;
	Real Vzz(int, int, int) const;
	Array3d<Vector<Real, STATE_NF> , GNX, GNX, GNX> D;
	Array3d<Vector<Real, STATE_NF> , GNX, GNX, GNX> F;
	Array3d<State, GNX, GNX, GNX> U0;
#ifdef SAVE_RECONSTRUCTIONS
	Array3d<State, GNX, GNX, GNX> Uxf;
	Array3d<State, GNX, GNX, GNX> Uyf;
	Array3d<State, GNX, GNX, GNX> Uzf;
#endif
	void add_to_dif(const Vector<Real, STATE_NF>&, int, int, int);
	void set_offset(const Vector<int, 3>&);
	void set_flux(const Vector<Real, STATE_NF>&, int, int, int);
	const Vector<Real, STATE_NF> get_flux(int, int, int) const;
	void incD(const Vector<Real, STATE_NF>&, int, int, int);
	Vector<Real, STATE_NF> getD(int, int, int) const;
	Vector<Real, STATE_NF>& differential(int, int, int);
	Vector<Real, STATE_NF> differential(int, int, int) const;
public:
	virtual void set_origin( const _3Vec& );
	virtual _3Vec get_origin() const {
		return origin;
	}
	void write(FILE*) const;
	void read(FILE*);
	Grid();
	const Vector<int, 3>& get_offset() const;
	virtual ~Grid();
	virtual State operator()(int, int, int) const;
	virtual State& operator()(int, int, int);
	virtual State operator()(const Vector<int,3>& i) const {
		return (*this)(i[0], i[1], i[2]);
	}
	virtual State& operator()(const Vector<int,3>& i){
		return (*this)(i[0], i[1], i[2]);
	}
	Real get_time() const;
	Real get_dx() const;
	Real xc(int) const;
	Real yc(int) const;
	Real zc(int) const;
	Real xf(int) const;
	Real yf(int) const;
	Real zf(int) const;
	virtual _3Vec X(int,int,int) const;
	virtual _3Vec V(int,int,int) const;
	virtual _3Vec Xfx(int,int,int) const;
	virtual _3Vec Vfx(int,int,int) const;
	virtual _3Vec Xfy(int,int,int) const;
	virtual _3Vec Vfy(int,int,int) const;
	virtual _3Vec Xfz(int,int,int) const;
	virtual _3Vec Vfz(int,int,int) const;
	Real rc(int, int, int) const;
	virtual void add_difs(Real, Real);
	virtual void clear_difs();
	virtual void store();
	virtual void restore();
	virtual void compute_x_flux();
	virtual void compute_y_flux();
	virtual void compute_z_flux();
	virtual void sum_x_difs();
	virtual void sum_y_difs();
	virtual void sum_z_difs();
	virtual void sum_x_flow_off();
	virtual void sum_y_flow_off();
	virtual void sum_z_flow_off();
	virtual void set_time(Real);
	virtual void set_dx(Real);
	virtual Real max_dt() const;
	State flow_off() const;
};

#ifdef USE_POTENTIAL_GRID
#include "grid_potential.h"
typedef GridPotential Grid;
#else
typedef Grid Grid;
#endif

#endif /* GRID_H_ */
