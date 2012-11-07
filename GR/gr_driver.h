#ifndef GR_DRIVER_H_
#define GR_DRIVER_H_

#include "grid/grid_node.h"

class GRDriver {
private:
	int rk;
	int int_output_freq;
	int ostep_cnt;
	int verbosity;
	Real time_output_freq;
protected:
	GridNode* root;
	int step_cnt;
	virtual void state_sums_out( FILE*, const State& );
public:
	void reset_time();
	virtual void sub_step(Real dt, Real beta,int);
	void set_output_frequency_by_time(Real);
	void set_output_frequency_by_step(int);
	void set_output_off();
	void set_verbosity(int);
	virtual void step(Real dt,bool=false);
	void step();
	void step_to_time(Real);
	void set_rk(int);
	int get_rk() const;
	Real get_time() const;
	GRDriver(GridNode*);
	virtual ~GRDriver();
};


#endif /* HYDRO_DRIVER_H_ */
