#include "../defs.h"
#ifdef BINARY
#ifndef POI4SSON_DRIVER_H_
#define POI4SSON_DRIVER_H_

#include "../poisson/poisson.h"
#include "binary.h"

class BinaryDriver {
private:
	int int_output_freq;
	int ostep_cnt;
	int verbosity;
	Real time_output_freq;
	Binary* root;
	Real hydrotime, poissontime;
	Real tolerance;
	Real last_dy, last_dx;
	_3Vec last_com;
protected:
	int step_cnt;
public:
	virtual Real step(Real,bool=false);
	BinaryDriver(Binary*);
	virtual ~BinaryDriver();
	virtual void sub_step(Real,Real,int);
	void solve(Real ell_toler=ELL_TOLER);
	void reset_time();
	void set_output_frequency_by_time(Real);
	void set_output_frequency_by_step(int);
	void set_output_off();
	void set_verbosity(int);
	void step();
	void step_to_time(Real);
	int get_rk() const;
	Real get_time() const;
	virtual void state_sums_out( FILE*, const State& );
};

#endif

#endif /* POISSON_DRIVER_H_ */
