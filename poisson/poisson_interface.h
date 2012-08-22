
#ifndef POISSONINTERFACE_H_
#define POISSONINTERFACE_H_

#include "real.h"

class PoissonInterface {
public:
	PoissonInterface() {

	}
	virtual ~PoissonInterface() {

	}
	virtual Real get_phi(int, int, int) const = 0;
	virtual Real get_dphi(int, int, int) const = 0;
};

#endif /* POISSONINTERFACE_H_ */
