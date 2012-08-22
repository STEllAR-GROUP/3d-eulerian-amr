#ifndef GRID_INTERFACE_H_
#define GRID_INTERFACE_H_

#include "defs.h"
#include "../state.h"


class GridInterface {
public:
	GridInterface() {
		return;
	}
	virtual ~GridInterface() {
		return;
	}
	virtual State operator()(int, int, int) const = 0;
protected:
	static int bw;
};

#endif /* GRID_INTERFACE_H_ */
