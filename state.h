/*
 * state.h
 *
 *  Created on: May 17, 2012
 *      Author: dmarce1
 */

#ifndef STATE_H_
#define STATE_H_

#ifdef SINGLE
#include "single/state.h"
#endif

#ifdef BINARY
#include "binary/state.h"
#endif

#ifdef EULER
#include "euler/state.h"
#endif

#endif /* STATE_H_ */
