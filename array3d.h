#ifndef ARRAY3D_H_
#define ARRAY3D_H_

#include "real.h"
#include "vector.h"
#include <stdio.h>

template<class T, int NX, int NY, int NZ>
class Array3d {
protected:
	T data[NZ][NY][NX];
public:
	Array3d() {
		return;
	}
	Array3d(const Array3d<T, NX, NY, NZ>& a) {
		*this = a;
	}
	Array3d<T, NX, NY, NZ>& operator=(const Array3d<T, NX, NY, NZ>& a) {
		if (&a != this) {
			for (int k = 0; k < NZ; k++) {
				for (int j = 0; j < NY; j++) {
					for (int i = 0; i < NX; i++) {
						data[k][j][i] = a.data[k][j][i];
					}
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator/=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] /= a.data[k][j][i];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator*=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] *= a.data[k][j][i];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator+=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] += a.data[k][j][i];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator-=(const Array3d<T, NX, NY, NZ>& a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] -= a.data[k][j][i];
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator+=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] += a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator-=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] -= a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] = a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator*=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] *= a;
				}
			}
		}
		return *this;
	}
	Array3d<T, NX, NY, NZ>& operator/=(const T a) {
		for (int k = 0; k < NZ; k++) {
			for (int j = 0; j < NY; j++) {
				for (int i = 0; i < NX; i++) {
					data[k][j][i] /= a;
				}
			}
		}
		return *this;
	}
	virtual T& operator()(int i, int j, int k) {
		assert( i < NX );
		assert( j < NY );
		assert( k < NZ );
		assert( i >= 0 );
		assert( j >= 0 );
		assert( k >= 0 );
		return data[k][j][i];
	}
	virtual const T operator()(int i, int j, int k) const {
		assert( i < NX );
		assert( j < NY );
		assert( k < NZ );
		assert( i >= 0 );
		assert( j >= 0 );
		assert( k >= 0 );
		return data[k][j][i];
	}
	virtual T divergence(int i, int j, int k) const {
		assert( i < NX-1 );
		assert( j < NY-1 );
		assert( k < NZ-1 );
		assert( i >= 1 );
		assert( j >= 1 );
		assert( k >= 1 );
		return (data[k][j][i + 1] + data[k][j][i - 1] + data[k][j + 1][i]
				+ data[k][j - 1][i] + data[k + 1][j][i] + data[k - 1][j][i])
				- (data[k][j][i] * 6.0);
	}
	virtual T oct_avg(int i, int j, int k) const {
		assert( i < NX-1 );
		assert( j < NY-1 );
		assert( k < NZ-1 );
		assert( i >= 1 );
		assert( j >= 1 );
		assert( k >= 1 );
		return (data[k][j][i] + data[k][j][i + 1] + data[k][j + 1][i]
				+ data[k][j + 1][i + 1] + data[k + 1][j][i] + data[k + 1][j][i
				+ 1] + data[k + 1][j + 1][i] + data[k + 1][j + 1][i + 1])
				* 0.125;
	}
	virtual ~Array3d() {
		return;
	}
};

#endif /* ARRAY3D_H_ */
