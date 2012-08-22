#include "indexer_2d_by2.h"
#include "assert.h"

Indexer2d_by2::Indexer2d_by2(int a, int xstop, int c, int ystop) :
	xstart(a), ystart(c), xspan(xstop - xstart + 1), yspan(ystop - ystart + 1), maximum((ystop - ystart + 1) * (xstop
			- xstart + 1) / 4 - 1) {
	assert( xspan > 0 );
	assert( yspan > 0 );
}

int Indexer2d_by2::x(int i) const {
	assert( i >= 0);
	assert( i <= maximum);
	return ((2*i) % xspan) + xstart;

}
int Indexer2d_by2::y(int i) const {
	assert( i >= 0);
	assert( i <= maximum);
	return 2*((2*i) / xspan) + ystart;
}

int Indexer2d_by2::max_index() const {
	return maximum;
}

Indexer2d_by2::~Indexer2d_by2() {
}
