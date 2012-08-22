#include "indexer_2d.h"
#include "assert.h"

Indexer2d::Indexer2d(int a, int xstop, int c, int ystop) :
	xstart(a), ystart(c), xspan(xstop - xstart + 1), yspan(ystop - ystart + 1), maximum((ystop - ystart + 1) * (xstop
			- xstart + 1) -1) {
	assert( xspan > 0 );
	assert( yspan > 0 );
}

int Indexer2d::x(int i) const {
	assert( i >= 0);
	assert( i <= maximum);
	return i % xspan + xstart;

}
int Indexer2d::y(int i) const {
	assert( i >= 0);
	assert( i <= maximum);
	return i / xspan + ystart;
}

int Indexer2d::max_index() const {
	return maximum;
}

Indexer2d::~Indexer2d() {
}
