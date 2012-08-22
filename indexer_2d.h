
#ifndef INDEXER_2D_H_
#define INDEXER_2D_H_

class Indexer2d {
private:
	const int xstart;
	const int ystart;
	const int xspan;
	const int yspan;
	const int maximum;
public:
	Indexer2d(int,int,int,int);
	int x(int) const;
	int y(int) const;
	int max_index() const;
	virtual ~Indexer2d();
};

#endif /* INDEXER_2D_H_ */
