
#ifndef INDEXER_2d_by2_BY2_H_
#define INDEXER_2d_by2_BY2__

class Indexer2d_by2 {
private:
	const int xstart;
	const int ystart;
	const int xspan;
	const int yspan;
	const int maximum;
public:
	Indexer2d_by2(int,int,int,int);
	int x(int) const;
	int y(int) const;
	int max_index() const;
	virtual ~Indexer2d_by2();
};

#endif /* INDEXER_2d_by2_H_ */
