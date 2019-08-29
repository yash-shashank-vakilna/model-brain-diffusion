#ifndef KDTREE_H
#define KDTREE_H

/*
 *
 * kdTree.h: kd-tree in 3D
 * BRIAN Software Package Version 3.0
 *
 * $Id: kdTree.h 438 2016-11-10 00:58:28Z frithjof $
 *
 * 0.10 (02/04/16): re-implementation based on nanoflann
 * v406 (28/09/16): bumped to version 3.0
 *
 */

//! Implements a kd-tree

template<typename T> class kdTree {

	//! Represents a search result from a neighbor query.

	struct kdResult {
		const bool nearest;		//!< true if a single neighbor is sought for.
		T	dmax;			//!< best distance.
		std::vector<unsigned int> ind;	//!< list of items found.

		kdResult(const bool nr, const T d = std::numeric_limits<T>::max())
			//! Allocates a context for collecting results in a nearest neighbor search.
			: nearest(nr), dmax(d), ind()
			{ }
		void	add(const T d, const unsigned int i)				//! adds item i at distance d to result list.
			{ if (nearest) dmax = d; ind.push_back(i); }
	};

	//! Represents a node in a kd-tree.

	class kdNode {									//! virtual base class for nodes
	public:
		virtual ~kdNode() { }
		virtual bool isLeaf() const = 0;					// distinguishes between leafs and internal nodes
		virtual kdNode* clone() const = 0;					// virtual constructor for copying nodes
		virtual kdNode* move() const = 0;					// virtual constructor for moving nodes
	};

	//! Represents a leaf in a kd-tree.

	struct kdLeaf : public kdNode {							//! represents leaf nodes
		const unsigned int left;	//!< index of left child node.
		const unsigned int right;	//!< index of right child node.

		kdLeaf(const unsigned int _left = 0, const unsigned int _right = 0)	//! Allocates a node in a kd-tree with children left and right.
		 : kdNode(), left(_left), right(_right) { }
		kdLeaf(const kdLeaf& k)							//! copy-constructs node from node k.
		 : left(k.left), right(k.right) { }
		kdLeaf(kdLeaf&& k)							//! move-constructs node from node k.
		 : left(k.left), right(k.right) { }
		bool	isLeaf() const							//! returns true if this is a leaf node.
			{ return true; }
		kdLeaf& operator=(const kdLeaf& k)					//! copy assigns from kdNode b.
			{ left = k.left; right = k.right; return *this; }
		kdLeaf* clone() const							//! returns a deep copy of this node.
			{ return new kdLeaf(*this); };
		kdLeaf* move() const							//! returns a move-copy of this node.
			{ return new kdLeaf(std::move(*this)); };
	};

	//! Represents an internal node in a kd-tree.

	struct kdInt : public kdNode {							//! represents internal nodes
		const unsigned int feat;	//!< index of distinguishing feature.
		const T	low;			//!< low feature limit.
		const T	high;			//!< high feature limit.
		const kdNode* c1;		//!< pointer to left child node.
		const kdNode* c2;		//!< pointer to right child node.

		kdInt(const unsigned int _feat = 0, const T _low = T(0), const T _high = T(0),
			const kdNode* _c1 = nullptr, const kdNode* _c2 = nullptr)
		: feat(_feat), low(_low), high(_high), c1(_c1), c2(_c2) { }
		kdInt(const kdInt& k)
		 : kdNode(), feat(k.feat), low(k.low), high(k.high), c1(nullptr), c2(nullptr)
		{ if (k.c1) c1 = k.c1->clone(); if (k.c2) c2 = k.c2->clone(); }
		kdInt(kdInt&& k)
		 : kdNode(), feat(k.feat), low(k.low), high(k.high), c1(nullptr), c2(nullptr)
		{ if (k.c1) c1 = k.c1->clone(); 
		  if (k.c2) c2 = k.c2->clone();
		  delete k.c1; delete k.c2; }
		~kdInt()								//! node destructor
			{ delete c1; delete c2; }
		bool	isLeaf() const
			{ return false; }
		kdInt& operator=(const kdInt& k)					//! copy assigns from kdNode b.
			{ delete c1; c1 = k.c1? k.c1->clone(): nullptr;
			  delete c2; c2 = k.c2? k.c2->clone(): nullptr;
			  feat = k.feat; low = k.low; high = k.high; return *this; }
		kdInt* clone() const
			{ return new kdInt(*this); };
		kdInt* move() const
			{ return new kdInt(std::move(*this)); };
	};

	//! Helper structure for intervals in kd-tree allocation.

	struct Interval {
		T	low;			//!< low feature limit.
		T	high;			//!< high feature limit.
		Interval(const T l = std::numeric_limits<T>::max(),
				const T h = std::numeric_limits<T>::lowest())
			//! Allocates a structure representing the interval between scalars l and h.
			: low(l), high(h)
			{ }
		void	extend(const T v)						//! extends interval by scalar v.
			{ low = std::min(low,v); high = std::max(high,v); }
		Interval join(const Interval& b)					//! joins this interval with interval b.
			{ return Interval(std::min(low,b.low),std::max(high,b.high)); }
	};

	const unsigned int maxLeaf = 10;	//!< maximum number of leafs in a terminal node.
	const T	eps = T(1e-6);			//!< a small number.
	unsigned int dim;			//!< dimensionality of each data point.
	kdNode* root;				//!< points to the root node of this tree.
	std::vector<Interval> bbox;		//!< bounding box of this tree.
	uvecD	vind;				//!< vector of data indices.
private:
	virtual T getData(const unsigned int i, const unsigned int d) const = 0;
		//! returns component d of data point i. Implemented in subclass.
	virtual T distance(const vecD<T>& p, const unsigned int i) const = 0;
		//! returns squared distance between p and data point i. Implemented in subclass.
	void	computeBB(std::vector<Interval>& box, const unsigned int n) const	//! computes the bounding box over all dimensions n.
		{ for (unsigned int d = 0; d < dim; d++) { Interval iv;
		  	for (unsigned int i = 0; i < n; i++) iv.extend(getData(vind(i),d));
			box[d] = iv; } }
	kdNode* divideTree(const unsigned int l, const unsigned int r,
			std::vector<Interval>& box)					//! divides tree at indices (l,r).
		{ if (r-l <= maxLeaf) {
			for (unsigned int d = 0; d < dim; d++) { Interval iv;
		  		for (unsigned int i = l; i < r; i++) iv.extend(getData(vind[i],d));
				box[d] = iv; };
			return new kdLeaf(l,r); }
		  else { unsigned int i, f; T v;
			middleSplit(&vind[0]+l, r-l, i, f, v, box);
			std::vector<Interval> lbox(box); lbox[f].high = v;
			kdNode* c1 = divideTree(l, l+i, lbox);
			std::vector<Interval> rbox(box); rbox[f].low = v;
			kdNode* c2 = divideTree(l+i, r, rbox);
			for (unsigned int d = 0; d < dim; d++)
				box[d] = lbox[d].join(rbox[d]);
			return new kdInt(f,lbox[f].high,rbox[f].low,c1,c2); } }
	Interval computeIv(const unsigned int* ind, const unsigned int n,
			const unsigned int d)						//! returns interval of dimension d in n data points indexed by ind.
		{ Interval iv; for (unsigned int i = 0; i < n; i++)
			iv.extend(getData(ind[i],d));
		  return iv; }
	Interval planeSplit(unsigned int* ind, const unsigned int n,
			const unsigned int f, const T v)				//! splits data points indexed by ind at feature f and value v.
		{ unsigned int l = 0, r = n-1; Interval iv;
		  for (;;) { while (l <= r && getData(ind[l],f) < v) l++;
			while (r && l <= r && getData(ind[r],f) >= v) r--;
			if (l > r || !r) break;
			std::swap(ind[l],ind[r]); l++; r--; }
		  iv.low = l; r = n-1;
		  for (;;) { while (l <= r && getData(ind[l],f) <= v) l++;
			while (r && l <= r && getData(ind[r],f) > v) r--;
			if (l > r || !r) break;
			std::swap(ind[l],ind[r]); l++; r--; }
		  iv.high = l; return iv; }
	void	middleSplit(unsigned int* ind, const unsigned int n, unsigned int& i,
			unsigned int& f, T& v, const std::vector<Interval>& box)	//! splits data points indexed by ind; returns split feature f and value v.
		{ T bmax = T(0); for (unsigned int d = 0; d < dim; d++)
			bmax = std::max(box[d].high-box[d].low,bmax);
		  T smax = T(-1.0); f = 0;
		  for (unsigned int d = 0; d < dim; d++) {
			if (box[d].high-box[d].low <= T(1.0-eps)*bmax) continue;
			Interval iv = computeIv(ind,n,f); const T s = iv.high-iv.low;
			if (s > smax) { f = d; smax = s; } }
		  const T s = T(0.5)*(box[f].low+box[f].high);
		  Interval iv = computeIv(ind, n, f); v = CLAMP(s,iv.low,iv.high);
		  iv = planeSplit(ind,n,f,v); i = CLAMP(T(0.5*n),iv.low,iv.high); }
	T	computeInitialDistances(const vecD<T>& p, vecD<T>& dists) const		//! computes distances of data point p to all dimensions of bounding box.
		{ T d2 = T(0); dists = T(0);
		  for (unsigned int d = 0; d < dim; d++) {
			if (p(d) < bbox[d].low) dists[d] = SQR(p(d)-bbox[d].low);
			if (p(d) > bbox[d].high) dists[d] = SQR(p(d)-bbox[d].high);
			d2 += dists[d]; }
		  return d2; }
	void	searchLevel(kdResult& res, const vecD<T>& p, const kdNode* node,
			T mind2, vecD<T>& dists) const
		//! searches for nearest neighbors to data point p at node *n; returns results in res.
		{ if (node->isLeaf()) {
			const kdLeaf* lf = reinterpret_cast<const kdLeaf*>(node);
			for (unsigned int i = lf->left; i < lf->right; i++) {
				const T d = distance(p,vind(i));
				if (d <= res.dmax) res.add(d,vind(i)); };
			return; }
		  const kdInt* nd = reinterpret_cast<const kdInt*>(node);
		  const unsigned int f = nd->feat; const T d1 = p(f)-nd->low, d2 = p(f)-nd->high;
		  const kdNode *c1, *c2; T cut;
		  if (d1+d2 < T(0)) { c1 = nd->c1; c2 = nd->c2; cut = SQR(d2); }
		  else { c1 = nd->c2; c2 = nd->c1; cut = SQR(d1); }
		  searchLevel(res, p, c1, mind2, dists);
		  const T dst = dists[f]; mind2 += cut-dst; dists[f] = cut;
		  if (mind2 <= res.dmax) searchLevel(res, p, c2, mind2, dists);
		  dists[f] = dst; }
	void	findNeighbors(kdResult& res, const vecD<T>& p) const
		//! searches for nearest neighbors to data point p; returns results in res.
		{ if (vind.N == 0 || root == nullptr) return;
		  vecD<T> dists(dim); const T d2 = computeInitialDistances(p, dists);
		  searchLevel(res, p, root, d2, dists); }
protected:
	void	build(const uvecD& ind)							//! constructs kd-tree for data indices ind.
		{ vind = ind; const unsigned int n = vind.N;
		  computeBB(bbox,n); root = divideTree(0,n,bbox); }
public:
	kdTree(const unsigned int _dim = 0)						//! Allocates an empty kd-tree of dimension dim.
		: dim(_dim), root(nullptr), bbox(dim), vind()
		{ }
	kdTree(const kdTree& b)								//! copies kd-tree from tree b.
		 : dim(b.dim), root(nullptr), bbox(b.bbox), vind(b.vind)
		{ if (b.root) root = b.root->clone(); }
	kdTree(kdTree&& b)								//! moves kd-tree from tree b.
		 : dim(b.dim), root(b.root), bbox(b.bbox), vind(b.vind)
		{ b.root = nullptr; }
	virtual ~kdTree()								//! de-allocates this kd-tree.
		{ delete root; root = nullptr; }
	kdTree& operator=(const kdTree& b)						//! assigns a tree from b.
		{ if (this != &b) { delete root; root = nullptr;
			if (b.root) root = b.root->clone();
		 	dim = b.dim; bbox = b.bbox; vind = b.vind; };
		  return *this; }
	kdTree& operator=(kdTree&& b)							//! move assigns a tree from b.
		{ assert(this != &b); delete root; root = b.root; b.root = nullptr;
		  dim = b.dim; bbox = std::move(b.bbox); vind = std::move(b.vind);
		  return *this; }
	unsigned int nearest(const vecD<T>& p) const					//! returns index of nearest neighbor to data point p.
		{ kdResult res(true); findNeighbors(res, p);
		  return res.ind.size()? res.ind.front(): 0; }
	std::vector<unsigned int> within(const vecD<T>& p, const T rad) const
		//! returns list of indices of nearest neighbors within radius rad of data point p.
		{ kdResult res(false,rad); findNeighbors(res, p);
		  return res.ind; }
};
#endif

