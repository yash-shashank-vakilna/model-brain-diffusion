#ifndef RANDOMFOREST_H
#define RANDOMFOREST_H

/*
 *
 * randomForest.h: classification and regression using random forests
 * BRIAN Software Package Version 3.0
 *
 * $Id: randomForest.h 422 2016-10-19 00:27:19Z frithjof $
 *
 * 0.10 (07/07/15): initial version
 * 0.11 (08/07/15): no missing data
 * 0.12 (08/07/15): discrete labels, posterior distribution introduced
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements a random forest classification algorithm.
*/

//! Represents a data sample.

template<typename T>
class sample {
public:
	vecD<T>	fv;				//!< feature vector of this sample
	unsigned int lb;			//!< class label (discrete, input)
	T	wt;				//!< weight of this sample (input)
	unsigned int idx;			//!< sort index (temporary)

	sample(const vecD<T>& f, const unsigned int l, const T w = T(1.0))
		//! allocates a sample from feature vector ft.
		 : fv(f), lb(l), wt(w), idx(std::numeric_limits<unsigned int>::max())
		{ }
	unsigned int nf() const
		{ return fv.N; }
};

using uvec = std::vector<unsigned int>;		//! a vector of unsigned integers

template<typename T>
using vsample = std::vector<sample<T>>;		//! a vector of samples

//! Represents a forest of regression trees.

template<typename T>
class randomForest {
public:
	enum lossfn { sqloss, entropy };
private:
	//! Implements a context to split nodes in a regression tree.
	class splitter {								//! local class, implements a context for splitting nodes in a regression tree.
		const vsample<T>& data;		//!< data sample for splitting decision
		const randomForest<T>& rf;	//!< reference to random forest
		uvecD	inv;			//!< maps index to sample number
		uvecD	dc;			//!< counts occurrence of sample with index i

		uvec	getLocations(const uvecD& sidx) const				//! initialize inv, dc, and loc from sorted indices.
			{ uvec loc; for (unsigned int i = 0; i < sidx.N; i++) {
       		         	const unsigned int z = sidx(i), l = inv(z);
        	        	for (unsigned int j = 0; j < dc(z); j++) loc.push_back(l); };
			  return loc; }
		T	entropyOf(const vecD<T>& v, const T w) const			//! returns entropy of data v
			{ if (w == T(0)) return T(0); T e{0};
    			  for (unsigned int i = 0; i < v.N; i++) { if (isSmall<T>(v(i))) continue;
				const T p = v(i)/w; e -= p*std::log(p); };
			  return e; }
		T	testSqLoss(const uvec& loc, const unsigned int f, unsigned int& fs, T& vs) const;
		T	testEntropy(const uvec& loc, const unsigned int f, unsigned int& fs, T& vs) const;
	public:
		splitter(const vsample<T>& _data, const randomForest<T>& _rf)
			//! allocates a context for splitting data in a regression tree.
			: data(_data), rf(_rf), inv(rf.nt), dc(rf.nt)
			{ inv = std::numeric_limits<unsigned int>::max(); dc = 0;
			  for (unsigned int i = 0; i < data.size(); i++) {
				unsigned int j = data[i].idx; inv[j] = i; dc[j]++; } }
		bool	findSplit(unsigned int& fs, T& vs) const;
	};
	//! Represents a node in a regression tree.
	class node {									//! local class, represents a node in a regression tree.
		node*	child[2];		//!< child nodes to the left and right
		unsigned int fi;		//!< feature id this node splits on
		T	vi;			//!< value this node splits on
		vecD<T>	post;			//!< posterior distribution
		bool	leaf;			//!< true if this node is a leaf

		bool	isPure(const vsample<T>& data) const				//! returns true if all samples data contain the same label.
			{ const unsigned int l = data[0].lb;
		 	  for (const auto& d: data) if (d.wt > T(0) && d.lb != l) return false;
			  return true; }
		unsigned int sideOf(const sample<T>& d) const
			{ return d.fv(fi) <= vi; }
		bool	splitData(std::vector<vsample<T>>& cd, const vsample<T>& data) const
			//! splits samples data into children (left, right) based on feature fi and value vi.
			{ for (const auto d: data) cd[sideOf(d)].push_back(d);		// split data
			  return cd[0].size() != 0 && cd[1].size() != 0; }		// check if node is pure
		void	initLeaf(const vsample<T>& data)
			{ leaf = true; for (const auto d: data) post[d.lb]++; }
		void	clear()
			{ delete child[0]; delete child[1]; }
		void	assign(const node& n)
			{ fi = n.fi; vi = n.vi; post = n.post; leaf = n.leaf;
			  child[0] = n.child[0]? new node(*n.child[0]): nullptr;
			  child[1] = n.child[1]? new node(*n.child[1]): nullptr; }
	public:
		node()
			: fi(0), vi(T(0)), post(), leaf(true)
			{ child[0] = nullptr; child[1] = nullptr; }
		node(const vsample<T>& data, const randomForest& ft, const unsigned int d = 1)
		//! allocates a node in a regression tree rt upto a maximum depth d, based on training samples data.
		 : fi(std::numeric_limits<unsigned int>::max()), vi(T(0)), post(ft.nl), leaf(false)
			{ child[0] = nullptr; child[1] = nullptr; post = 0;		// initialize fields
			  if (d == ft.depth || isPure(data)) { initLeaf(data); return; } // check if leaf node
			  splitter sp(data,ft); std::vector<vsample<T>> cd(2);
			  if (sp.findSplit(fi,vi) == false || splitData(cd,data) == false) { // try to split data
				initLeaf(data); return; };				// not successful, so mark as leaf
			  child[0] = new node(cd[0], ft, d+1); child[1] = new node(cd[1], ft, d+1); } // else set children and recurse
		node(const node& n)							//! copy constructor
			: fi(0), vi(0), post(nullptr), leaf(nullptr)
			{ assign(n); }
		node(node&& n)								//! move constructor
			: fi(n.fi), vi(n.vi), post(n.post), leaf(n.leaf)
			{ child[0] = n.child[0]; n.child[0] = nullptr;
			  child[1] = n.child[1]; n.child[1] = nullptr; }
		~node()	{ clear(); }
		node&	operator=(const node& b)					//! copy assigns from n.
			{ if (this != &b) { clear(); assign(b); }; return *this; }
		node&	operator=(node&& b)						//! move assigns from n.
			{ assert(this != &b);
			  delete child[0]; child[0] = n.child[0]; n.child[0] = nullptr;
			  delete child[1]; child[1] = n.child[1]; n.child[1] = nullptr;
			  fi = n.fi; vi = n.vi; post = std::move(n.post); leaf = n.leaf;
			  return*this; }
		vecD<T>	pdf(const sample<T>& d) const					//! returns posterior label pdf for sample d based on this decision tree.
			{ if (leaf) { const T s = post.sum(); return s? post/s: post; }	// at leaf, return pdf 
			  return child[sideOf(d)]->pdf(d); }				// else recurse
	};
// random forest local data and functions follow.
	const unsigned int depth = 1000;	//!< maximum depth of the resulting tree
	const unsigned int nl;			//!< number of labels (classes)
	const unsigned int nf;			//!< number of features in each sample
	const lossfn loss;			//!< loss function for data splitting
	std::vector<node> trees;		//!< sample indices sorted by each feature
	std::vector<uvecD> sidx;		//!< sample indices sorted by each feature
	unsigned int nt;			//!< number of training samples (for each tree)

	using spair = std::pair<unsigned int, T>;
	static bool sfless(const spair& a, const spair& b) { return a.second < b.second; }
	void	sortByFeatures(vsample<T>& s)						//! sorts training sample by features. Side effect: adds index to samples.
		{ nt = s.size(); sidx.clear();
		  for (unsigned int i = 0; i < nt; i++) s[i].idx = i;
        	  for (unsigned int f = 0; f < nf; f++) { std::vector<spair> sf(nt);	// sample indices and feature f
                	for (unsigned int i = 0; i < nt; i++) sf[i] = spair{i,s[i].fv[f]};
                	std::sort(sf.begin(), sf.end(), sfless); uvecD idx(nt); 	// sort by feature f
                	for (unsigned int i = 0; i < nt; i++) idx[i] = sf[i].first; 	// and keep indices
			sidx.push_back(idx); } }
	node	makeTree(vsample<T>& train, const T fr)					//! allocates a tree based on a random subsample.
		{ if (fr == T(1.0)) { sortByFeatures(train); return node(train,*this); } // use full set
		  const unsigned int nt = train.size(), n = FTOU(nt*fr);
		  uniformDist<unsigned int> ud(0,nt); vsample<T> s;			// get a random sample of the training data
		  for (unsigned int i = 0; i < n; i++) s.push_back(train[ud()]);
		  sortByFeatures(s); return node(s,*this); }
public:
	randomForest(vsample<T>& train, const unsigned int ntr, const unsigned int _nl,
		const T fr = T(1.0), const lossfn _loss = sqloss)
		//! allocates a context for random forest classification, given training data train,
		// number of trees ntr, number of labels nl, subsample factor fr, and loss function loss.
		: nl(_nl), nf(train[0].nf()), loss(_loss), trees(ntr), sidx(), nt(0)
		{ for (unsigned int i = 0; i < ntr; i++) trees[i] = makeTree(train,fr); }
	uvecD	work(const vsample<T>& data) const					//! classifies samples in data.
		{ uvecD preds(data.size()); preds = 0; unsigned int i = 0;		// allocate memory for predictions
		  for (const auto& d: data) { vecD<T> post(nl); post = T(0);		// for all test samples
		  	for (const auto& tr : trees) post += tr.pdf(d);			// sample posterior from all trees
			preds[i] = post.amax(); i++; };					// get label with maximum probability
		  return preds; }
	vecD<T>	pdfOf(const sample<T>& d) const						//! returns PDF for sample d.
		{ vecD<T> post(nl); post = T(0);
		  for (const auto& tr : trees) post += tr.pdf(d);			// sample posterior from all trees
		  post.normalize(); return post; }					// return normalized posterior probabilities
	T	msr(const uvecD& preds, const vsample<T>& data) const			//! given preds and labels, returns mislabeled rate
		{ const unsigned int n = data.size(); if (n == 0) return T(0); T r{0}; 
		  for (unsigned int i = 0; i < n; i++) r += (data[i].lb != preds(i));	// sum mismatch
		  return r/T(n); }
};

template<typename T>
T randomForest<T>::splitter::testSqLoss(const uvec& loc, const unsigned int f, unsigned int& fs, T& vs) const
//! returns criterion for splitting data according to feature f based on squared loss.
{	const unsigned int n = loc.size(); T min{std::numeric_limits<T>::max()};
	T r{0}, s{0}, cl{0}, cr{0}, wl{0}, wr{0};
	for (const auto l : loc) { const T w = data[l].wt, t = T(data[l].lb);		// put all data into right side R
    		r += w*SQR(t); cr += w*t; wr += w; }
	for (unsigned int i = 0; i < n-1; i++) { const auto l = loc[i], l1 = loc[i+1];	// for all features
		const T t = T(data[l].lb), w = data[l].wt;
		s += w*SQR(t); r -= w*SQR(t); r = std::max(r,T(0));			// shift label l from right to left side
		cr -= w*t; cl += w*t; wl += w; wr -= w;
		const T L = std::max(s-SQR(cl)/wl,T(0)), R = std::max(r-SQR(cr)/wr,T(0)); // calculate (left, right) squared loss
		if (data[l].fv(f) == data[l1].fv(f)) continue;				// do not split here if feature is the same as next
		const T I = L+R; if (I < min) { min = I; fs = f;			// check criterion, if we improved...
			vs = T(0.5)*(data[l].fv(f)+data[l1].fv(f)); } };		// ...report value vs and feature fs to calling function
	return min;									// return criterion
}

template<typename T>
T randomForest<T>::splitter::testEntropy(const uvec& loc, const unsigned int f, unsigned int& fs, T& vs) const
//! returns criterion for splitting data according to feature f based on entropy loss.
{	const unsigned int n = loc.size(); T min{std::numeric_limits<T>::max()};
	vecD<T> cl(rf.nl), cr(rf.nl); cl = T(0); cr = T(0); T wr{0}, wl{0};	
	for (const auto l : loc) { const T w = data[l].wt; const auto t = data[l].lb;  	// put all data into right side R
		cr[t] += w; wr += w; }
	for (unsigned int i = 0; i < n-1; i++) { const auto l = loc[i], l1 = loc[i+1];	// for all features
		const auto t = data[l].lb; const T w = data[l].wt; 
		wr -= w; wl += w; cl[t] += w; cr[t] -= w;				// shift label l from right to left side
		if (data[l].fv(f) == data[l1].fv(f)) continue;				// do not split here if feature is the same as next
		const T L = entropyOf(cl,wl), R = entropyOf(cr,wr);			// calculate (left, right) entropy
		const T I = (L*wl+R*wr)/(wl+wr); if (I < min) { min = I; fs = f;	// check criterion, if we improved...
			vs = T(0.5)*(data[l].fv(f)+data[l1].fv(f)); } };		// ...report value vs and feature fs to calling function
	return min;									// return criterion
}

template<typename T>
bool randomForest<T>::splitter::findSplit(unsigned int& fs, T& vs) const
//! finds optimal split according to loss function fn; returns feature fs and value vs.
{	fs = std::numeric_limits<unsigned int>::max(); T min{std::numeric_limits<T>::max()};
	for (unsigned int f = 0; f < rf.nf; f++) {
        	const auto loc = getLocations(rf.sidx[f]); unsigned int fi = 0; T vi{0};
		const T I = rf.loss == lossfn::entropy? testEntropy(loc,f,fi,vi): testSqLoss(loc,f,fi,vi);
		if (I < min) { min = I; fs = fi; vs = vi; } };
	return min != std::numeric_limits<T>::max();
}

#endif
