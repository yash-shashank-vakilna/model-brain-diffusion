#ifndef SVM_H
#define SVM_H

/*
 *
 * svm.h: linear and kernel support vector machines 
 * BRIAN Software Package Version 3.0
 *
 * $Id$
 *
 * 0.10 (20/11/10): initial implementation - incomplete
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 * v412 (17/10/16): revised, introduced templates
 * v426 (25/10/16): rewritten
 *
 */

/*! \file
    \brief Definitions and classes for least-squares support vector machines.
*/

//! Base class of all support vector machine kernels.

/* for a list of kernels, refer to:
   http://crsouza.com/2010/03/17/kernel-functions-for-machine-learning-applications/#linear
*/ 

enum class selectLandmark { random, kmeans, kmediods };

template<typename T> class svmKernel {
protected:
	T	c;				//!< constant offset
public:
	svmKernel(const T _c = T(0))
		: c(_c) { }
	virtual ~svmKernel() { }
	virtual void setParams(const vecD<T>& par)					//! sets parameter vector.
		{ c = par(0); }
	virtual vecD<T> getParams() const						//! returns current parameters.
		{ vecD<T> p(1); p[0] = c; return p; }
	virtual T operator()(const vecD<T>& a, const vecD<T>& b) const = 0;		//! computes distance between (a,b).
};

//! Linear kernel.

template<typename T> class linKernel : public svmKernel<T> {
public:
	linKernel(const T c = T(0))
		: svmKernel<T>(c) { }
	T	operator()(const vecD<T>& a, const vecD<T>& b) const			//! computes distance between (a,b).
		{ return dot(a,b)+svmKernel<T>::c; }
};

//! Polynomial kernel.

template<typename T> class polyKernel : public svmKernel<T> {
	T	alpha;				//!< slope alpha.
	T	d;				//!< polynomial degree d.
public:
	polyKernel(const T _c, const T _alpha, const T _d)
		: svmKernel<T>(_c), alpha(_alpha), d(_d) { }
	T	operator()(const vecD<T>& a, const vecD<T>& b) const			//! computes distance between (a,b).
		{ return std::pow(alpha*dot(a,b)+svmKernel<T>::c,d); }
	void	setParams(const vecD<T>& par)						//! sets parameter vector.
		{ svmKernel<T>::c = par(0); alpha = par(1); d = par(2); }
	vecD<T> getParams() const							//! returns current parameters.
		{ vecD<T> p(3); p[0] = svmKernel<T>::c; p[1] = alpha; p[2] = d;
		  return p; }
};

//! Radial basis (Gaussian) function kernel.

template<typename T> class rbfKernel : public svmKernel<T> {
public:
	rbfKernel(const T c)
		: svmKernel<T>(c) { }
	T	operator()(const vecD<T>& a, const vecD<T>& b) const			//! computes distance between (a,b).
		{ return std::exp(-svmKernel<T>::c*norm2(a-b)); }
};

//! Multilayer perceptron (sigmoid) kernel.

template<typename T> class mlpKernel : public svmKernel<T> {
	T	alpha;				//!< slope alpha.
public:
	mlpKernel(const T _c, const T _alpha)
		: svmKernel<T>(_c), alpha(_alpha) { }
	T	operator()(const vecD<T>& a, const vecD<T>& b) const			//! computes distance between (a,b).
		{ return std::tanh(alpha*dot(a,b)+svmKernel<T>::c); }
	void	setParams(const vecD<T>& par)						//! sets parameter vector.
		{ svmKernel<T>::c = par(0); alpha = par(1); }
	vecD<T> getParams() const							//! returns current parameters.
		{ vecD<T> p(2); p[0] = svmKernel<T>::c; p[1] = alpha; return p; }
};

//! Inverse multiquadratic kernel.

template<typename T> class imkKernel : public svmKernel<T> {
public:
	imkKernel(const T c)
		: svmKernel<T>(c) { }
	T	operator()(const vecD<T>& a, const vecD<T>& b) const			//! computes distance between (a,b).
		{ return T(1)/std::sqrt(norm2(a-b)+SQR(svmKernel<T>::c)); }
};

//! Implements support vector machines with kernel support using a low-rank linearization approach.

/*
 * modeled after:
 * Djuric, N., Lan, L., Vucetic, S., & Wang, Z. (2014).
 * BudgetedSVM: A Toolbox for Scalable SVM Approximations.
 * Journal of Machine Learning Research, 14, 3813-3817.
 */

template<typename T> class kernelSVM {
	const unsigned int maxit = 1000;	//!< maximum number of solver iterations.
	const T lambda = T(0.1);		//!< regularization parameter.
	const std::vector<vecD<T>>& data;	//!< points to set of samples.
	const uvecD& lbl;			//!< points to true labels of samples.
	const unsigned int N;			//!< number of samples.
	uniformDist<T> ud;			//!< local random number generator.
	svmKernel<T>& kn;			//!< kernel operator.
	const unsigned int M;			//!< number of landmarks.
	const std::vector<vecD<T>> lm;		//!< set of landmarks.
	matD<T>	W;				//!< model: weight matrix.
	vecD<T>	wt;				//!< model: weight vector.

	unsigned int predict(const vecD<T>& d) const					//! predicts label for observation d.
		{ vecD<T> v(M); for allY(m) v[m] = kn(d,lm[m]);
		  return dot(wt,v*W) > T(0); }
	T	trainWith(const vecD<T>& par)						//! determines performance for kernel parameter vector par, returns false positive rate.
		{ kn.setParams(par); train(); return fpRate(); }
	uvecD	closestLM(const std::vector<vecD<T>>& med) const			//! returns assignment of data to landmark centers.
		{ uvecD lm(N); for allX(i) { T dmin = std::numeric_limits<T>::max();
		  	for allY(m) { const T d = norm(data[i]-med[m]);
				if (d < dmin) { dmin = d; lm[i] = m; } } };
		  return lm; }
	void	kMediods(std::vector<vecD<T>>& med, const unsigned int nit = 100,
			const bool verbose = false) const;
	void	kMeans(std::vector<vecD<T>>& med, const unsigned int nit = 100,
			const bool verbose = false) const;
	std::vector<vecD<T>> selectLandmarks(const selectLandmark lt, const unsigned int M) const;
	vecD<T>	solve(const std::vector<vecD<T>>& e) const;
	void	train();
public:
	kernelSVM(const std::vector<vecD<T>>& _data, const uvecD& _lbl,	svmKernel<T>& _kn, 
		const selectLandmark lt = selectLandmark::kmeans, const unsigned int _M = 100)
		: data(_data), lbl(_lbl), N(data.size()), ud(), kn(_kn), M(_M),
		lm(selectLandmarks(lt,M)), W(), wt() { train(); }
	T	fpRate() const								//! returns the false positive rate for this classifier.
		{ const int def[2] = { -1, 1 }; unsigned int e = 0;	
		  for allX(i) { const bool p = predict(data[i]), t = def[lbl(i)] > 0;
			if (p != t) e++; }
		  return T(e)/T(N); }
	void	optimizeKernel(const vecD<T>& min, const vecD<T>& max,
			const bool verbose = false);
};

template <typename T>
void kernelSVM<T>::kMeans(std::vector<vecD<T>>& med, const unsigned int nit, const bool verbose) const
//! starting from landmarks med, computes k-means clustering of data; returns updated landmarks.
{	for (unsigned int it = 0; it < nit; it++) {
		const uvecD lm = closestLM(med); T d{0};
		for allY(m) { vecD<T> v(med[m].N); v = T(0); unsigned int n = 0;
			for allX(i) { if (lm(i) == m) { v += data[i]; n++; } };
			v = n? v/T(n): data[FTOU(ud()*N)];
			d += norm(med[m]-v); med[m] = v; };
		if (verbose) { printf("[%3d] %.4e\r", it, d); fflush(stdout); };
		if (d == T(0)) break; };
	if (verbose) printf("\n");
}

template <typename T>
void kernelSVM<T>::kMediods(std::vector<vecD<T>>& med, const unsigned int nit, const bool verbose) const
//! starting from landmarks med, computes k-mediods clustering of data; returns updated landmarks.
// see Park & Jun (2009) Expert Systems with Applications 36, 3336-3341.
{	for (unsigned int it = 0; it < nit; it++) {
		const uvecD lm = closestLM(med); unsigned int n = 0;
		for allY(m) { T dmin = std::numeric_limits<T>::max(); unsigned int k = 0;
			for allX(i) { T d{0}; if (lm(i) != m) continue;
				for allX(j) { if (lm(j) != m) continue;
					d += norm2(data[i]-data[j]); };
				if (d < dmin) { dmin = d; k = i; } };
			const T d = norm(med[m]-data[k]);
			if (d > T(0)) { med[m] = data[k]; n++; } };
		if (verbose) { printf("[%3d] %u\r", it, n); fflush(stdout); };
		if (n == 0) break; };
	if (verbose) printf("\n");
}

template <typename T>
std::vector<vecD<T>> kernelSVM<T>::selectLandmarks(const selectLandmark md, const unsigned int M) const
//! returns a set of M landmarks according to selection mode md
{	std::vector<vecD<T>> med(M); for allY(m) med[m] = data[FTOU(ud()*N)];
	switch (md) {
	case selectLandmark::random:
		break;									// use random vectors above.
	case selectLandmark::kmeans:
		kMeans(med); break;
	case selectLandmark::kmediods:
		kMediods(med); break;
	default: throw optException("unknown landmark selection method"); };
	return med;
}

template <typename T>
vecD<T> kernelSVM<T>::solve(const std::vector<vecD<T>>& e) const
//! implements a coordinate descent algorithm for L2-loss SVM dual problems.
// See Algorithm 3 of Hsieh et al., ICML 2008
{	T def[2] = { T(-1), T(1) }, INF = std::numeric_limits<T>::max();
	T PGmaxOld = INF, PGminOld = -INF; uvecD id(N); vecD<T> alpha(N), Q(N);
	for allX(i) { Q[i] = dot(e[i],e[i]); id[i] = i; alpha[i] = T(0); }; 
	unsigned int act = N; vecD<T> w(M); w = T(0);
	for (unsigned int it = 0; it < maxit; it++) {
		T PGminNew = INF, PGmaxNew = -INF;
		for (unsigned int s = 0; s < act; s++)
			std::swap(id[s],id[s+FTOU(ud()*(act-s))]);
		for (unsigned int s = 0; s < act; s++) { const unsigned int i = id[s];
			T l = def[lbl(i)], g = dot(w,e[i])*l-T(1), pg = T(0);
			if (alpha[i] == T(0)) {
				if (g > PGmaxOld) {
					act--; std::swap(id[s],id[act]); s--; continue; }
				else if (g < T(0)) pg = g; }
			else if (alpha[i] == T(1)/lambda) {
				if (g < PGminOld) {
					act--; std::swap(id[s],id[act]); s--; continue; }
				else if (g > T(0)) pg = g; }
			else pg = g;
			PGmaxNew = std::max(PGmaxNew,pg); PGminNew = std::min(PGminNew,pg);
			if (std::abs(pg) > safeMin<T>()) { const T ai = alpha[i];
				alpha[i] = std::min(std::max(alpha[i]-g/Q(i),T(0)),T(1)/lambda);
				w += ((alpha[i]-ai)*l)*e[i]; } };					
		if (PGmaxNew-PGminNew <= T(0.01)) {
			if (act == N) break;
			else { act = N; PGmaxOld = INF; PGminOld = -INF; continue; } };
		PGmaxOld = PGmaxNew; if (PGmaxOld <= T(0)) PGmaxOld = INF;
		PGminOld = PGminNew; if (PGminOld >= T(0)) PGminOld = -INF; };
	return w;
}

template <typename T>
void kernelSVM<T>::train()
//! computes the W matrix and wt vector; done once per training.
{	matD<T> V(M,M); 
	for (unsigned int i = 0; i < M; i++) {
		for (unsigned int j = 0; j < M; j++) {
			const T v = kn(lm[i],lm[j]); V(i,j) = v; V(i,j) = v; } };
	W = invSqrt(V); std::vector<vecD<T>> E; vecD<T> e(M);
	for allX(i) { for allY(j) { e[j] = kn(data[i],lm[j]); }; E.push_back(e*W); };
	wt = solve(E);
}

template <typename T>
void kernelSVM<T>::optimizeKernel(const vecD<T>& min, const vecD<T>& max, const bool verbose)
//! determines optimal parameters for the classification kernel.
{	class optimizer : public powell<T> {
		kernelSVM<T>& l;
		T	goal(const vecD<T>& par) { return l.trainWith(par); }		// minimizes mis-classification
	public:
		optimizer(kernelSVM<T>& _l, const unsigned int d)
		: powell<T>(d), l(_l) { }
	};
	vecD<T> par = (min+max)*T(0.5);
	const T m = optimizer(*this, par.N).optimize(par, min, max);			// optimize parameter set
	kn.setParams(par);
	if (verbose) {  for (unsigned int i = 0; i < par.N; i++) printf("%f ", par[i]);
		printf("\nperf %f\n", T(1)-m); };
}

template<typename T> class linearSVM : public ntrust<T> {
	const unsigned int N;			//!< number of data samples
	const unsigned int D;			//!< size of data point
	std::vector<vecD<T>> data;		//!< dense data vectors
	uvecD	tl;				//!< true label
	std::vector<unsigned int> lbl;		//!< label types
	vecD<T>	wt;				//!< weights
	vecD<T>	wi;				//!< initial weight
	uvecD	ss;				//!< subset of observations
	uvecD	ts;				//!< subset of labels
	T	c;				//!< regularization parameter

	T	funcAt(const vecD<T>& x) const						//! returns function value at x.
		{ T f = T(0.5)*dot(x,x);
		  for (unsigned int i = 0; i < ss.N; i++) {
			const T w = dot(x,data[ss(i)]), d = ts(i)? T(1)-w: T(1)+w;
			if (d > T(0)) f += c*SQR(d); }
		  return f; }
	vecD<T>	gradAt(const vecD<T>& x) const						//! returns gradient at x.
		{ vecD<T> g(x.N); g = T(0);
		  for (unsigned int i = 0; i < ss.N; i++) {
			const T w = dot(x,data[ss(i)]);
			const T q = ts(i)? w-T(1): -w-T(1); if (q >= T(0)) continue;
			const T f = ts(i)? c*q: -c*q; g += f*data[ss(i)]; };
		  return g+g+x; }
	uvecD	permute() const								//! returns vector of randomly permuted indices.
		{ uvecD pm(N); for (unsigned int i = 0; i < pm.N; i++) pm[i] = i;
		  for (unsigned int i = 0; i < pm.N; i++) std::swap(pm[i],pm[i+rand()%(N-i)]);
		  return pm; }
	uvecD	partition(const unsigned int np) const					//! returns vector of np partitions of set L.
		{ uvecD st(np+1); for (unsigned int i = 0; i <= np; i++) st[i] = i*N/np;
		  return st; }
	uvecD	getSubset(const uvecD& pm, const unsigned int s, const unsigned int e) const
		{ uvecD ss(pm.N-e+s); unsigned int k = 0;
		  for (unsigned int i = 0; i < s; i++, k++) ss[k] = pm(i);
		  for (unsigned int i = e; i < pm.N; i++, k++) ss[k] = pm(i);
		  return ss; }
	T	accuracy(const uvecD& pd) const						//! returns accuracy of predictions pd.
		{ unsigned int cp = 0;
		  for (unsigned int i = 0; i < N; i++) { if (pd(i) == tl(i)) cp++; };
		  return T(cp)/T(N); }
	T	accuracy(const T c, const unsigned int nf) const			//! returns accuracy from nf rounds of cross-validation.
		{ return accuracy(validate(c,nf)); }
	T	startC() const								//! calculates initial C for parameter selection.
		{ T nm = 0.0; for (unsigned int i = 0; i < N; i++) {
			const T n = norm2(data[i]); if (n > nm) nm = n; }
		  const T c = T(1)/T(2*N*nm);
		  return std::pow(T(2),std::floor(std::log(c)/std::log(T(2)))); }
	vecD<T>	solve(const uvecD& _ss, const uvecD& _ts)				//! computes weights for this classification.
		{ vecD<T> w(D); if (wi.N) w = wi; else w = T(0.01);			// NB must be non-zero
		  ss = _ss; ts = _ts; return ntrust<T>::cg(w); }
	uvecD	validate(const T c, const unsigned int nf);				//! calculates a labeling for cross-validation.
public:
	linearSVM(const std::vector<vecD<T>>& _data, uvecD& _tl)			//! allocates a context for SVM classification.
		: ntrust<T>(T(1e-6)), N(_data.size()), D(_data[0].N), data(_data),
		tl(_tl), lbl(), wt(), wi(), ss(), ts(), c(T(0)) { }
	unsigned int predict(const vecD<T>& x) const					//! predicts label for sample x.
		{ const unsigned int nc = lbl.size(); vecD<T> v(nc); v = T(0);
		  if (nc == 2) return dot(x,wt) > T(0)? lbl[0]: lbl[1];
		  for (unsigned int i = 0; i < x.N; i++) {
			for (unsigned int c = 0; c < nc; c++) v[c] += wt(i*nc+c)*x(i); }
		  return lbl[v.imax()]; }
	void	train(const uvecD& pm, const unsigned int s, const unsigned int e);
	T	optimize(const unsigned int nf,	T cmin, T cmax);
};

template <typename T>
void linearSVM<T>::train(const uvecD& pm, const unsigned int s, const unsigned int e)
//! trains a classifier for this linSVM. NB: sets lbl and wt.
{	uvecD ss = getSubset(pm,s,e);
	std::vector<unsigned int> cnt; uvecD dt(ss.N); lbl.clear();
	for (unsigned int i = 0; i < ss.N; i++) { unsigned int j = 0, t = tl(ss(i));	// find unique labels and their count
		for (j = 0; j < lbl.size(); j++) {
			if (lbl[j] == t) { cnt[j]++; break; } };
		dt[i] = j;
		if (j == lbl.size()) { lbl.push_back(t); cnt.push_back(1); } };
	uvecD st(lbl.size()), si(ss.N); st[0] = 0;					// determine class start points
	for (unsigned i = 1; i < st.N; i++) st[i] = st[i-1]+cnt[i-1];
	for (unsigned int i = 0; i < ss.N; i++) {
		const unsigned int s = st[dt[i]]++; si[s] = ss(i); }
	st[0] = 0; for (unsigned i = 1; i < st.N; i++) st[i] = st[i-1]+cnt[i-1];
	if (lbl.size() == 2) { const unsigned int e = st[0]+cnt[0];			// two-class linSVM
		for (unsigned int i = 0; i < ss.N; i++) dt[i] = i < e? 1: 0;
		wt = solve(si,dt); }
	else {	wt.resize(D*lbl.size());						// multi-class linSVM
		for (unsigned int i = 0; i < st.N; i++) { dt = 0;
			for (unsigned int j = st[i]; j < st[i]+cnt[i]; j++) dt[j] = +1;
			const vecD<T> w = solve(si,dt);
			for (unsigned int j = 0; j < w.N; j++) wt[j*lbl.size()+i] = w(c); } }
}

template <typename T>
uvecD linearSVM<T>::validate(const T cv, const unsigned int nf)
//! returns labeling for a cross-validation with nf groups.
{	c = cv; wi.resize(0); uvecD pm = permute(), st = partition(nf); uvecD pd(N);	// set permutation and partition vector
	for (unsigned int i = 0; i < nf; i++) {	train(pm,st(i),st(i+1));		// train all partitions
		for (unsigned int j = st[i]; j < st[i+1]; j++) {			// predict labels
			const unsigned int k = pm[i]; pd[k] = predict(data[k]); } };
	return pd;
}

template <typename T>
T linearSVM<T>::optimize(const unsigned int nf, T cmin, T cmax)
//! returns optimal classification parameter c, given range from cs to ce.
{	if (cmin <= T(0)) cmin = startC(); T br{0}, bc{0};				// determine minimal c from data
	const uvecD pm = permute(), st = partition(nf);					// set permutation and partition vector
	std::vector<vecD<T>> wo(nf); uvecD pd(N);
	for (c = cmin; c <= cmax; c *= T(2)) { int uc = 0;				// for increasing c
		for (unsigned int i = 0; i < nf; i++) {					// for all partitions
			wi = wo[i]; train(pm,st(i),st(i+1));				// train partition i
			if (wo[i].N && uc >= 0 && norm(wt-wo[i]) > T(1e-6)) uc = -1;	// check weights for changes
			wo[i] = wt;
			for (unsigned int j = st(i); j < st(i+1); j++) {		// predict labels
				unsigned int k = pm(j); pd[k] = predict(data[k]); } }
		const T r = accuracy(pd); if (r > br) { bc = c; br = r; };		// save best accuracy
		printf("c %6.2f r %5.2f bc %.4e\n", std::log2(c), r*100.0, bc); fflush(stdout);
		if (++uc == 3) break; };
	return bc;									// return best c
}

template <typename T>
std::vector<vecD<T>> pca(const std::vector<vecD<T>>& data, const bool scale, const unsigned int nc)
//! performs a principal component analysis of observations data with optional scaling to unit variance; returns projected data.
{	const unsigned int N = ITOU(data.size()), D = data[0].N; assert(N > 0);
	matD<T> X(N,D); for (unsigned int i = 0; i < N; i++) X.setRow(i,data[i]);
	const matD<T> Xc = scale? normalizeCols(X): centerCols(X);
	const matD<T> XtX = (trp(Xc)*Xc)/T(N);
	matD<T> V; vecD<T> v = sev<T>(XtX,V); sortEV<T>(v,V,true);			// perform eigenvalue decomposition with descending sort
	const underscore _; const matD<T> Xt = Xc*V(_,_(0,nc)); std::vector<vecD<T>> tr;
	for (unsigned int i = 0; i < N; i++) tr.push_back(Xt.getRow(i));
	return tr;
}

template <typename T>
std::vector<vecD<T>> kpca(const std::vector<vecD<T>>& data, const svmKernel<T>& kn, const unsigned int nc)
//! performs a kernel PCA of observations data with optional scaling to unit variance; returns projected data.
{	const unsigned int N = data.size(); assert(N > 0); matD<T> K(N,N);
	for (unsigned int i = 0; i < N; i++)
		for(unsigned int j = i; j < N; j++) {
			const T k = kn(data[i],data[j]); K(i,j) = k; K(j,i) = k; };
	matD<T> V; vecD<T> v = sev(K,V); sortEV<double>(v,V,true);			// perform eigenvalue decomposition with descending sort
	const underscore _; const matD<T> Xt = K*V(_,_(0,nc)); std::vector<vecD<T>> tr;
	for (unsigned int i = 0; i < N; i++) tr.push_back(Xt.getRow(i));
	return tr;
}

//! Provides a class for a multi-class linear discriminant analysis of point sets.

template<typename T> class LDA {
	const std::vector<vecD<T>>& data;	//!< reference to data points
	const unsigned int N;			//!< number of data points
	const uvecD& lbl;			//!< cluster membership
	const unsigned int nc;			//!< number of clusters
	std::vector<vecD<T>> cen;		//!< vector of class-wise centroids

	void	computeCentroids()							//! recomputes cluster centers.
		{ vecD<T> cc(data[0].N); cc = T(0); uvecD cnt(nc); cnt = 0;
		  for allClasses(c) cen[c] = cc;					// set centroids to zero
		  for allX(n) {	cen[lbl(n)] += data[n]; cnt[lbl(n)]++; };
		  for allClasses(c) if (cnt(c)) cen[c] /= cnt(c); }			// divide by count
	matD<T> withinScatterMatrix() const						//! returns the sum of the within-class scatter matrices.
		{ const unsigned int D = cen[0].N; matD<T> S(D,D); S = T(0);
		  for allX(n) { const vecD<T> d = data[n]-cen[lbl(n)]; S += outer(d,d); };
		  return S; }
	matD<T> betweenScatterMatrix() const						//! returns the between-class scatter matrix.
		{ const unsigned int D = cen[0].N; vecD<T> m(D); m = T(0);
		  for allClasses(c) { m += cen[c]; }; m /= T(nc);
		  matD<T> S(D,D); S = T(0);
		  for allX(n) { const vecD<T> d = cen[n]-m; S += outer(d,d); };
		  return S; }
	unsigned int closestCentroid(const vecD<T>& v) const				//! returns index of centroid closest to vector v.
		{ T dmin = std::numeric_limits<T>::max(); unsigned int cmin = 0;
		  for allClasses(c) { const T d = norm2(cen[c]-v);
			if (d < dmin) { dmin = d; cmin = c; } };
		  return cmin; }
public:
	LDA(std::vector<vecD<T>>& _data, const uvecD _lbl)				//! allocates a context for kmeans segmentation of data into nc clusters.
		: data(_data), N(data.size()), lbl(_lbl), nc(lbl(lbl.imax())+1), cen(nc)
		{ }
	std::vector<vecD<T>> work(const unsigned int r)					//! performs linear discriminant analysis; returns mapped data.
		{ computeCentroids();
		  const matD<T> Sw = withinScatterMatrix(), Sb = betweenScatterMatrix();
		  const matD<T> K = inv(Sw)*Sb;
		  matD<T> V; vecD<T> v = sev(K,V); sortEV<double>(v,V,true);		// perform eigenvalue decomposition with descending sort
		  const underscore _; const matD<T> Vr = V(_,_(0,r));
		  std::vector<vecD<T>> tr; for allX(n) tr.push_back(data[n]*Vr);
		  return tr; }
};

#endif
