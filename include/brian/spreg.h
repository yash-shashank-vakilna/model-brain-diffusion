#ifndef SPREG_H
#define SPREG_H

/*
 *
 * spreg.h: functions for sparse regression problems and independent component analysis
 * BRIAN Software Package Version 3.0
 *
 * $Id: spreg.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (05/02/15): initial version
 * v406 (28/09/16): bumped to version 3.0
 * 
 * [1] B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani (2004). Least Angle Regression. Ann. Statist. 32, 407-499.
 * [2] H. Zou and T. Hastie (2005). Regularization and variable selection via the elastic net. J. Royal Stat. Soc. B. 67, 301-320. 
 * [3] S. Rosset, and Ji Zhu (2007). Piecewise Linear Regularized Solution Paths. Ann. Statist. 35, 1012-1030.
 *
 * loosely modeled after:
 * SpasSM - a Matlab toolbox for performing sparse regression (http://www.imm.dtu.dk/projects/spasm/).
 *
 * [4] A. Hyvaerinen (1999) Fast and Robust Fixed-Point Algorithms for Independent Component Analysis. IEEE Trans. Neur. Netw. 10.
 * [5] T. Lee, M. Girolami and T. Sejnowski (1999). Independent component analysis using an extended infomax algorithm 
 *     for mixed sub-Gaussian and super-Gaussian sources. Neur. Comp. 11, 409-433.
 *
 */

/*! \file
    \brief Provides functions for sparse regression problems.
*/

//! Provides a context for performing sparse linear regression.

template <typename T> 
class spreg {
	const matD<T>& X;			//!< regressor matrix
	const vecD<T>& y;			//!< left-hand side (data)
	const unsigned int n;			//!< number of observations (rows in X)
	const unsigned int p;			//!< number of variables (columns in X)
	std::list<unsigned int> act;		//!< set of active regressors
	std::list<unsigned int> inact;		//!< set of inactive regressors
	matD<T>	G;				//!< Gram matrix Xt*X

	matD<T>	cholinsert(const matD<T>& R, const vecD<T>& x, const matD<T>& X,
			const T delta = T(0)) const					//! returns the Cholesky factorization of [X x]'*[X x].
		{ const T d = dot(x,x)+delta;
		  if (R.M == 0) { matD<T> Q(1,1); Q(0,0) = std::sqrt(d); return Q; }
		  const vecD<T> c = inv(trp(R))*(X.tm(x));
		  vecD<T> r(R.N+1); r = 0; r[r.N-1] = std::sqrt(d-dot(c,c)); 
		  matD<T> Q = R; Q.addCol(c); Q.addRow(r); return Q; }
	matD<T>	choldelete(const matD<T>& R, const unsigned int j) const		//! returns the Cholesky factorization of X'*X where column j of X has been removed.
		{ matD<T> Q = R; Q.delCol(j); const unsigned int n = Q.N;
		  for (unsigned int i = j; i < n; i++) { matD<T> G(2,2);
			if (isSmall<T>(Q(i+1,i))) { Q(i+1,i) = T(0); continue; };
			const T q0 = Q(i,i), q1 = Q(i+1,i), r = std::sqrt(SQR(q0)+SQR(q1));
			G(0,0) = q0/r; G(0,1) = q1/r; G(1,0) = -G(0,1); G(1,1) = G(0,0);
			Q(i,i) = r; Q(i+1,i) = 0; if (i >= n-1) continue;
			const underscore _; matD<T> GQ = G*Q(_(i,i+2),_(i+1,n));
			Q.set(GQ,i,i+1); };
		  Q.delRow(Q.M-1); return Q; }
	T	maxCorrelation(const vecD<T>& r, unsigned int& c) const			//! find variable c with max correlation cmax.
		{ T cmax = T(0); c = 0;
		  for (const unsigned int i: inact) {
			const T ci = std::abs(dot(X.getCol(i),r));
			if (ci > cmax) { cmax = ci; c = i; } };
		  return cmax; };
	matD<T>	Xact() const								//! returns sub-matrix of active variables.
		{ unsigned int c = 0; matD<T> Xa(n,act.size()); 
		  for (const unsigned int i: act) Xa.setCol(c++, X.getCol(i));
		  return Xa; }
	vecD<T>	bols(const matD<T>& Xa, const matD<T>& R) const				//! finds partial OLS solution.
		{ const vecD<T> xy = Xa.tm(y);
		  if (G.M) { matD<T> Ga(act.size(),act.size()); unsigned int ai = 0;
			for (const unsigned int i: act) { unsigned int aj = 0;
				for (const unsigned int j: act) Ga(ai++,aj) = G(i,j);
				aj++; };
			return inv(Ga)*xy; }
		  else return inv(R)*(inv(trp(R))*xy); }
	T	gamma(const vecD<T>& r, const vecD<T>& d, const T cmax) const		//! determines step length.
		{ T gm = T(1.0); if (inact.size() == 0) return gm;
		  for (const unsigned int i: inact) { vecD<T> xi = X.getCol(i); 
			const T cd = dot(xi,d), cr = dot(xi,r);
			const T g0 = (cr-cmax)/(cd-cmax);
			if (g0 > T(0)) gm = std::min(g0,gm);
			const T g1 = (cr+cmax)/(cd+cmax);
			if (g1 > T(0)) gm = std::min(g1,gm); };
		  return gm; }
	bool	updateB(vecD<T>& b, const vecD<T>& bo, const T gm, const T stop) const
		{ const vecD<T> bp = b; unsigned int c = 0;	
		  for (const unsigned int i: act) { b[i] += gm*(bo(c++)-b(i));		// update beta
			if (isSmall<T>(std::abs(b[i]))) b[i] = T(0); };
		  if (stop > 0) { const T t1 = norm1(bp), t2 = norm1(b);		// early stopping at specified bound on L1 norm of beta
			if (t2 >= stop) { const T s = (stop-t1)/(t2-t1);		// interpolation factor 0 < s < 1
				b = bp+s*(b-bp); return true; } }
		  return false; } 
public:
	spreg(const matD<T>& _X, const vecD<T>& _y)
		: X(_X), y(_y), n(X.M), p(X.N), act(), inact(), G()
		{ for (unsigned int i = 0; i < p; i++) inact.push_back(i);
		  if (n/p > 10 && p < 1000) G = trp(X)*X; }
	vecD<T>	lars(const T stop = 0, const bool verbose = false);
	vecD<T>	larsen(const T delta = T(0), T stop = T(0), const bool verbose = false);
};

template <typename T>
vecD<T> spreg<T>::lars(const T stop, const bool verbose)
//! computes least angle regression solution using the LARS algorithm (see Sjostrand, Alg 2)
{	const unsigned int maxVar = std::min(n-1,p);					// maximum number of active variables
	vecD<T> b(p), yh(n); b = T(0); yh = T(0); matD<T> R, Xa;			// set up the LAR coefficient vector
	for (unsigned int it = 0; act.size() < maxVar; it++) {				// step 3
		const vecD<T> r = y-yh; unsigned int c = 0;				// step 4: update residual	
		const T cmax = maxCorrelation(r,c); 					// step 5: find max correlation c
		if (G.M == 0) R = cholinsert(R,X.getCol(c),Xa);
		act.push_back(c); inact.remove(c); 				 	// step 6: move variable c to active set
		if (verbose) printf("[%3d] add %u %zd\n", it, c, act.size());
		Xa = Xact(); const vecD<T> bo = bols(Xa,R), d = Xa*bo-yh;		// step 7, 8: find partial beta_OLS & current direction d
		const T gm = gamma(r,d,cmax); yh += gm*d;				// step 9: get step length gamma & update yh (step 11)
		if (updateB(b,bo,gm,stop)) break;					// step 10: update beta & check convergence
		if (stop < 0 && T(act.size()) >= -stop) break; };			// else stop at specified number of variables
	return b;
}

template <typename T>
vecD<T> spreg<T>::larsen(const T delta, T stop, const bool verbose)
//! computes the LARS-EN algorithm for estimating Elastic Net solutions.
{	const unsigned int maxVar = delta == T(0)? std::min(n,p): p;			// maximum number of active variables
	vecD<T> b(p), mu(n); b = T(0); mu = T(0); matD<T> R, Xa; bool cond = false;	// set up the LASSO coefficient vector
	if (delta > T(0) && stop > T(0)) stop = stop/(T(1.0)+delta);			// correction of stopping criterion to fit naive Elastic Net
	for (unsigned int it = 0; act.size() < maxVar && it < 8*maxVar; it++) {		// LARS main loop
		const vecD<T> r = y-mu; unsigned int c = 0, drop = 0;
		const T cmax = maxCorrelation(r,c);					// find max correlation
		if (cond == false) {
			if (G.M == 0) R = cholinsert(R,X.getCol(c),Xa,delta);
			act.push_back(c); inact.remove(c); 				// add variable
			if (verbose) printf("[%3d] add %u %zd\n", it, c, act.size()); }
		else cond = false;							// if a variable has been dropped, do one step with this configuration
		Xa = Xact(); const vecD<T> bo = bols(Xa,R), d = Xa*bo-mu;		// find partial OLS solution
	  	T gt = std::numeric_limits<T>::max(); c = 0;
		for (const unsigned int i: act) {					// compute length of walk along equiangular direction
			const T g = b(i)/(b(i)-bo(c++)); if (g <= T(0)) continue;
			if (g < gt) { gt = g; drop = i; } };
		T gm = gamma(r,d,cmax); if (gt < gm) { cond = true; gm = gt; };		// check if variable should be dropped
		mu += gm*d; if (updateB(b,bo,gm,stop)) break;				// early stopping at specified number of variables
		if (cond) { if (G.M == 0) { c = 0; 					// if LASSO condition satisfied...
			for (const unsigned int i: act) {
				if (i == drop) { R = choldelete(R,c); break; }; c++; } };
			inact.push_back(drop); act.remove(drop);			// ...move dropped variable to inactive set
			if (verbose) printf("[%3d] drop %u %zd\n", it, drop, act.size()); }
		if (stop < T(0) && T(act.size()) >= -stop) break; } 			// early stopping at specified number of variables
	return b;
}

template <typename T>
vecD<T> lasso(const matD<T>& X, const vecD<T>& y, T stop = T(0), const bool verbose = false)
{	spreg<T> sp(X,y); return sp.larsen(T(0),stop,verbose);				// adjust to avoid T shrinkage
}  

template <typename T>
vecD<T> elasticNet(const matD<T>& X, const vecD<T>& y, const T delta = T(0),
	T stop = T(0), const bool verbose = false)
{	spreg<T> sp(X,y); return sp.larsen(delta,stop,verbose)*T(1.0+delta);		// adjust to avoid T shrinkage
}  

template <typename T>
vecD<T> spca(matD<T>& B, const matD<T>& X, unsigned int K, const T delta = 0,
	const T stop = T(0), const bool verbose = false)
{	const unsigned int n = X.M, p = X.N; K = std::min(K,std::min(p,n-1));		// number of principal components
	const T eps = T(1e-6); const unsigned int maxit = 100; const underscore _;
	matD<T> A; B.resize(p,K); B = T(0);						// setup SPCA matrices A and B
	{ matD<T> U,V; vecD<T> d = svd(X,U,V,true); A = V(_,_(0,K)); };			// use standard PCA as starting condition
	for (unsigned int k = 0; k < K; k++) {						// for each component
		for (unsigned int it = 0; it < maxit; it++) {  
			vecD<T> bo = B.getCol(k), ak = X*A.getCol(k), bk;
			if (delta == std::numeric_limits<T>::max()) {			// soft thresholding, calculate beta directly
				vecD<T> ax = X.tm(ak), aax = ax; bk.resize(ax.N);
				for (unsigned int i = 0; i < ax.N; i++) aax[i] = std::abs(ax[i]);
				if (stop < 0 && -stop < p) {
        				vecD<T> sa = aax; sort(sa,true);		// ascending sort
					const T off = sa(-std::floor(stop));
					for (unsigned int i = 0; i < bk.N; i++) {
						const T b = std::max(T(0),aax(i)-off);
						bk[i] = std::copysign(b,ax(i)); } }
				else { 	for (unsigned int i = 0; i < bk.N; i++) {
						const T b = std::max(T(0),aax(i)-stop);
						bk[i] = std::copysign(b,ax(i)); } } }
			else { spreg<T> sr(X,ak); bk = sr.larsen(delta,stop,verbose); }; // find beta by elastic net regression
			T n = norm(bk); if (n > 0) bk /= n;    				// normalize to Euclidean length 1
			matD<T> Ak = A(_,_(0,k)); vecD<T> t = X.tm(X*bk);		// update A
			ak = t-Ak*(Ak.tm(t)); n = norm(ak); if (n > 0) ak /= n;    	// normalize to Euclidean length 1
			A.setCol(k,ak); B.setCol(k,bk);
			const T d = bo.chisq(bk);
			if (verbose) { printf("%u [%2d] %e\n", k, it, d); }
			if (d < eps) break; } };
	const matD<T> XB = X*B; vecD<T> sd(K); sd = T(0);
	if (K == 1) {
		for (unsigned int i = 0; i < XB.M*XB.N; i++) sd[0] += SQR(XB.x[i]); }
	else { matD<T> R,Q = qr(XB,R);
		for (unsigned int i = 0; i < K; i++) sd[i] = SQR(R(i,i)); };
	sd /= n; return sd;
}

template <typename T>
matD<T> slda(matD<T>& B, const matD<T>& X, const matD<T>& Y, unsigned int K, const T delta = 0,
	const T stop = 0, const bool verbose = false)
//! performs sparse linear disciminant analysis on variables in X and the dummy encoded class-memberships in Y.
{	const unsigned int Q = std::min(K-1,Y.N-1), maxit = 100; const T eps = 1e-6;	// number of classes
	const underscore _; vecD<T> dpi(Y.N); matD<T> Ydpi = Y;
	for (unsigned int j = 0; j < Y.N; j++) { T s = T(0);
		for (unsigned int i = 0; i < Y.M; i++) { s += Y(i,j); }; dpi[j] = s/Y.M; // diagonal "matrix" of class priors
		for (unsigned int i = 0; i < Y.M; i++) Ydpi(i,j) /= dpi[j]; };		// class belongings scaled according to priors
	B.resize(X.N,Q); matD<T> Theta(K,Q); B = T(0); Theta.id();			// coefficients of discriminative directions and optimal scores
	for (unsigned int i = 0; i < Q; i++) {						// for all directions
		for (unsigned int it = 0; it < maxit; it++) {  
			vecD<T> bo = B.getCol(i), bi;
    			{ vecD<T> yi = Y*Theta.getCol(i); spreg<T> sr(X,yi); 
			  bi = sr.larsen(delta,stop,verbose); };			// find beta by elastic net regression
			vecD<T> t = Ydpi.tm(X*bi), u = t; T s = T(0);			// estimate theta
			if (i) { matD<T> Ts = Theta(_,_(0,i)); vecD<T> dt = dpi*t; u -= Ts*(Ts.tm(dt)); }
			for (unsigned int j = 0; j < u.N; j++) s += dpi(j)*SQR(u(j));
			u *= T(1.0/std::sqrt(s)); Theta.setCol(i,u); B.setCol(i,bi);
			s = bo.chisq(bi)/dot(bo,bi);
			if (verbose) { printf("%u [%2d] %e\n", i, it, s); fflush(stdout); }
			if (s < eps) break; } };    
	return Theta;
}

enum nlin { ica_pow3, ica_tanh, ica_gauss, ica_skew };

//! Implements a class for performing independent component analysis.
template <typename T>
class ica {
	const unsigned int maxit = 1000;
	const matD<T>& X;			//!< input data matrix (variables in columns)
	const unsigned int ns;			//!< number of samples
	const unsigned int nt;			//!< number of traces
	unsigned int pc;			//!< number of principal components
	const unsigned int ic;			//!< number of independent components
	const T	eps;				//!< residual tolerance (for convergence criterion)
	matD<T>	W;				//!< whitening matrix
	matD<T>	D;				//!< de-whitening matrix
	const matD<T> Xw;			//!< whitened data

	matD<T>	permute(const matD<T>& Xi) const						//! returns data Xi permuted in time.
		{ matD<T> Xp = Xi; uniformDist<unsigned int> ud(0,Xp.M);
		  for (unsigned int i = 0; i < Xp.M; i++) Xp.setRow(ud(),Xi.getRow(i));
		  return Xp; }
	matD<T>	mpower(const matD<T>& A) const							//! returns square root of matrix A.
		{ matD<T> Q; vecD<T> d = sev(A,Q); matD<T> S(d.N,d.N); S = 0;
		  for (unsigned int i = 0; i < d.N; i++) S(i,i) = T(1.0/std::sqrt(d(i)));
		  for (unsigned int i = 0; i < Q.N; i++) {
			vecD<T> qi = Q.getCol(i); qi /= norm(qi); Q.setCol(i,qi); }
		  return Q*S*trp(Q); };
	vecD<T>	pca(matD<T>& P, const bool verbose = false)					//! computes PCA of data X, returns projector P and loadings ld.
		{ const auto S = sev(inner(X,X)); const T vtot = S.d.sum();			// compute eigen-decomposition
		  pc = std::min(S.d.nnz(),pc); vecD<T> ld(pc); P.resize(S.U.N,pc);		// find non-zero eigenvalues
		  for (unsigned int i = 0; i < pc; i++) { const unsigned int j = S.d.N-i-1;	// sort by decreasing magnitude
			ld[i] = S.d(j); P.setCol(i,S.U.getCol(j)); };
		  if (verbose) printf("retained variance %.2f%%\n", 100.0*ld.sum()/vtot);
		  return ld; }
	void	whiten(const matD<T>& P, const vecD<T>& ld)					//! computes whitening and de-whitening matrix
		{ W.resize(P.M,P.N); D.resize(P.N,P.M);
		  for (unsigned int i = 0; i < P.N; i++) {					// for all principal components
			const vecD<T> p = P.getCol(i); const T d = std::sqrt(ld(i));		// assumes loadings are non-zero
			W.setCol(i,p/d); D.setRow(i,p*d); } }					// set scaled eigenvectors
	matD<T>	project()
		{ normalizeCols(X); matD<T> P; vecD<T> ld = pca(P,true);
		  whiten(P,ld); return X*W; }
	void	updateSymmetric(matD<T>& B, const nlin nl, const T a) const;
	void	updateDeflate(vecD<T>& b, const nlin nl, const T a) const;
	void	learn(matD<T>& B, const matD<T>& Xp, const vecD<T>& signs, const T lr);
	void	setSigns(vecD<T>& signs, const matD<T>& B) const;
public:
	ica(const matD<T>& _X, const unsigned int _pc, const unsigned int _ic)
		: X(_X), ns(X.M), nt(X.N), pc(_pc), ic(_ic), eps(T(1e-6)), W(), D(), Xw(project())
		{ }
	matD<T>	symmetric(const nlin nl = ica_pow3, const T a = T(1.0), const bool verbose = false);
	matD<T>	deflate(const nlin nl = ica_pow3, const T a = T(1.0), const bool verbose = false);
	matD<T>	extended(vecD<T>& signs, const T rate = T(0.0001), const bool anneal = false,
			const bool verbose = false);
	matD<T>	pcaOnly()
		{ return Xw; }
};

template <typename T>
void ica<T>::updateSymmetric(matD<T>& B, const nlin nl, const T a) const
{	matD<T> Xt = trp(Xw), U = Xw*B; vecD<T> s(U.N); s = T(0); T *u = U.x;
	switch (nl) {
	case ica_pow3:
		for (unsigned int i = 0; i < U.M*U.N; i++) u[i] = u[i]*u[i]*u[i];
		B = (Xt*U)/T(ns)-B*T(3.0); break;
	case ica_tanh:
		for (unsigned int j = 0; j < U.N; j++) {
			for (unsigned int i = 0; i < U.M; i++) {
				T t = tanh(a*u[i]); U(i,j) = t; s[j] += T(1.0)-SQR(t); } };
		for (unsigned int j = 0; j < B.N; j++)
			for (unsigned int i = 0; i < B.M; i++) B(i,j) *= a*s[j];
		B = (Xt*U-B)/T(ns); break;
	case ica_gauss:
		for (unsigned int j = 0; j < U.N; j++) {
			for (unsigned int i = 0; i < U.M; i++) {
				T u2 = SQR(U(i,j)), ex = T(std::exp(-0.5*a*u2));
				U(i,j) *= ex; s[j] += ex*T(1.0-a*u2); } };
		for (unsigned int j = 0; j < B.N; j++)
			for (unsigned int i = 0; i < B.M; i++) B(i,j) *= s[j];
		B = (Xt*U-B)/T(ns); break;
	case ica_skew:
		for (unsigned int i = 0; i < U.M*U.N; i++) u[i] = u[i]*u[i];
		B = (Xt*U)/T(ns);
	default: break; };
}

template <typename T>
matD<T> ica<T>::symmetric(const nlin nl, const T a, const bool verbose)
{	matD<T> B(pc,ic), Bo(pc,ic); random(B); Bo = T(0);
	for (unsigned int it = 0; it < maxit; it++) {					// iterate...
		B = B*mpower(inner(B,B)); T bmin = std::numeric_limits<T>::max();
		for (unsigned int i = 0; i < ic; i++) {					// for all independent components
			T bi = dot(B.getCol(i),Bo.getCol(i));
			bmin = std::min(bmin,std::abs(bi)); };
		if (verbose) { printf("[%4d] %.4e\r", it, bmin); fflush(stdout); };
		if (1.0-bmin < eps) break;						// check convergence
		Bo = B; updateSymmetric(B,nl,a); }					// else update from data and go on
	matD<T> A = trp(D)*B, M = W*B; return X*M;
}

template <typename T>
void ica<T>::updateDeflate(vecD<T>& b, const nlin nl, const T a) const
{	vecD<T> u = Xw*b; T s = T(0);
	switch (nl) {
	case ica_pow3:
		for (unsigned int i = 0; i < u.N; i++) u[i] = u[i]*u[i]*u[i];
		b = Xw.tm(u)/Xw.M-b*T(3.0); break;
	case ica_tanh:
		for (unsigned int i = 0; i < u.N; i++) {
			T t = std::tanh(a*u[i]); u[i] = t; s += T(1.0)-SQR(t); };
		b = Xw.tm(u)-b*(a*s); break;
	case ica_gauss:
		for (unsigned int i = 0; i < u.N; i++) {
			T u2 = SQR(u[i]), ex = T(std::exp(-0.5*a*u2));
			u[i] = u[i]*ex; s += ex*(T(1.0)-a*u2); };
		b = Xw.tm(u)-b*s; break;
	case ica_skew:
		for (unsigned int i = 0; i < u.N; i++) u[i] = u[i]*u[i];
		b = Xw.tm(u); break;
	default: break; };
}

template <typename T>
matD<T> ica<T>::deflate(const nlin nl, const T a, const bool verbose)
{	matD<T> B(pc,ic); random(B); const unsigned int flimit = 5; unsigned int nf = 0; 
	for (unsigned int c = 0; c < ic; c++) { T e = T(0);				// for all independent components
		vecD<T> b = B.getCol(c), bo(pc); b /= norm(b); bo = T(0); 		// get column
		for (unsigned int it = 0; it < maxit; it++) {				// iterate...
			b -= B*(B.tm(b)); b /= norm(b); 				// orthonormalize column
			e = std::sqrt(b.chisq(bo)); if (e > T(1.0) || e < eps) break;	// ...until no change: vector converged
			if (verbose) { printf("%u [%4d] %.4e\r", c, it, e); fflush(stdout); };
			bo = b; updateDeflate(b,nl,a); b /= norm(b); };			// else update from data and go on
		if (e < eps) { B.setCol(c,b); nf = 0; }					// converged: save vector
		else if (nf < flimit) { nf++; random(b); B.setCol(c,b); c--; }		// else restart with new random vector
		else break; }			 					// else fail completely					
	matD<T> A = trp(D)*B, M = W*B; return X*M;
}

template <typename T>
void ica<T>::setSigns(vecD<T>& signs, const matD<T>& B) const
//! sets signs from estimated kurtosis.
{	const unsigned int ns = Xw.M, pc = Xw.N; matD<T> S = Xw*trp(B);
	vecD<T> m2(pc), m4(pc); m2 = T(0); m4 = T(0);
	for (unsigned int j = 0; j < pc; j++)						// for all components
		for (unsigned int i = 0; i < ns; i++) {					// for all samples
			T v = SQR(S(i,j)); m2[j] += v; m4[j] += SQR(v); }		// update 2nd and 4th moment
	m2 /= T(ns); m4 /= T(ns);
	for (unsigned int i = 0; i < pc; i++) {						// for all components
		const T v = m4[i]/SQR(m2[i])-T(3.0);					// compute kurtosis and set sign
		signs[i] = v > 0? T(1.0): T(-1.0); }
}

template <typename T>
void ica<T>::learn(matD<T>& B, const matD<T>& Xp, const vecD<T>& signs, const T lr)
//! learns from permuted block Xp.
{	const unsigned int ns = Xw.M, pc = Xw.N, bs = std::min(500u,ns), nb = ns/bs;	// set batch parameters
	for (unsigned int t = 0; t < nb; t++) {						// for all blocks...
		matD<T> U(bs,pc), Ut(bs,pc), Q(pc,pc); U = T(0); Q = T(0);
		for (unsigned int i = 0; i < bs; i++) {					// calculate partial activation on blocks
			const unsigned int l = (t*bs+i)%ns;
			for (unsigned int j = 0; j < U.N; j++) {
				for (unsigned int k = 0; k < U.N; k++)
					U(i,j) += Xp(l,k)*B(j,k);
				Ut(i,j) = signs(j)*std::tanh(U(i,j))+U(i,j); } }
		for (unsigned int k = 0; k < bs; k++) { 				// learning rule
			for (unsigned int i = 0; i < Q.N; i++)
				for (unsigned int j = 0; j < Q.N; j++)
					Q(i,j) -= Ut(k,i)*U(k,j); }
		Q.addToD(bs); B += Q*B*lr; }						// update weights from this block
}

template <typename T>
matD<T> ica<T>::extended(vecD<T>& signs, const T rate, const bool anneal, const bool verbose)
//! implements the extended infomax algorithm.
{	matD<T> B(pc,pc), Bo(pc,pc), Dt(pc,pc); B.id(); Bo = T(0); Dt = T(0);		// initialize B
	matD<T> Xp = permute(Xw); signs.resize(pc); signs = T(1.0);			// initialize signs vector
	T co = T(1.0), lr = rate;
	for (unsigned int it = 0; it < maxit; it++) {					// iterate...
		learn(B, Xp, signs, lr); setSigns(signs, B);				// learn from permuted blocks
		T dt = T(0), ch = T(0);							// compute convergence statistics
		for (unsigned int i = 0; i < pc; i++)
			for (unsigned int j = 0; j < pc; j++) {
				const T b = B(i,j)-Bo(i,j), d = Dt(i,j);
				dt += b*d; ch += SQR(b); Dt(i,j) = b; };
		const T dc = co*ch, ad = T(std::acos(dt/std::sqrt(dc))*180.0/M_PI);	// change in B and angle B,Bo
		if (verbose) { printf("[%4d] change %.4e angle %.4e\r", it, ch, ad); fflush(stdout); };
		if (anneal && ad < T(2.0)) lr *= T(0.7);				// if annealing, decrease learning rate
		Bo = B; co = ch; if (ch < eps) break; };				// convergence check
	return Xw*B;
}

#endif

