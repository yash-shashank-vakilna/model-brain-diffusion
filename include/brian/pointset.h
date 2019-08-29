#ifndef POINTSET_H
#define POINTSET_H

/*
 *
 * pointset.h: Classification and registration of n-dimensional point sets
 * BRIAN Software Package Version 3.0
 *
 * $Id: pointset.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (12/05/16): initial version FK
 * v406 (28/09/16): bumped to version 3.0
 *
 */

#define allDim(d)	(unsigned int d = 0; d < D; d++) 
#define allX(n)		(unsigned int n = 0; n < N; n++) 
#define allY(m)		(unsigned int m = 0; m < M; m++) 
#define allClasses(c)	(unsigned int c = 0; c < nc; c++)

//! Provides a class for k-means clustering of point sets.

template<typename T> class kmeans {
	const std::vector<vecD<T>>& data;	//!< reference to data points
	const unsigned int N;			//!< number of data points
	const unsigned int nc;			//!< number of clusters

	std::vector<vecD<T>> computeCentroids(const uvecD& lbl)	const			//! recomputes cluster centers.
		{ std::vector<vecD<T>> cen(nc); vecD<T> cc(data[0].N); cc = T(0); 
		  for allClasses(c) { cen[c] = cc; }; uvecD cnt(nc); cnt = 0;		// set centroids to zero
		  for allX(n) {	cen[lbl(n)] += data[n]; cnt[lbl(n)]++; };
		  for allClasses(c) { if (cnt(c)) cen[c] /= cnt(c); }; return cen; }	// divide by count
	unsigned int closestCentroid(const vecD<T>& v, const std::vector<vecD<T>>& cen) const
		//! returns index of centroid closest to vector v.
		{ T dmin = std::numeric_limits<T>::max(); unsigned int cmin = 0;
		  for allClasses(c) { const T d = norm2(cen[c]-v);
			if (d < dmin) { dmin = d; cmin = c; } };
		  return cmin; }
public:
	kmeans(const std::vector<vecD<T>>& _data, const unsigned int _nc)		//! allocates a context for kmeans segmentation of data into nc clusters.
		: data(_data), N(data.size()), nc(_nc) { }
	uvecD work() const								//! performs kmeans classification, returns classification.
		{ uniformDist<float> ud;
		  uvecD lbl(N); for allX(n) lbl[n] = FTOU(ud()*nc);			// initial partition of points
		  for (bool cont = true; cont; ) { cont = false; 
			const std::vector<vecD<T>> cen = computeCentroids(lbl);
			for allX(n) { const unsigned int l = closestCentroid(data[n],cen);
				if (lbl[n] != l) { lbl[n] = l; cont = true; } } }
		  return lbl; }
};

//! Provides a class for parameter estimation of a Gaussian mixture model using variational Bayes 

/*
 * see Wikipedia: http://en.wikipedia.org/wiki/Variational_Bayesian_methods
 *
 * M. Kuusela, Algorithms for variational learning of mixture of Gaussians
 * BS thesis, Helsinki University of Technology, 2009
 *
 */

template<typename T> class gmmVB {
	const T eps = T(1e-8);
	const std::vector<vecD<T>>& data;	//!< reference to data points
	const unsigned int N;			//!< number of data points
	const unsigned int D;			//!< dimensionality of data points
	unsigned int nc;			//!< number of classes
	T	alpha0, beta0, nu0;
	vecD<T>	alpha, beta, nu;
	vecD<T>	mu0;
	matD<T>	Wi0;
	std::vector<vecD<T>> mu, xm;
	std::vector<matD<T>> S, W;
	vecD<T>	pi, lambda, nk;
	matD<T>	R;
	const T removalCriterion;

	T	lnDirichletC(const vecD<T>& alpha, const unsigned int nc)		//! returns log of the Dirichlet normalization constant C.
		{ T as = T(0), cn = T(0);
		  for allClasses(c) { as += alpha(c); cn -= std::lgamma(alpha(c)); }
		  return cn+std::lgamma(as); }
	T	lnWishartB(const matD<T>& W, const T n)					//! returns log of the Wishart normalization constant B.
		{ unsigned int p = W.M; T b = T(-0.25)*T(std::log(M_PI))*p*(p-T(1));
		  for (unsigned int i = 0; i < p; i++) b -= std::lgamma(T(0.5)*(n-i));	// \Gamma_p(n/2)
		  b -= T(0.5)*n*std::log(det(W));
		  b += T(0.5)*n*p*std::log(T(2)); return b; }				// |W|^(n/2) * 2^((n*p)/2)
	void	computeNk()
		{ for allClasses(c) { T s = T(0); for allX(n) { s += R(n,c); };
			nk[c] = s; } }							// compute Nk from r
	void	init();
	void	computeLambda();
	bool	removeClasses();
	void	eStep();
	void	mStep();
	T	costFunction();
	void	print(const vecD<T>& m, const matD<T>& w, const T pp);
public:
	gmmVB(const std::vector<vecD<T>>& _data, const unsigned int _nc, const T rc)
		 : data(_data), N(ITOU(data.size())), D(data[0].N), nc(_nc), alpha0(T(1)),
		beta0(T(1)), nu0(D), alpha(nc), beta(nc), nu(nc), mu0(D), Wi0(), mu(nc),
		xm(nc), S(nc), W(nc), pi(nc), lambda(nc), nk(nc), R(N,nc), removalCriterion(rc)
		{ init();
		  computeLambda(); }
	void	work(const unsigned int nit = 100, const bool verbose = false);
};

template <typename T>
void gmmVB<T>::init()
{	R = T(1)/nc; mu0 = T(0); nu = D;						// init responsibilities
	uniformDist<unsigned int> ud(0,N); for allClasses(c) mu[c] = data[ud()];  	// select random points as mean
	Wi0.resize(D,D); Wi0 = T(0); Wi0.addToD(T(4)/D); for allClasses(c) W[c] = Wi0;	// init precision matrices
	Wi0 = T(0); Wi0.addToD(D/T(4)); alpha = T(1); beta = T(10); 			// init helper variables
}

template <typename T>
void gmmVB<T>::computeLambda()
{	T as = T(0); for allClasses(c) as += alpha(c); 
	for allClasses(c) pi[c] = std::exp(digamma(alpha[c])-digamma(as));		// compute pi
	for allClasses(c) { T t = T(0);
		for allDim(d) t += digamma(T(0.5)*(nu(c)-d));
		lambda[c] = std::exp(t+D*std::log(T(2.0))+std::log(det(W[c]))); };	// compute lambda
}

template <typename T>
bool gmmVB<T>::removeClasses()
//! checks if any classes should be removed.
{	unsigned int l = 0, r = 0;
	for allClasses(c) { if (nk(c) < removalCriterion*N) { r++; continue; };		// remove this column
		if (l != c) { for allX(n) R(n,l) = R(n,c); }; l++; };			// else move column
	assert(nc == l+r); nc = l; return r != 0;
}

template <typename T>
void gmmVB<T>::mStep()
//! estimates auxilliary variables.
{	computeNk(); if (removeClasses()) computeNk();					// compute Nk from r
	for allClasses(c) { vecD<T> s(D); s = T(0);					// compute xm from Nk & r
		for allX(n) s += R(n,c)*data[n];
		xm[c] = s/nk(c); };
	for allClasses(c) { matD<T> s(D,D); s = T(0);					// compute S from xm & Nk
		for allX(n) { vecD<T> d = data[n]-xm[c]; s += R(n,c)*outer(d,d); };
		S[c] = s/nk[c]; };
	for allClasses(c) alpha[c] = alpha0+nk(c);					// compute alpha from Nk
	for allClasses(c) beta[c] = beta0+nk(c);					// compute beta from Nk
	for allClasses(c) mu[c] = (beta0*mu0+nk(c)*xm[c])/beta(c);			// compute mu from xm & beta
	for allClasses(c) nu[c] = nu0+nk(c);						// compute nu from Nk
	for allClasses(c) { vecD<T> d = xm[c]-mu0; matD<T> v = outer(d,d);		// compute W from S, xm & Nk
		v = Wi0+nk[c]*S[c]+v*((beta0*nk(c))/(beta0+nk(c))); W[c] = inv(v); };
}

template <typename T>
void gmmVB<T>::eStep()
//! maximizes responsibilities.
{	for allX(n) { T s = T(0); bool nm = false;					// compute responsibilites
		for allClasses(c) { vecD<T> d = data[n]-mu[c];
			const T t = D/beta(c)+nu(c)*dot(d,W[c]*d);
			T r = pi(c)*std::sqrt(lambda(c))*std::exp(T(-0.5)*t);
			if (r == T(0)) r = T(1);
			s += r; R(n,c) = r; };
		for allClasses(c) { T r = R(n,c)/s;					// normalization
			if (r < eps) { r = eps; nm = true; }; R(n,c) = r; };		// check for small values
		if (nm) { T t = T(0); for allClasses(c) t += R(n,c);			// re-normalize if a small value was found
			for allClasses(c) R(n,c) /= t; } };
}

template <typename T>
T gmmVB<T>::costFunction()
//! estimates log q.
{	T cn = T(0); for allX(n) { for allClasses(c) cn += R(n,c)*std::log(R(n,c)); };	// Eq. 61
	for allClasses(c) cn += (alpha(c)-alpha0)*std::log(pi(c));			// Eq. 62
	cn += lnDirichletC(alpha,nc); vecD<T> al(nc); al = alpha0;
	cn -= lnDirichletC(al,nc);							// Eq. 66
	for allClasses(c) cn += lnWishartB(W[c],nu(c));
	cn -= nc*lnWishartB(inv(Wi0),nu0);
	T t = T(0); for allClasses(c) { const vecD<T> d = mu[c]-mu0; 			// Eq. 67
		t += D*(beta0/beta(c)-std::log(beta0/beta(c))-nu(c)-T(1));
		t += beta0*nu(c)*dot(d,W[c]*d);
		t += nu(c)*(Wi0*W[c]).trace()+(nu(c)-nu0)*std::log(lambda(c)); };
	cn += T(0.5)*t; t = T(0); for allClasses(c) { vecD<T> d = xm[c]-mu[c]; 		// Eq. 64
		T t1 = std::log(lambda(c))+T(2)*std::log(pi(c))-D/beta(c)-D*std::log(T(2)*M_PI);
		T t2 = (S[c]*W[c]).trace()+dot(d,W[c]*d);
    		t += nk(c)*(t1-nu(c)*t2); };
	cn -= T(0.5)*t; return cn;	
}

template <typename T>
void gmmVB<T>::print(const vecD<T>& m, const matD<T>& w, const T pp)
{	printf("mean: %e ", pp); for allDim(d) printf("%7.4e ", m(d)); 
	printf("\nsigma:"); const matD<T> s = inv(w*N);
	for allDim(i) { for allDim(j) printf(" %7.4e", s(i,j)); printf(";"); };
	printf("\n");
}

template <typename T>
void gmmVB<T>::work(const unsigned int nit, const bool verbose)
{	T prevCost = std::numeric_limits<T>::max();
	for (unsigned int it = 0; it < nit; it++) {
		eStep();
		mStep();
		computeLambda();
		const T cost = costFunction();
		if (std::abs(prevCost-cost) < eps) break;
		prevCost = cost;
		if (verbose) { printf("[%u] cost %6.2f %u\r", it, cost, nc); fflush(stdout); } };
	if (verbose) { for allClasses(c) print(mu[c],W[c],nk[c]/N); };
}

//! Provides a class for mean shift classification 

template<typename T> class meanShift {
	const T eps = T(1e-8);
	const std::vector<vecD<T>>& data;	//!< reference to data points
	const unsigned int N;			//!< number of data points
	const unsigned int D;			//!< dimensionality of data points
        const T h;				//!< window radius
        vecD<T>	xm;				//!< window center

        T   kernel(const vecD<T>& di, const T beta, const T lambda) const		//! truncated Gaussian kernel.
                { const T nm = norm(di);
		  return nm < lambda? T(0): std::exp(-beta*SQR(nm)); };
        vecD<T>	shift() const
                { vecD<T> nom(D); nom = T(0); T den = T(0);
                  for allX(n) { const T gi = T(-2)*kernel(n);
			den += gi; nom += data[n]*gi; };
                  return isSmall<T>(den)? nom: (nom/den)-xm; }
        T   density() const
                { T s = T(0); for allX(n) s += kernel(n);
		  return s/(N*std::pow(h,D)); }

public:
        meanShift(const std::vector<vecD<T>>& _data, const T _h)
        : data(_data), N(ITOU(data.size())), D(data[0].N), h(_h), xm(D) { }
        vecD<T>	work(const unsigned int nit = 100, const bool verbose = false)
		{ xm = data[0];
        	  for (unsigned int it = 0; it < nit; it++) {
                	const vecD<T> dx = shift(); const T nm = norm(dx);
                	if (verbose) { printf("[%2d] %e %e\n", it, nm, density()); fflush(stdout); };
                	if (nm < eps) break;
			xm += dx; };
        	  return xm; }
};


//! Provides a class for affine registration of point sets.

template<typename T> class cpdAffine {
	const std::vector<vecD<T>>& X;		//!< set of data vectors X
	const std::vector<vecD<T>>& Y;		//!< set of data vectors Y
	const unsigned int M;			//!< number of observations in Y
	const unsigned int N;			//!< number of observations in X
	const unsigned int D;			//!< dimensionality of { X, Y }
	const matD<T> Xm;			//!< data X as matrix
	const matD<T> Ym;			//!< data Y as matrix
	matD<T>	B;				//!< affine transformation matrix
	vecD<T>	t;				//!< shift vector

	T	prob(const unsigned int n, const unsigned int m, const T s2) const	//! returns probability of match { x,T(y) }.
		{ return std::exp(T(-0.5)*norm2(X[n]-(B*Y[m]+t))/s2); }
	matD<T>	matX() const								//! converts data X into matrix.
		{ matD<T> A(N,D); for allX(n) for allDim(d) A(n,d) = X[n](d); return A; }
	matD<T>	matY() const								//! converts data Y into matrix.
		{ matD<T> A(N,D); for allY(m) for allDim(d) A(m,d) = Y[m](d); return A; }
	T	initSigma() const							//! returns initial standard deviation.
		{ T s2 = T(0); for allX(n) for allY(m) s2 += norm2(X[n]-Y[m]);
		  return s2/T(N*M*D); }
	matD<T>	eStep(const T w, const T s2) const					//! E-step: computes probabilities P.
		{ const T ws = w*M*std::pow(T(2*M_PI)*s2,T(0.5)*D)/((T(1)-w)*N);
		  matD<T> P(M,N); for allX(n) { T pn = ws;
			for allY(m) { const T p = prob(n,m,s2); pn += p; P(m,n) = p; }
		  	for allY(m) P(m,n) /= pn; };
		  return P; }
	T	mStep(const matD<T>& P)							//! M-step: computes B,t,s2.
		{ vecD<T> P1 = P.rowsum(), Pt1 = P.colsum(); const T np = P.sum();
		  vecD<T> mux = trp(Xm)*Pt1/np, muy = trp(Ym)*P1/np;			// get mu_x, mu_y
		  matD<T> Xh(N,D); for allX(n) for allDim(d) Xh(n,d) = Xm(n,d)-mux(d);	// subtract mu_x from X
		  matD<T> Yh(M,D); for allY(m) for allDim(d) Yh(m,d) = Ym(m,d)-muy(d);	// subtract mu_y from Y
		  matD<T> XP(N,D); for allX(n) for allDim(d) XP(n,d) = Xh(n,d)*Pt1(n);	// X^T P^T 1
		  matD<T> YP(M,D); for allY(m) for allDim(d) YP(m,d) = Yh(m,d)*P1(m);	// Y^T P 1
		  const matD<T> A = trp(Xh)*trp(P)*Yh;
		  B = A*inv(trp(Yh)*YP); t = mux-B*muy;					// update B and t
		  const T t1 = (trp(XP)*Xh).trace(), t2 = (A*trp(B)).trace();
		  return std::max((t1-t2)/(np*D),T(0)); }				// return new s2
	std::vector<vecD<T>> transform() const						// returns transformed data
		{ std::vector<vecD<T>> Yt(M); for allY(m) Yt[m] = B*Y[m]+t; return Yt; }
public:
	cpdAffine(const std::vector<vecD<T>>& _X, const std::vector<vecD<T>>& _Y)	//! allocates a context for affine point set matching.
		: X(_X), Y(_Y), M(Y.size()), N(X.size()), D(X[0].N),
		Xm(matX()), Ym(matY()), B(D,D), t(D)
		{ assert(M > 0 && N > 0); assert(X[0].N == Y[0].N); B.id(); t = T(0); }
	std::vector<vecD<T>> work(const T w, const unsigned int nit = 100, 
			const bool verbose = false)					//! returns transformed Y given uniformity weight w.
		{ const T wc = CLAMP(w,T(1e-6),T(1)); T s2 = initSigma();		// avoid zero weight
		  for (unsigned int it = 0; it < nit; it++) {
			matD<T> P = eStep(wc,s2); s2 = mStep(P); assert(std::isfinite(s2));
			if (verbose) { printf("[%2d] %.4e\r", it, s2); fflush(stdout); };
			if (s2 < T(1e-8)) break; };
		  return transform(); }
	matD<T>	getMatrix() const
		{ return B; }
	vecD<T>	getTranslation() const
		{ return t; }
};

//! Provides a class for non-rigid registration of point sets.

template<typename T> class cpdNonRigid {
	const std::vector<vecD<T>>& X;		//!< set of data vectors X
	const std::vector<vecD<T>>& Y;		//!< set of data vectors Y
	const unsigned int M;			//!< number of observations in Y
	const unsigned int N;			//!< number of observations in X
	const unsigned int D;			//!< dimensionality of { X,Y }
	const matD<T> Xm;			//!< data X as matrix
	const matD<T> Ym;			//!< data Y as matrix
	matD<T>	G;				//!< kernel matrix
	matD<T>	W;				//!< weights matrix
	matD<T>	TY;				//!< transformed Y

	T	prob(const unsigned int n, const unsigned int m, const T s2) const	//! returns probability of match { x,T(y) }.
		{ return std::exp(T(-0.5)*norm2(X[n]-TY.getRow(m))/s2); }
	matD<T>	matX() const								//! converts data X into matrix.
		{ matD<T> A(N,D); for allX(n) for allDim(d) A(n,d) = X[n](d); return A; }
	matD<T>	matY() const								//! converts data Y into matrix.
		{ matD<T> A(N,D); for allY(m) for allDim(d) A(m,d) = Y[m](d); return A; }
	T	initSigma() const							//! returns initial standard deviation.
		{ T s2 = T(0); for allX(n) for allY(m) s2 += norm2(X[n]-Y[m]);
		  return s2/T(N*M*D); }
	void	initKernel(const T beta)						//! initializes kernel matrix G.
		{ for allY(i) for allY(j) {
			if (i == j) G(i,i) = T(1);
			else { const T g = std::exp(T(-0.5)*SQR(norm(Y[i]-Y[j])/beta));
				G(i,j) = g; G(j,i) = g; } } }
	matD<T>	eStep(const T w, const T s2) const					//! E-step: computes probabilities P.
		{ const T ws = w*M*pow(T(2.0*M_PI)*s2,T(0.5)*D)/((T(1)-w)*N);
		  matD<T> P(M,N); for allX(n) { T pn = ws;
			for allY(m) { const T p = prob(n,m,s2); pn += p; P(m,n) = p; }
		  	for allY(m) P(m,n) /= pn; };
		  return P; }
	T	mStep(const matD<T>& P, const T lambda, const T s2)			//! M-step: computes weight matrix
		{ const vecD<T> P1 = P.rowsum(), Pt1 = P.colsum();
		  matD<T> A = G, B = P*Xm;						// setup equation system
		  for allY(m) { A(m,m) += lambda*s2/P1(m); const T p = P1(m);
			for allDim(d) B(m,d) = B(m,d)/p-Ym(m,d); };
		  W = inv(A)*B;	TY = Ym+G*W;						// solve for W and transform Y using new weights
		  matD<T> XtPt(D,N); for allX(n) for allDim(d) XtPt(d,n) = X[n](d)*Pt1(n); // X^T P^T 1
		  matD<T> Tt(D,M); for allY(m) for allDim(d) Tt(d,m) = TY(m,d)*P1(m);	// T^T P 1
		  const T t1 = (XtPt*Xm).trace();					// trace X^T P^T 1 X
		  const T t2 = T(-2)*(trp(P*Xm)*TY).trace();				// trace (P X)^T T
		  const T t3 = (Tt*TY).trace(), np = P1.sum();				// trace T^T P 1 T
		  return std::max((t1+t2+t3)/(np*D),T(0)); }				// return estimate for s2.
	std::vector<vecD<T>> transform() const						//! returns transformed data
		{ std::vector<vecD<T>> Yt(M); for allY(m) Yt[m] = TY.getRow(m);
		  return Yt; }
public:
	cpdNonRigid(const std::vector<vecD<T>>& _X, const std::vector<vecD<T>>& _Y)	//! allocates a context for nonlinear point set matching.
		: X(_X), Y(_Y), M(Y.size()), N(X.size()), D(X[0].N),
		Xm(matX()), Ym(matY()), G(M,M), W(M,D), TY()
		{ assert(M > 0 && N > 0); assert(X[0].N == Y[0].N); W = T(0); }
	std::vector<vecD<T>> work(const T w, const T beta = T(2), const T lambda = T(2),
			const unsigned int nit = 100, const bool verbose = false)	//! returns transformed Y given uniformity weight w, smoothness beta and lambda.
		{ const T wc = CLAMP(w,T(1e-6),T(1)); T s2 = initSigma();		// avoid zero weight
		  initKernel(beta); TY = Ym+G*W;
		  for (unsigned int it = 0; it < nit; it++) {
			matD<T> P = eStep(wc,s2);
			s2 = mStep(P,lambda,s2); assert(std::isfinite(s2));
			if (verbose) { printf("[%2d] %.4e\r", it, s2); fflush(stdout); };
			if (s2 < T(1e-8)) break; };
		  return transform(); }
	matD<T>	getTransform() const
		{ return G*W; }
};

//! Provides a class for finding k nearest neighbors in a point set

template<typename T> class approxNN {
	using uvec = std::vector<unsigned int>;
	using rvec = std::vector<std::pair<T,unsigned int>>;

	//! Represents a node in a approxNN tree.

	struct Node {
		uvec	ch;			//!< children of this node
		T	a;			//!< plane offset
		vecD<T>	v;			//!< plane normal
		unsigned int type;		//!< 0: undefined, 1: terminal, 2: internal

		Node()									//! allocates an undetermined node.
			: ch(), a(0), v(), type(0) { }
		Node(const uvec& ix)							//! allocates a terminal node.
			: ch(ix), a(0), v(), type(1) { }
		Node(const vecD<T>& w, T _a = 0)					//! allocates an internal node.
			: ch(2), a(_a), v(w), type(2) { ch[0] = 0; ch[1] = 0; }
		T	margin(const vecD<T>& y) const					//! returns the split limit of this node.
			{ assert(type == 2); return a+dot(v,y); }
		bool	side(const vecD<T>& v) const					//! returns the side of v with respect to the split limit.
			{ const T d = margin(v); return d != T(0)? d > T(0): rand()&1; }
	};

	const std::vector<vecD<T>>& data;	//!< data items
	const unsigned int B;			//!< bucket size of terminal node
	std::vector<Node> nodes;		//!< nodes
	unsigned int root;			//!< index to root node

	bool	normalize(vecD<T>& v) const						//! scales vector v to unit length.
		{ const T n = norm(v); if (n == T(0)) return false;
		  v /= n; return true; }
	Node	internalNode(const uvec& ix)						//! builds an internal node.
		{ const unsigned int N = ix.size(), i = rand()%N, j = rand()%(N-1);
		  vecD<T> vi = data[ix[i]]; unsigned int ni = 1; 
		  vecD<T> vj = data[ix[j+(j >= i)]]; unsigned int nj = 1;
		  for (unsigned int it = 0; it < 200; it++) {
			const vecD<T>& v = data[ix[rand()%N]];
			const T di = ni*norm(vi-v), dj = nj*norm(vj-v);
			if (di < dj) { vi = (vi*T(ni)+v)/T(ni+1); ni++; }
			else { vj = (vj*T(nj)+v)/T(nj+1); nj++; } };
		  vecD<T> v = vi-vj; if (normalize(v) == false) return Node(ix);	// if data in ix is homogeneous, return terminal node
		  const T a = dot(T(-0.5)*(vi+vj),v);
		  assert(std::isfinite(a)); Node n = Node(v,a); uvec ch[2];
		  for (const auto i : ix) {						// sort nodes by side
			const unsigned int s = n.side(data[i]); ch[s].push_back(i); }
		  while (ch[0].size() == 0 || ch[1].size() == 0) {			// if either side has no children
			ch[0].clear(); ch[1].clear(); n.v = T(0);			// clear previous assignment
			for (const auto i : ix)	ch[rand()&1].push_back(i); }		// and assign randomly
		  const unsigned int f = (ch[0].size() > ch[1].size());
		  for (unsigned int s = 0; s < 2; s++) n.ch[s^f] = makeTree(ch[s^f]);
		  return n; }
	unsigned int makeTree(const uvec& ix)						//! builds node tree recursively.
		{ const Node nd = ix.size() < B? Node(ix): internalNode(ix);		// make terminal or internal node
		  nodes.push_back(nd); return nodes.size()-1; }				// and add to list
	void	search(rvec& res, const vecD<T>& v, const unsigned int n) const
		{ std::priority_queue<std::pair<T,unsigned int>> q;			// initialize queue
		  q.push({std::numeric_limits<T>::infinity(),root});
		  while (res.size() < n*B && q.empty() == false) {			// while insufficient items found
			const std::pair<T,unsigned int> t = q.top(); q.pop();
			const T d = t.first; const Node& nd = nodes[t.second];
			switch (nd.type) {						// check node
			case 1: for (const auto i : nd.ch)				// terminal: add all items to result set
					res.push_back({norm(v-data[i]),i});
				break;
			case 2: { const T m = nd.margin(v);				// internal: add both sides to queue
				  q.push({std::min(d,+m),nd.ch[1]});
				  q.push({std::min(d,-m),nd.ch[0]}); };
				break;
			default: throw rtException("undefined node %u", t.second); } };	// should not happen
		  std::sort(res.begin(), res.end());					// sort by increasing distance
		  std::unique(res.begin(), res.end()); }				// remove duplicates
	uvecD	convertResults(vecD<T>* dist, const unsigned int n, const rvec& res) const
		//! converts list of result pairs into vectors
		{ const unsigned int m = res.size(), p = std::min(n,m);			// determine size of result set
		  uvecD ids(p); if (dist) dist->resize(p);				// save ids and distances
		  for (unsigned int i = 0; i < p; i++) {
			(*dist)[i] = res[i].first; ids[i] = res[i].second; };
		  return ids; }
public:
	approxNN()									//! allocates an empty tree for nearest neighbor search.
		 : data(), B(), nodes(), root(0) { }
	approxNN(const std::vector<vecD<T>>& _data, const unsigned int _B = 16)		//! allocates a context for approximate nearest neighbor search.
		 : data(_data), B(_B), nodes(), root(0)
		{ const unsigned int n = data.size(); nodes.reserve(n*2/B);		// allocate approximate number of nodes
		  uvec ix(n); for (unsigned int i = 0; i < ix.size(); i++) ix[i] = i;
		  root = makeTree(ix); }
	uvecD	search(const vecD<T>& v, const unsigned int n, vecD<T>* dist = nullptr) const
		//! returns at most n nearest neighbors to v and their distances.
		{ rvec res; search(res,v,n);
		  for (auto it = res.begin(); it != res.end(); )			// remove items with 0 distance from set
			if (it->first == T(0)) res.erase(it); else it++;
		  return convertResults(dist,n,res); }
	uvecD	search(const unsigned int id, const unsigned int n, vecD<T>* dist = nullptr) const
		//! returns at most n nearest neighbors to data item id and their distances.
		{ rvec res; search(res,data[id],n);
		  for (auto it = res.begin(); it != res.end(); it++)			// remove item id from set
			if (it->second == id) { res.erase(it); break; };
		  return convertResults(dist,n,res); }
	T	dist(const vecD<T>& v, const unsigned int n) const			//! returns mean distance to k nearest neighbors.
		{ vecD<T> dist; search(v,n,&dist); const T d = dist.mean(); 
		  assert(std::isfinite(d) && d >= T(0)); return d; }
	T	dist(const unsigned int id, const unsigned int n) const			//! returns mean distance to k nearest neighbors.
		{ vecD<T> dist; search(id,n,&dist); const T d = dist.mean(); 
		  assert(std::isfinite(d) && d >= T(0)); return d; }
};

/*
 * Graph-based similarity measures
 *
 * refer to:
 * Neemuchwala H (2005), Entropic graphs for image registration.
 * PhD thesis UoM.
 *
 */

template<typename T>
T alphaGA(const std::vector<vecD<T>>& a, const std::vector<vecD<T>>& b,
	const T alpha = T(0.99), const unsigned int nn = 1)
//! returns the alpha geometric arithmetic mean average of feature sets (a,b).
// alphaGA == 0 if distribution(a) == distribution(b)
{	const unsigned int N = a.size(); assert(b.size() == N);				// both sample must have the same size
	if (alpha <= T(0) || alpha >= T(1)) throw optException("alpha is outside of (0,1)");
	const T D = a[0].N, gamma = D*(T(1)-alpha), eps = T(1e-6); assert(D > T(0));
	approxNN<T> nna(a), nnb(b); T p = T(0); 						// allocate structures for nearest neighbor search
	for (unsigned int i = 0; i < N; i++) {						// sum p over all features
		const T ea = std::max(nna.dist(i,nn),eps);				// edge length of nn nearest neighbors in a
		const T eb = std::max(nnb.dist(i,nn),eps);				// edge length of nn nearest neighbors in b
		const T q = ea < eb? ea/eb: eb/ea; p += std::pow(q,T(0.5)*gamma); };	// see Eq 4.11
	return p > T(0)? std::log(p/N)/(T(1)-alpha):
		std::log(std::numeric_limits<T>::min());
}	

template<typename T>
T alphaMI(const std::vector<vecD<T>>& a, const std::vector<vecD<T>>& b,
	const T alpha = T(0.99), const unsigned int nn = 1)
//! returns the Renyi alpha entropy of feature sets (a,b).
{	const unsigned int N = a.size(); assert(b.size() == N);				// both sample must have the same size
	if (alpha <= T(0) || alpha >= T(1)) throw optException("alpha is outside of (0,1)");
	const T D = a[0].N, gamma = D*(T(1)-alpha); assert(D > T(0));
	std::vector<vecD<T>> c(N);
	for (unsigned int i = 0; i < N; i++) c[i] = stack(a[i],b[i]);			// concatenate a and b
	approxNN<T> nna(a), nnb(b), nnc(c); T p = T(0); 					// allocate structures for nearest neighbor search
	for (unsigned int i = 0; i < N; i++) {						// sum p over all features
		const T ea = nna.dist(i,nn);						// edge length of nn nearest neighbors in a
		const T eb = nnb.dist(i,nn);						// edge length of nn nearest neighbors in b
		const T ec = nnc.dist(i,nn);						// edge length of nn nearest neighbors in a x b
		const T q = ec/std::max(std::sqrt(ea*eb),T(1e-6));
		p += std::pow(q,T(2)*gamma); };						// see Eq 4.12
	return p > T(0)? std::log(p/std::pow(N,alpha))/(T(1)-alpha):
		std::log(std::numeric_limits<T>::min());
}	

template<typename T>
T nlcc(const std::vector<vecD<T>>& a, const std::vector<vecD<T>>& b, const unsigned int nn = 1)
//! returns nonlinear correlation coefficient of feature sets (a,b).
{	const unsigned int N = a.size(); assert(b.size() == N);				// both sample must have the same size
	std::vector<vecD<T>> c(N);
	for (unsigned int i = 0; i < N; i++) c[i] = stack(a[i],b[i]);			// concatenate a and b
	approxNN<T> nna(a), nnb(b), nnc(c); T p = T(0); 					// allocate structures for nearest neighbor search
	for (unsigned int i = 0; i < N; i++) {						// sum p over all features
		const T ea = nna.dist(i,nn);						// edge length of nn nearest neighbors in a
		const T eb = nnb.dist(i,nn);						// edge length of nn nearest neighbors in b
		const T ec = nnc.dist(i,nn);						// edge length of nn nearest neighbors in a x b
		p += ec/std::max(std::sqrt(ea*eb),T(1e-6)); };				// see Eq 4.24
	const T rho = p/N; return rho/(rho+T(1));					// NB: optimum is a minimum
}	


template<typename T>
T cc(const std::vector<vecD<T>>& a, const std::vector<vecD<T>>& b)
//! returns linear correlation coefficient of feature sets (a,b).
{	const unsigned int N = a.size(), D = a[0].N; assert(b.size() == N);		// both sample must have the same size
	vecD<T> ma(D); for (unsigned int i = 0; i < N; i++) ma += a[i]; ma /= N;	// mean of set a
	vecD<T> mb(D); for (unsigned int i = 0; i < N; i++) mb += b[i]; mb /= N;	// mean of set b
	T eab = T(0), ea = T(0), eb = T(0);
	for (unsigned int i = 0; i < N; i++) {
		const vecD<T> da = a[i]-ma, db = b[i]-mb;
		eab += dot(da,db); ea += dot(da,da); eb += dot(db,db); };
	return T(1)-eab/std::max(std::sqrt(ea*eb),T(1e-6));				// NB: use difference s.t. optimum is a minimum
}	

#endif

