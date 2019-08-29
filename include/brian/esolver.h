#ifndef ESOLVER_H
#define ESOLVER_H

/*
 *
 * esolver.h: classes for solving sparse eigenvalue problems
 * BRIAN Software Package Version 3.0
 *
 * $Id: esolver.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (25/03/17) for BRIAN3
 *
 */

/*! \file
    \brief Implements classes for solving sparse eigenvalue problems.
*/

//! Helper class for sorting eigenvalues

enum SELECT_EIGENVALUE
{	LARGEST_MAGN,				///< largest magnitude
	LARGEST_REAL,				///< largest real part (for general eigensolvers)
	LARGEST_IMAG,				///< largest imaginary part (in magnitude, for general eigen solvers)
	LARGEST_ALGE,      			///< largest algebraic value
	SMALLEST_MAGN,     			///< smallest magnitude
	SMALLEST_REAL,				///< smallest real part (for general eigensolvers)
	SMALLEST_IMAG,				///< smallest imaginary part (in magnitude, for general eigen solvers)
	SMALLEST_ALGE,				///< smallest algebraic value
};

template <typename X> struct ElemType {
	typedef X type;
};
template <typename X> struct ElemType< std::complex<X>> {
	typedef X type;
};
template <typename X, int SelectionRule> struct SortingTarget {
	static typename ElemType<X>::type get(const X& val)
		{ throw std::invalid_argument("incompatible selection rule");
        	  return -std::abs(val); }
};
template <typename X> struct SortingTarget<X,LARGEST_MAGN> {
	static typename ElemType<X>::type get(const X& val)
		{ return -std::abs(val); }
};
template <typename X> struct SortingTarget<std::complex<X>,LARGEST_REAL> {
	static X get(const std::complex<X>& val)
		{ return -val.real(); }
};
template <typename X> struct SortingTarget<std::complex<X>,LARGEST_IMAG> {
	static X get(const std::complex<X>& val)
		{ return -std::abs(val.imag()); }
};
template <typename X> struct SortingTarget<X,LARGEST_ALGE> {
	static X get(const X& val)
		{ return -val; }
};
template <typename X> struct SortingTarget<X,SMALLEST_MAGN> {
	static typename ElemType<X>::type get(const X& val)
		{ return std::abs(val); }
};
template <typename X> struct SortingTarget<std::complex<X>,SMALLEST_REAL> {
	static X get(const std::complex<X>& val)
		{ return val.real(); }
};
template <typename X> struct SortingTarget<std::complex<X>,SMALLEST_IMAG> {
	static X get(const std::complex<X>& val)
		{ return std::abs(val.imag()); }
};
template <typename X> struct SortingTarget<X,SMALLEST_ALGE> {
	static X get(const X& val)
		{ return val; }
};
template <typename PairType> struct PairComparator {
	bool operator()(const PairType& v1, const PairType& v2)
		{ return v1.first < v2.first; }
};
template <typename X, unsigned int SelectionRule> class SortEigenvalue {
private:
	using TargetType = typename ElemType<X>::type;					// type of the sorting target, e.g. "double"
	using PairType = std::pair<TargetType,unsigned int>;				// type of the sorting pair
	std::vector<PairType> v;
public:
	SortEigenvalue(const X* start, const unsigned int size)
		: v(size)
		{ for (unsigned int i = 0; i < size; i++)
			v[i] = std::make_pair(SortingTarget<X,SelectionRule>::get(start[i]), i);
          	  PairComparator<PairType> comp;
          	  std::sort(v.begin(),v.end(),comp); }
	std::vector<unsigned int> index()
		{ std::vector<unsigned int> ix(v.size());
        	  for (unsigned int i = 0; i < v.size(); i++) ix[i] = v[i].second;
		  return ix; }
};

//! Implements a context for performing a QR decomposition of a tridiagonal matrix.

template<typename T> class deflateTri {
	unsigned int n;				//!< matrix dimension
	matD<T>	M;				//!< copy of the matrix to be factorized
	T	t;				//!< shift constant
	vecD<T> rcos;				//!< rotation cosine
	vecD<T> rsin;				//!< rotation sine

	T	setRot(const T a, const T b, unsigned int i)				//! computes a Givens rotation
		{ T r = hypot(a,b);
        	  if (r <= std::numeric_limits<T>::epsilon()) {
			r = T(0); rcos[i] = T(1); rsin[i] = T(0); }
        	  else { rcos[i] = a/r; rsin[i] = -b/r; }
		  return r; }
	void	rot(T *p, const unsigned int i)						//! applies a Givens rotation
            	{ T t = *p;
            	  p[0] = rcos(i)*t-rsin(i)*p[1];
        	  p[1] = rsin(i)*t+rcos(i)*p[1]; }
public:
	deflateTri(const matD<T>& _M, const T _t)
		: n(_M.M), M(n,n), t(_t), rcos(n-1), rsin(n-1)
		{ M = T(0); T *Tii = M.x, *p;
		  for (unsigned int i = 0; i < n; i++) M(i,i) = _M(i,i)-t;
		  for (unsigned int i = 0; i < n-1; i++) M(i,i+1) = _M(i+1,i);
		  for (unsigned int i = 0; i < n-1; i++) M(i+1,i) = _M(i+1,i);
        	  for (unsigned int i = 0; i < n-2; i++) {
        		T r = setRot(Tii[0], Tii[1], i);
            		Tii[0] = r; Tii[1] = T(0); p = Tii+n; rot(Tii+n,i);
        		p += n; p[0] = -rsin[i] * p[1]; p[1] *= rcos[i];
        		Tii += n+1; }
        	  T r = setRot(Tii[0], Tii[1], n-2);
        	  Tii[0] = r; Tii[1] = T(0); rot(Tii+n,n-2); }
	void	YQ(matD<T>& Q) const							//! returns decomposition matrix Q.
		{ for (unsigned int i = 0; i < n-1; i++) {
        	  	T *Qi = Q.x+i*Q.M, *Qi1 = Q.x+(i+1)*Q.M;
        	  	for (unsigned int j = 0; j < Q.M; j++) { const T t = Qi[j];
                		Qi[j]  = rcos(i)*t-rsin(i)*Qi1[j];
                		Qi1[j] = rsin(i)*t+rcos(i)*Qi1[j]; } } }
	matD<T>	QtHQ() const								//! returns deflated matrix.
		{ matD<T> RQ(n,n); RQ = T(0);
		  for (unsigned int i = 0; i < n; i++) RQ(i,i) = M(i,i);
		  for (unsigned int i = 0; i < n-1; i++) RQ(i,i+1) = M(i,i+1);
        	  T *m11 = RQ.x, *m12, *m21, *m22;
		  for (unsigned int i = 0; i < n-1; i++) {
			m21 = m11+1; m12 = m11+n; m22 = m12+1; const T t = *m21;
			*m11 = rcos(i) * (*m11) - rsin(i) * (*m12);
			*m21 = rcos(i) * t - rsin(i) * (*m22);
			*m22 = rsin(i) * t + rcos(i) * (*m22); m11 = m22; }
		  for (unsigned int i = 0; i < n-1; i++) RQ(i,i+1) = RQ(i+1,i);
		  for (unsigned int i = 0; i < n; i++) RQ(i,i) += t;
		  return RQ; }
};

//! Implements a context for performing a single shift QR decomposition.

template <typename T> class deflateSingle {
	unsigned int n;				//!< matrix dimension
	matD<T>	H;				//!< copy of the matrix to be factorized
	T	t;				//!< shift constant
	vecD<T> rcos;				//!< rotation cosine
	vecD<T> rsin;				//!< rotation sine
public:
	deflateSingle(const matD<T>& _H, const T _t)					//! allocates a context for a QR decomposition of matrix H with shift t.
		: n(_H.N), H(_H), t(_t), rcos(n-1), rsin(n-1)
		{ T c, s, eps = std::numeric_limits<T>::epsilon();
		  for (unsigned int i = 0; i < n; i++) H(i,i) -= t;
		  for (unsigned int i = 0; i < n-1; i++) { T *Tii = H.x+i*n+i;
			for (unsigned int j = i+2; j < n; j++) H(j,i) = T(0);
			T xi = Tii[0], xj = Tii[1], r = std::hypot(xi,xj);
			if (r <= eps) { r = T(0); rcos[i] = c = T(1); rsin[i] = s = T(0); }
			else { rcos[i] = c = xi/r; rsin[i] = s = -xj/r; }
			Tii[0] = r; Tii[1] = T(0); T *p = Tii+n;
			for (unsigned int j = i+1; j < n; j++, p += n) { T t = p[0];
				p[0] = c*t-s*p[1]; p[1] = s*t+c*p[1]; } } }
	void	YQ(matD<T>& Q)								//! returns decomposition matrix Q.
		{ for (unsigned int i = 0; i < n-1; i++) {
			T *Qi = Q.x+i*Q.M, *Qi1 = Q.x+(i+1)*Q.M;
			for (unsigned int j = 0; j < Q.M; j++) { T t = Qi[j];
                		Qi[j]  = rcos(i)*t-rsin(i)*Qi1[j];
                		Qi1[j] = rsin(i)*t+rcos(i)*Qi1[j]; } } }
	matD<T>	QtHQ() const								//! returns deflated matrix.
		{ matD<T> RQ = H; 
		  for (unsigned int i = 0; i < n-1; i++) {
			T *Yi = RQ.x+i*n, *Yi1 = Yi+n;
			for (unsigned int j = 0; j < i + 2; j++) { T t = Yi[j];
				Yi[j]  = rcos(i)*t-rsin(i)*Yi1[j];
				Yi1[j] = rsin(i)*t+rcos(i)*Yi1[j]; } }
		  for (unsigned int i = 0; i < n; i++) RQ(i,i) += t;
		  return RQ; }
};

//! Implements a context for performing a double shift QR decomposition.

template <typename T> class deflateDouble {
	const T	eps = std::numeric_limits<T>::epsilon();
	unsigned int n;				//!< matrix dimension
	matD<T> H;				//!< copy of the matrix to be factorized
	T	s;	 			//!< shift constant
	T	t;				//!< shift constant
	matD<T>	R;	 			//!< Householder reflectors
	uvecD	nr;				//!< # rows each reflector affects

	void	computeReflector(const T& x1, const T& x2, const T& x3,
			const unsigned int i)						//! computes reflectors, stores in column i of R.
		{ nr[i] = 3; T x2x3 = T(0);						// in general case the reflector affects 3 rows
		  if (std::abs(x3) < eps) {						// If x3 is zero, decrease nr by 1
			if (std::abs(x2) < eps) { nr[i] = 1; return; }			// if x2 is also zero, nr will be 1.
			nr[i] = 2; x2x3 = std::abs(x2); }
		  else x2x3 = std::hypot(x2,x3);
		  const T x1n = x1-((x1 <= 0)-(x1 > 0))*std::hypot(x1,x2x3);
		  const T xn = std::hypot(x1n,x2x3);
		  if (xn < eps) { nr[i] = 1; return; }					// double check the norm of new x
		  R(0,i) = x1n/xn; R(1,i) = x2/xn; R(2,i) = x3/xn; }
	void	doPX(matD<T>& M, const unsigned int is, const unsigned int ie, const unsigned int js, 
			const unsigned int je, const unsigned int r)			//! applies PX rotation
		{ if (nr[r] == 1) return;
		  if (nr[r] == 2 || ie-is == 2) {
			for (unsigned int j = js; j < je; j++) {
				T t = T(2)*R(0,r)*M(is,j)+T(2)*R(1,r)*M(is+1,j);
				M(is,j) -= t*R(0,r); M(is+1,j) -= t*R(1,r); } }
		  else { for (unsigned int j = js; j < je; j++) {
				T t = T(2)*R(0,r)*M(is,j)+T(2)*R(1,r)*M(is+1,j)+T(2)*R(2,r)*M(is+2,j);
				M(is,j) -= t*R(0,r); M(is+1,j) -= t*R(1,r);
				M(is+2,j) -= t*R(2,r); } } }
	void	doXP(matD<T>& M, const unsigned int is, const unsigned int ie, const unsigned int js, 
			const unsigned int je, const unsigned int r)			//! applies XP rotation
		{ if (nr[r] == 1) return;
		  if (nr[r] == 2 || je-js == 2) {
			for (unsigned int i = is; i < ie; i++) {
				T t = T(2)*R(0,r)* M(i,js) + T(2)*R(1,r)*M(i,js+1);
				M(i,js) -= t*R(0,r); M(i,js+1) -= t*R(1,r); } }
		  else { for (unsigned int i = is; i < ie; i++) {
				T t = T(2)*R(0,r)*M(i,js)+T(2)*R(1,r)*M(i,js+1)+T(2)*R(2,r)*M(i,js+2);
				M(i,js) -= t*R(0,r); M(i,js+1) -= t*R(1,r);
				M(i,js+2) -= t*R(2,r); } } }
	void	updateBlock(const unsigned int il, const unsigned int iu)		//! updates the block X = H(il:iu, il:iu)
		{ const unsigned int bs = iu-il+1;					// block size
		  if (bs == 1) { nr[il] = 1; return; }
		  else if (bs == 2) {							// for block size == 2, do a Givens rotation
			const T m00 = H(il,il)*(H(il,il)-s)+H(il,il+1)*H(il+1,il)+t;
			const T m10 = H(il+1,il)*(H(il,il)+H(il+1,il+1)-s);
			computeReflector(m00,m10,T(0),il);
			doPX(H,il,il+2,il,n,il); doXP(H,0,il+2,il,il+2,il);
			nr[il+1] = 1; return; }
		  const T m00 = H(il,il)*(H(il,il)-s)+H(il,il+1)*H(il+1,il)+t;		// for block size >= 3, use the regular strategy
		  const T m10 = H(il+1,il)*(H(il,il)+H(il+1,il+1)-s);
		  const T m20 = H(il+2,il+1)*H(il+1,il);
		  computeReflector(m00,m10,m20,il);
		  doPX(H,il,il+3,il,n,il); doXP(H,0,il+std::min(bs,4u),il,il+3,il);
		  for (unsigned int i = 1; i < bs-2; i++) {
			const unsigned int is = il+i, js = il+i-1;
			computeReflector(H(is,js),H(is+1,js),H(is+2,js), il+i);
			doPX(H,il+i,il+i+3,il+i-1,n,il+i);
			doXP(H,0,il+std::min(bs,i+4),il+i,il+i+3,il+i); }
		  computeReflector(H(iu-1,iu-2), H(iu,iu-2),T(0),iu-1);
		  doPX(H,iu-1,iu+1,iu-2,n,iu-1); doXP(H,0,il+bs,iu-1,iu+1,iu-1);
		  nr[iu] = 1; }
public:
	deflateDouble(const matD<T>& _H, const T _s, const T _t)			//! allocates a context for a QR decomposition of matrix H with shifts s and t.
		: n(_H.M), H(_H), s(_s), t(_t), R(3,n), nr(n)
		{ std::vector<unsigned int> b; b.push_back(0);
		  for (unsigned int i = 0; i < n-2; i++) { const T h = std::abs(H(i+1,i));
			if (h <= eps || h <= eps*(std::abs(H(i,i))+std::abs(H(i,i+1)))) {
				H(i+1,i) = T(0); b.push_back(i+1); };
			for (unsigned int j = i+2; j < n; j++) H(j,i) = T(0); }
		  b.push_back(n);
		  for (unsigned i = 0; i < b.size()-1; i++)
				updateBlock(b[i],b[i+1]-1); }				// compute reflectors and update each block
	void	YQ(matD<T>& Q)								//! returns decomposition matrix Q.
		{ for (unsigned int i = 0; i < n-2; i++)
			doXP(Q,0,Q.M,i,i+3,i);
		  doXP(Q,0,Q.M,n-2,n,n-2); }
	matD<T>	QtHQ() const								//! returns deflated matrix.
		{ return H; }
};

//! Implements a context to solve a sparse symmetric eigenvalue problem.

template<typename T> class ssev {
	static constexpr T eps = std::numeric_limits<T>::epsilon();
	const T prec = std::pow(eps,T(2.0/3.0));				//!< a number that is approximately zero
	const unsigned int maxit = 10000;
protected:
	const matS<T>& A;			//!< sparse matrix
	const unsigned int n;			//!< dimension of matrix A
	const unsigned int nev;			//!< number of eigenvalues requested
	const unsigned int dim;			//!< number of ritz values
	const unsigned int rule;		//!< sorting rule
	std::vector<vecD<T>> V;			//!< V matrix in the Arnoldi factorization
	matD<T> H;				//!< H matrix in the Arnoldi factorization
	vecD<T> f;				//!< residual in the Arnoldi factorization
	vecD<T> ritzVal;			//!< ritz values
	matD<T> ritzVec;			//!< ritz vectors
	vecD<T> ritzEst;			//!< last row of ritzVec
	bvecD	ritzConv;			//!< indicator of the convergence of ritz values

	void	sortEV(const vecD<T>& m, const matD<T>& M)
		{ SortEigenvalue<T,LARGEST_MAGN> sorting(m.x,m.N);
		  auto ix = sorting.index();
		  switch (rule) {
 		  case LARGEST_MAGN: break;
		  case LARGEST_ALGE: {
			SortEigenvalue<T,LARGEST_ALGE> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  case SMALLEST_MAGN: {
			SortEigenvalue<T,SMALLEST_MAGN> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  case SMALLEST_ALGE: {
			SortEigenvalue<T,SMALLEST_ALGE> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  default: throw rtException("unsupported sorting rule"); }
		  for (unsigned int i = 0; i < m.N; i++) {
            		ritzEst[i] = M(dim-1,ix[i]);
			ritzVal[i] = m(ix[i]); ritzVec.setCol(i,M.getCol(ix[i])); } }
	unsigned int nConv(const T tol)							//! calculates the number of converged Ritz values.
		{ const T fn = norm(f); unsigned int n = 0;
		  for (unsigned int i = 0; i < dim; i++) {
			const T thr = tol*std::max(prec,std::abs(ritzVal(i)));
			const T res = std::abs(ritzEst(i))*fn;
			ritzConv[i] = (res < thr); if (ritzConv(i)) n++; }
		  return n; }
	unsigned int nAdj(const unsigned int nconv)					//! returns the adjusted # of eigenvalues for restarting
		{ unsigned int nev_new = nev;
		  for (unsigned int i = nev; i < dim; i++)
			if (std::abs(ritzEst[i]) < eps) nev_new++;
		  nev_new += std::min(nconv,(dim-nev_new)/2);
		  if (nev_new == 1 && dim >= 6) return dim / 2;
		  else if (nev_new == 1 && dim > 2) return 2;
		  return nev_new > dim-1? dim-1: nev_new; }
	void	getPair()								//! retrieves and sorts ritz values & vectors
		{ const auto D = sev(H); sortEV(D.d,D.U); }				// sort by criterion
	matD<T>	Vmult(const matD<T>& M) const						//! multiplies subspace V by dense matrix M.
		{ vecD<T> c(M.N); matD<T> A(n,c.N);
		  for (unsigned int i = 0; i < n; i++) { c = T(0);
			for (unsigned int j = 0; j < c.N; j++)
				for (unsigned int l = 0; l < M.m(); l++) c[j] += V[l](i)*M(l,j);
			for (unsigned int j = 0; j < c.N; j++) A(i,j) = c[j]; };
		  return A; }
	T	corrResidual(const vecD<T>& q)						//! corrects residual v using new vector.
		{ for (unsigned int i = 0; i < q.N; i++) f -= q(i)*V[i];		// side effects: updates f, returns norm of f.
		  return norm(f); }
	void	init()									//! initializes the Krylov space
		{ vecD<T> v(n); random(v); v /= norm(v); V[0] = v; 
		  vecD<T> w = op(v); H(0,0) = dot(v,w); f = w-v*H(0,0); }		// op is virtual, so should not be called in the constructor		  
	void	factorize(const unsigned int st, const unsigned int end)		//! extends Arnoldi factorization from st to end.
		{ T beta = norm(f); if (end <= st) return;
		  for (unsigned int i = 0; i < H.M; i++)				// keep the upper left k x k submatrix of H and set other elements to 0
			for (unsigned int j = 0; j < H.N; j++)
				if (i >= st || j >= st) H(i,j) = T(0);
		  for (unsigned int i = st; i < end; i++) { bool rs = false;
			if (beta < prec) { random(f); 
				const vecD<T> q = innerProduct(f,i); 			// generate a new residual vector
				beta = corrResidual(q); rs = true; }			// beta <- ||f||
			const vecD<T> v = f/beta, w = op(v); V[i] = v;			// w <- A * v, v = V.col(i)
			H(i,i-1) = rs? T(0): beta; H(i-1,i) = H(i,i-1);
			H(i,i) = innerProduct(v,w);
			f = rs? w-H(i,i)*v: w-H(i,i-1)*V[i-1]-H(i,i)*v;
			beta = norm(f); vecD<T> q = innerProduct(f,i+1);
			for (unsigned int it = 0; it < 5; it++) {
				if (std::abs(q(q.amax())) < prec*beta) break;		// test whether V' * (f/||f||) ~= 0
				H(i-1,i) += q(i-1); H(i,i-1) = H(i-1,i); H(i,i) += q(i);// h <- h + Vf
				beta = corrResidual(q); q = innerProduct(f,i+1); } } }		// correct residual: f <- f - V * Vf
	void	restart(const unsigned int k)						//! restarts Arnoldi factorization
		{ matD<T> Q(dim,dim); Q.id(); if (k >= dim) return;
		  for (unsigned int i = k; i < dim; i++) { 				// shift by all eigenvalues
			deflateTri<T> D(H,ritzVal[i]); D.YQ(Q); H = D.QtHQ(); }		// perform QR shift and revert shift
		  vecD<T> v(n); std::vector<vecD<T>> Vs(k+1);
		  for (unsigned int i = 0; i <= k; i++) {
			const unsigned int nz = i < k? dim-k+i+1: dim; v = T(0);
			for (unsigned int j = 0; j < nz; j++) v += V[j]*Q(j,i);
			Vs[i] = v; }
		  for (unsigned int i = 0; i <= k; i++) V[i] = Vs[i];
		  f = f*Q(dim-1,k-1)+V[k]*H(k,k-1); factorize(k, dim); getPair(); }
	virtual vecD<T> innerProduct(const vecD<T>& v, const unsigned int l) const	//! multiplies transposed subspace V by vector v.
		{ vecD<T> c(l); for (unsigned int i = 0; i < l; i++) c[i] = dot(V[i],v);
		  return c; }
	virtual T innerProduct(const vecD<T>& v, const vecD<T>& w) const		//! multiplies vector v by vector w.
		{ return dot(v,w); }
	virtual T norm(const vecD<T>& v) const						//! returns L2 norm
		{ return ::norm(v); }
	virtual vecD<T> op(const vecD<T> v) const					//! returns matrix-vector product
		{ return A*v; }
	virtual void sortPairs(vecD<T>& val, matD<T>& vec)				//! sorts Ritz pairs by criterion 
		{ (void)val; (void)vec; }
public:
	ssev(const matS<T>& _A, const unsigned int _nev, const unsigned int _dim, const unsigned int _rule)
		//! allocates a context for solving a sparse symmetric eigenvalue problem.
		: A(_A), n(A.M), nev(_nev), dim(_dim > n? n: _dim), rule(_rule), V(dim),
		H(dim,dim), f(n), ritzVal(dim), ritzVec(dim,dim), ritzEst(dim), ritzConv(dim)
		{ H = T(0); ritzVal = T(0); ritzVec = T(0); ritzEst = T(0); ritzConv = false; }
	virtual ~ssev() { }
	retMV<T> solve(const T tol = T(1e-10), const bool verbose = false)
		//! solves the eigenvalue problem; returns eigenvectors and values.
		{ init(); factorize(1,dim); getPair(); unsigned int nc = 0;
		  for (unsigned int it = 0; it < maxit; it++) { nc = nConv(tol); 	// for maxit iterations
			if (verbose) { printf("[%4d] %d\r", it, nc); fflush(stdout); };
			if (nc < nev) restart(nAdj(nc)); else break; }			// if not converged, restart
		  nc = 0; for (unsigned int i = 0; i < dim; i++) if (ritzConv(i)) nc++;
		  nc = std::min(nc,nev); vecD<T> val(nc); matD<T> vec(dim,nc);
		  if (nc == 0) return { vec,val };
		  for (unsigned int i = 0, j = 0; i < dim && j < nc; i++) {		// collect converged pairs
			if (ritzConv(i) == false) continue;
			val[j] = ritzVal[i]; vec.setCol(j,ritzVec.getCol(i)); j++; }
		  sortPairs(val,vec);
		  return { Vmult(vec), val }; }						// transform and return
};

//! Implements a context to solve a sparse symmetric eigenvalue problem using shift-invert transformation.

template <typename T> class ssevs : public ssev<T> {
	T	sigma;
	matS<T> M;
	LDL<T> sv;				//!< direct solver

	vecD<T>	op(const vecD<T> v) const
		{ return sv.solve(v); }							// use LDL decomposition for solve
	void	sortPairs(dvecD& val, dmatD& vec)
		{ for (unsigned int i = 0; i < val.N; i++) val(i) = sigma+T(1)/val(i);
		  this->sortEV(val,vec); }
public:
	ssevs(const matS<T>& A, const unsigned int nev, const unsigned int dim,
		const T _sigma, const unsigned int rule)
		//! allocates a context for solving a sparse symmetric eigenvalue problem using shift sigma.
        	: ssev<T>(A,nev,dim,rule), sigma(_sigma), M(A), sv(M.addToD(-sigma))
		{ }
};

//! Implements a context to solve a sparse general symmetric eigenvalue problem.

template <typename T> class sgsev : public ssev<T> {
	const matS<T>& B;
	LDL<T> sv;				//!< direct solver

	vecD<T>	op(const vecD<T> v) const
		{ vecD<T> t = this->A*v; return sv.solve(t); }							
	vecD<T> innerProduct(const vecD<T>& v, const unsigned int l) const		//! multiplies subspace V by dense vector v.
		{ vecD<T> t = B*v, c(l);
		  for (unsigned int i = 0; i < l; i++) c[i] = dot(this->V[i],t);
		  return c; }
	T	innerProduct(const vecD<T>& v, const vecD<T>& w) const			//! multiplies dense vector v by dense vector Bw.
		{ return dot(v,B*w); }
	virtual T norm(const vecD<T>& v) const						//! returns B norm
		{ return std::sqrt(innerProduct(v,v)); }
public:
	sgsev(const matS<T>& A, const matS<T>& _B, const unsigned int nev,
		const unsigned int dim, const unsigned int rule)
		//! allocates a context for solving a sparse general symmetric eigenvalue problem.
		: ssev<T>(A,nev,dim,rule), B(_B), sv(_B) { }
};

//! Implements a context to solve a sparse non-symmetric eigenvalue problem.

template<typename T> class snev {
	using U = std::complex<T>;
	static constexpr T eps = std::numeric_limits<T>::epsilon();
	const T prec = std::pow(eps,T(2.0/3.0));				//!< a number that is approximately zero
	const unsigned int maxit = 10000;
protected:
	const matS<T>& A;			//!< sparse matrix
	const unsigned int n;			//!< dimension of matrix A
	const unsigned int nv;			//!< number of eigenvalues requested
	const unsigned int dim;			//!< number of ritz values
	const unsigned int rule;		//!< sorting criterion
	std::vector<vecD<T>> V;			//!< V matrix in the Arnoldi factorization
	matD<T> H;				//!< H matrix in the Arnoldi factorization
	vecD<T> f;				//!< residual in the Arnoldi factorization
	vecD<U> ritzVal;			//!< ritz values
	matD<U> ritzVec;			//!< ritz vectors
	vecD<U> ritzEst;			//!< last row of ritzVec
	bvecD	ritzConv;			//!< indicator of the convergence of ritz values

	void	sortEV(const vecD<U>& m, const matD<U>& M)
		{ SortEigenvalue<U,LARGEST_MAGN> sorting(m.x,m.N);
		  auto ix = sorting.index();
		  switch (rule) {
 		  case LARGEST_MAGN: break;
		  case LARGEST_REAL: {
                	SortEigenvalue<U,LARGEST_REAL> sorting(m.x,m.N);
                	ix = sorting.index(); }
			break;
		  case LARGEST_IMAG: {
			SortEigenvalue<U,LARGEST_IMAG> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  case SMALLEST_MAGN: {
			SortEigenvalue<U,SMALLEST_MAGN> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  case SMALLEST_REAL: {
			SortEigenvalue<U,SMALLEST_REAL> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  case SMALLEST_IMAG: {
			SortEigenvalue<U,SMALLEST_IMAG> sorting(m.x,m.N);
			ix = sorting.index(); }
			break;
		  default: throw rtException("unsupported sorting rule"); }
		  for (unsigned int i = 0; i < m.N; i++) {
            		ritzEst[i] = M(dim-1,ix[i]);
			ritzVal[i] = m(ix[i]); ritzVec.setCol(i,M.getCol(ix[i])); } }
	unsigned int nConv(const T tol)							//! calculates the number of converged Ritz values.
		{ const T fn = norm(f); unsigned int n = 0;
		  for (unsigned int i = 0; i < dim; i++) {
			const T thr = tol*std::max(prec,std::abs(ritzVal(i)));
			const T res = std::abs(ritzEst(i))*fn;
			ritzConv[i] = (res < thr); if (ritzConv(i)) n++; }
		  return n; }
	unsigned int nAdj(const unsigned int nconv)					//! returns the adjusted # of eigenvalues for restarting
		{ unsigned int nv_new = nv;
		  for (unsigned int i = nv; i < dim; i++)
			if (std::abs(ritzEst[i]) < eps) nv_new++;
		  nv_new += std::min(nconv,(dim-nv_new)/2);
		  if (nv_new == 1 && dim >= 6) return dim / 2;
		  else if (nv_new == 1 && dim > 3) return 2;
		  if (nv_new > dim-2) nv_new = dim-2;
		  if (isComplex(ritzVal[nv_new-1]) &&
			isConj(ritzVal[nv_new-1], ritzVal[nv_new])) nv_new++;
		  return nv_new; }
	void	getPair()								//! retrieves and sorts ritz values & vectors
		{ const auto D = nev(H); sortEV(D.d,D.U); }				// sort by criterion
	matD<U>	Vmult(const matD<U>& M) const						//! multiplies subspace V by dense matrix M.
		{ vecD<U> c(M.N); matD<U> A(n,c.N);
		  for (unsigned int i = 0; i < n; i++) { c = U(0);
			for (unsigned int j = 0; j < c.N; j++)
				for (unsigned int l = 0; l < M.m(); l++) c[j] += V[l](i)*M(l,j);
			for (unsigned int j = 0; j < c.N; j++) A(i,j) = c[j]; };
		  return A; }
	T	corrResidual(const vecD<T>& q)						//! corrects residual v using new vector.
		{ for (unsigned int i = 0; i < q.N; i++) f -= q(i)*V[i];		// side effects: updates f, returns norm of f.
		  return norm(f); }
	void	init()									//! initializes the Krylov space
		{ vecD<T> v(n); random(v); v /= norm(v); V[0] = v; 
		  vecD<T> w = op(v); H(0,0) = dot(v,w); f = w-v*H(0,0); }		// op is virtual, so should not be called in the constructor		  
	void	factorize(const unsigned int st, const unsigned int end)		//! extends Arnoldi factorization from st to end.
		{ T beta = norm(f); if (end <= st) return;
		  for (unsigned int i = 0; i < H.M; i++)				// keep the upper left k x k submatrix of H and set other elements to 0
			for (unsigned int j = 0; j < H.N; j++)
				if (i >= st || j >= st) H(i,j) = T(0);
		  for (unsigned int i = st; i < end; i++) { bool rs = false;
			if (beta < prec) { random(f); 
				const vecD<T> q = innerProduct(f,i); 			// generate a new residual vector
				beta = corrResidual(q); rs = true; }			// beta <- ||f||
			const vecD<T> v = f/beta; V[i] = v; f = op(v); 			// f <- A * v, v = V.col(i)
			H(i,i-1) = rs? T(0): beta; 
			for (unsigned int j = 0; j <= i; j++) H(j,i) = dot(V[j],f);
			for (unsigned int j = 0; j <= i; j++) f -= V[j]*H(j,i);
			beta = norm(f); if (beta < 0.717*::norm(H.getCol(i))) continue;
			vecD<T> q = innerProduct(f,i+1);
			for (unsigned int it = 0; it < 5; it++) {
				if (std::abs(q(q.amax())) < prec*beta) break;		// test whether V' * (f/||f||) ~= 0
				beta = corrResidual(q); 
				for (unsigned int j = 0; j < q.N; j++) H(j,i) += q(j);
				q = innerProduct(f,i+1); } } }				// correct residual: f <- f - V * Vf
	bool	isComplex(const U& v) const
		{ return v.imag() != T(0); }
	bool	isConj(const U& v1, const U& v2) const
		{ return v1 == conj(v2); }
	void	restart(const unsigned int k)						//! restarts Arnoldi factorization
		{ matD<T> Q(dim,dim); Q.id(); if (k >= dim) return;
		  for (unsigned int i = k; i < dim; i++) { 
			if (isComplex(ritzVal[i]) && isConj(ritzVal[i],ritzVal[i+1])) {
                		const T s = T(2)*ritzVal[i].real(), t = std::abs(ritzVal[i]);
                		deflateDouble<T> ds(H,s,t); ds.YQ(Q); H = ds.QtHQ(); i++; }
			else { deflateSingle<T> ss(H,ritzVal[i].real());
				ss.YQ(Q); H = ss.QtHQ(); } }
		  vecD<T> v(n); std::vector<vecD<T>> Vs(k+1);
		  for (unsigned int i = 0; i <= k; i++) {
			const unsigned int nz = i < k? dim-k+i+1: dim; v = T(0);
			for (unsigned int j = 0; j < nz; j++) v += V[j]*Q(j,i);
			Vs[i] = v; }
		  for (unsigned int i = 0; i <= k; i++) V[i] = Vs[i];
		  f = f*Q(dim-1,k-1)+V[k]*H(k,k-1); factorize(k, dim); getPair(); }
	virtual vecD<T> innerProduct(const vecD<T>& v, const unsigned int l) const	//! multiplies transposed subspace V by vector v.
		{ vecD<T> c(l); for (unsigned int i = 0; i < l; i++) c[i] = dot(V[i],v);
		  return c; }
	virtual T innerProduct(const vecD<T>& v, const vecD<T>& w) const		//! multiplies vector v by vector w.
		{ return dot(v,w); }
	virtual T norm(const vecD<T>& v) const						//! returns L2 norm
		{ return ::norm(v); }
	virtual vecD<T> op(const vecD<T> v)						//! returns matrix-vector product
		{ return A*v; }
	virtual void sortPairs(vecD<U>& val, matD<U>& vec)				//! sorts Ritz pairs by criterion 
		{ (void)val; (void)vec; }
public:
	snev(const matS<T>& _A, const unsigned int _nev, const unsigned int _dim, const unsigned int _rule)
		//! allocates a context for solving a sparse symmetric eigenvalue problem.
		: A(_A), n(A.M), nv(_nev), dim(_dim > n? n: _dim), rule(_rule), V(dim),
		H(dim,dim), f(n), ritzVal(dim), ritzVec(dim,dim), ritzEst(dim), ritzConv(dim)
		{ H = T(0); ritzVal = U(0); ritzVec = U(0); ritzEst = U(0); ritzConv = false; }
	virtual ~snev() { }
	retMcVc<T> solve(const T tol = T(1e-10), const bool verbose = false)
		//! solves the eigenvalue problem; returns eigenvectors and values.
		{ init(); factorize(1,dim); getPair(); unsigned int nc = 0;
		  for (unsigned int it = 0; it < maxit; it++) { nc = nConv(tol); 	// for maxit iterations
			if (verbose) { printf("[%4d] %d\r", it, nc); fflush(stdout); };
			if (nc < nv) restart(nAdj(nc)); else break; }			// if not converged, restart
		  nc = 0; for (unsigned int i = 0; i < dim; i++) if (ritzConv(i)) nc++;
		  nc = std::min(nc,nv); vecD<U> val(nc); matD<U> vec(dim,nc);
		  if (nc == 0) return { vec,val };
		  for (unsigned int i = 0, j = 0; i < dim && j < nc; i++) {		// collect converged pairs
			if (ritzConv(i) == false) continue;
			val[j] = ritzVal[i]; vec.setCol(j,ritzVec.getCol(i)); j++; }
		  sortPairs(val,vec);
		  return { Vmult(vec), val }; }						// transform and return
};

//! Implements a context to solve a sparse non-symmetric eigenvalue problem using shift-invert transformation.

template <typename T> class snevs : public snev<T> {
	using U = std::complex<T>;
	T	sigma;
	matS<T>	M; 
        lsmr<T> sv;				//!< sparse general solver

	vecD<T>	op(const vecD<T> v)
		{ return sv.solve(v,false); }							// for spare matrices use CholDec
	void	sortPairs(vecD<U>& val, matD<U>& vec)
		{ for (unsigned int i = 0; i < val.N; i++) val(i) = U(sigma+T(1)/val(i));
		  this->sortEV(val,vec); }
public:
	snevs(const matS<T>& A, const unsigned int nev, const unsigned int dim,
		const T _sigma, const unsigned int rule)
		//! allocates a context for solving a sparse symmetric eigenvalue problem using shift sigma.
        	: snev<T>(A,nev,dim,rule), sigma(_sigma), M(A),
		sv(M.addToD(-sigma),T(0),T(1e-10))
		{ }
};

//! Implements a context to solve a sparse general non-symmetric eigenvalue problem.

template <typename T> class sgnev : public snev<T> {
	const matS<T>& B;
        lsmr<T> sv;				//!< solver

	vecD<T>	op(const vecD<T> v) const
		{ vecD<T> t = this->A*v; return sv.solve(t); }							
	vecD<T> innerProduct(const vecD<T>& v, const unsigned int l) const		//! multiplies subspace V by dense vector v.
		{ vecD<T> t = B*v, c(l);
		  for (unsigned int i = 0; i < l; i++) c[i] = dot(this->V[i],t);
		  return c; }
	T	innerProduct(const vecD<T>& v, const vecD<T>& w) const			//! multiplies dense vector v by dense vector Bw.
		{ return dot(v,B*w); }
	T	norm(const vecD<T>& v) const						//! returns B norm
		{ return std::sqrt(innerProduct(v,v)); }
public:
	sgnev(const matS<T>& A, const matS<T>& _B, const unsigned int nev,
		const unsigned int dim, const unsigned int rule)
		//! allocates a context for solving a sparse general symmetric eigenvalue problem.
		: snev<T>(A,nev,dim,rule), B(_B), sv(B,T(0),T(1e-10)) { }
};

#endif
