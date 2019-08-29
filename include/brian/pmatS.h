#ifndef PMATS_H
#define PMATS_H

/*
 *
 * pmatS.h: partitioned vector and matrix classes
 * BRIAN Software Package Version 3.0
 *
 * $Id: pmatS.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (08/12/11): for BRIAN2.2 by FK
 * 0.20 (16/12/13): documented
 * 0.30 (02/11/14): revised
 * 0.50 (17/11/14): released for BRIAN2.7
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements partitioned vector, matrix and solver classes.
*/

#define allElem(i) (unsigned int i = 0; i < n; i++) 

#define unaryOpScalarPV(sym) \
	pvecD<T>& operator sym (const T b) \
	{ const unsigned int n = p.len(); \
	  for allElem(i) x[i] sym b; return *this; }
#define binaryOpScalarPV(sym) \
	pvecD<T> operator sym (const T b) const \
	{ const unsigned int n = p.len(); \
	  pvecD<T> a(p); for allElem(i) a[i] = x[i] sym b; return a; }
#define unaryOpVecPV(sym) \
	pvecD<T>& operator sym (const pvecD<T> &b) \
	{ const unsigned int n = p.len(); assert(n == b.p.len()); \
	  for allElem(i) x[i] sym b.x[i]; return *this; }
#define binaryOpVecPV(sym) \
	pvecD<T> operator sym (const pvecD<T> &b) const \
	{ const unsigned int n = p.len(); assert(n == b.p.len()); \
	  pvecD<T> a(p); for allElem(i) a[i] = x[i] sym b.x[i]; return a; }

//! Implements a partitioned dense vector.

template<typename T> class pvecD {
	void	alloc()									//! allocates space.
		{ N = p.tlen(); x = new T [N]; }
	void	clear()									//! frees up space.
		{ delete [] x; x = nullptr; N = 0; }
	void	assign(const pvecD<T>& b)						//! assign from partitioned vector b.
		{ if (N != b.N) { clear(); alloc(); };
		  for (unsigned int i = 0; i < p.tlen(); i++) x[i] = b(i); }
public:
	unsigned int N;				//!< number of local elements
	T*	x;				//!< array of elements
	partition p;				//!< partitioning info

	pvecD<T>(): N(0), x(nullptr), p()						//! allocates an empty partitioned vector.
		{ }			
	pvecD<T>(const partition& _p)							//! allocates a vector for partition p.
		 : N(0), x(nullptr), p(_p) { alloc(); }
	pvecD<T>(const pvecD<T>& b)							//! copies from vector b.
		 : N(0), x(nullptr), p(b.p) { assign(b); }
	pvecD<T>(pvecD<T>&& b)								//! copies from vector b.
		 : N(b.N), x(b.x), p(b.p) { b.x = nullptr; }
	~pvecD<T>()
		{ clear(); }
	pvecD<T>& operator=(const pvecD<T>& b)						//! assigns from vector b.
		{ if (this != &b) { p = b.p; assign(b); }; return *this; }
	pvecD<T>& operator=(pvecD<T>&& b)						//! move assigns from vector b.
		{ assert(this != &b); delete [] x; x = b.x; b.x = nullptr;
		  N = b.N; p = std::move(b.p); return *this; }
        T&	operator[](const unsigned int i)					//! returns LHS reference to element i.
		{ return x[i]; }
        T	operator()(const unsigned int i) const					//! returns RHS element i.
		{ return x[i]; }
        T&	operator()(const unsigned int i)					//! returns RHS reference to element i.
		{ return x[i]; }
	unsigned int n() const								//! returns length of partitioned vector.
		{ return p.dim(); };
	void	sync() const								//! synchronizes this vector across all partitions.
		{ p.sync(x); }
		unaryOpScalarPV(=)
		unaryOpScalarPV(+=)
		unaryOpScalarPV(-=)
		unaryOpScalarPV(*=)
		unaryOpScalarPV(/=)
		binaryOpScalarPV(+)
		binaryOpScalarPV(-)
		binaryOpScalarPV(*)
		binaryOpScalarPV(/)
		unaryOpVecPV(+=)
		unaryOpVecPV(-=)
		unaryOpVecPV(*=)
		unaryOpVecPV(/=)
		binaryOpVecPV(+)
		binaryOpVecPV(-)
		binaryOpVecPV(*)
		binaryOpVecPV(/)
	void	repartition(const partition& _p)					//! resizes using a new partition.
		{ const unsigned int l = p.len(); assert(l == 0 || l >= _p.len());	// my partition must be zero or match exactly
		  const T *a = x; x = nullptr; clear(); p = _p; alloc(); 		// save my elements and re-allocate
		  for (unsigned int i = 0; i < l; i++) x[i] = a[i];			// copy my elements to new buffer and clean up
		  delete [] a; }
	vecD<T>	gather() const								//! gathers partitioned vector into dense vector at the root process.
		{ const unsigned int np = p.nprocs, n = p.dim(); vecD<T> a(n);
		  if (np == 1) {
			for (unsigned int i = 0; i < n; i++) a[i] = x[i]; }		// single process: just copy x to a
        	  else { ivecD sz(np); for (unsigned int i = 0; i < np; i++)
				sz[i] = p.off(i+1)-p.off(i);
                	mpiGather(x, p.len(), a.x, sz.x, p.off.x);
			if (p.isRoot() == false) {
				a.N = 0; delete [] a.x; a.x = nullptr; } };
		 return a; }
	void	scatter(const vecD<T>& a)						//! scatters this vector to all processes.
		{ const unsigned int np = p.nprocs;
		  if (np == 1) {
			for (unsigned int i = 0; i < a.N; i++) x[i] = a(i); }		// single process: just copy a to x
        	  else { ivecD sz(np); for (unsigned int i = 0; i < np; i++)
				sz[i] = p.off(i+1)-p.off(i);
                	mpiScatter(a.x, sz.x, p.off.x, x, p.len()); } }
	T	max() const;
};

//! \relates pvecD
template <typename T>
T pvecD<T>::max() const
//! returns maximum across all partitions in all partitions.
{	T m = std::numeric_limits<T>::min(); const unsigned int np = p.nprocs;
	for (unsigned int i = 0; i < p.len(); i++) m = std::max(m,x[i]);
	vecD<T> mall(np); mall = 0; ivecD sz(np), off(np); sz = 1;
	for (unsigned int i = 0; i < np; i++) off[i] = i;
	mpiAllgather(&m, 1, mall.x, sz.x, off.x);
	m = std::numeric_limits<T>::min();
	for (unsigned int i = 0; i < np; i++) m = std::max(m,mall(i));
	return m;
}

//! \relates pvecD
template <typename T>
pvecD<T> operator*(const T a, pvecD<T>& b)
//! multiplies scalar a with real vector b.
{	return b*a;
}

//! \relates pvecD
template <typename T>
T norm2(const pvecD<T>& b)
//! returns infinite norm of real vector b.
{	return dot(b,b);
}

//! \relates pvecD
template <typename T>
T norm(const pvecD<T>& b)
//! returns norm of real vector b.
{	return std::sqrt(norm2(b));
}

//! \relates pvecD
template <typename T> 
T dot(const pvecD<T>& a, const pvecD<T>& b)
//! returns dot product of real vectors (a,b).
{	assert(a.N == b.N); const unsigned int n = a.p.len();
	T c = 0; for (unsigned int i = 0; i < n; i++) c += a(i)*b(i);
	return mpiSum(c, 1);
}

//! \relates pvecD
template <typename T>
void print(const pvecD<T>& a)
//! prints real vector b to stdout.
{	const unsigned int s = a.p.st(), l = a.p.len();
	for (unsigned int i = 0; i < l; i++) {
		printf("%5d ", i+s); printT(a.x[i]); printf("\n"); }
}

//! \relates pvecD
template <typename T>
vecD<T>	asDense(const pvecD<T>& a)
//! gathers partitioned vector into dense vector.
{	const unsigned int np = a.p.nprocs, n = a.p.dim(); vecD<T> b(n);
	if (np == 1) {
		for (unsigned int i = 0; i < n; i++) b[i] = a(i); }			// single process: just copy x to a
        else {	ivecD sz(np);
		for (unsigned int i = 0; i < np; i++)
			sz[i] = a.p.off(i+1)-a.p.off(i);
		mpiAllgather(a.x, a.p.len(), b.x, sz.x, a.p.off.x); };
	return b;
}

#undef allElem

using pivecD = pvecD<int>;
using puvecD = pvecD<unsigned int>;
using pfvecD = pvecD<float>;
using pdvecD = pvecD<double>;

#define unaryOpScalarPM(sym) \
	pmatS<T>& operator sym (const T b) { for (unsigned int i = 0; i < nz; i++) x[i] sym b; return *this; }

#define allRows(r)	(unsigned int r = 0; r < p.len(); r++)
#define allCols(t)	(unsigned int t = ri[r]; t < ri[r+1]; t++)

//! Implements a partitioned sparse matrix.

template<typename T> class pmatS {
public:
	unsigned int	M;			//!< global # of rows
	unsigned int	N;			//!< global # of columns
	unsigned int	nz;			//!< local # of non-zero entries
	unsigned int*	ci;			//!< col indices, size nz
	unsigned int*	ri;			//!< row indices, size M+1
	T		*x;			//!< array of values, size nz
	partition	p;			//!< partition info
	bool		assembled;		//!< whether this matrix' element configuration was frozen.
private:
	void	alloc(const unsigned int _M, const unsigned int _N,
			const unsigned int _nz)						//! allocates space for a matrix partition.
		{ M = _M; N = _N; nz = _nz;
		  ri = new unsigned int [p.len()+1];
		  ci = new unsigned int [nz]; x = new T [nz]; ri[0] = 0; }
	void	clear()  								//! frees up space.
		{ M = 0; N = 0; nz = 0; delete [] ci; ci = nullptr;
		  delete [] ri; ri = nullptr; delete [] x; x = nullptr; }
	void	assign(const pmatS<T>& B)						//! assigns from partitioned sparse matrix B.
		{ if (M != B.M || N != B.N || nz != B.nz || compatible(p,B.p) == false) {
			clear(); p = B.p; alloc(B.M,B.N,B.nz); };
		  p = B.p; assembled = B.assembled; 
		  if (M) { for allRows(i) ri[i] = B.ri[i]; ri[p.len()] = B.nz; };
		  if (nz) { for (unsigned int i = 0; i < nz; i++) {
				ci[i] = B.ci[i]; x[i] = B.x[i]; } } }
	void	mv(const T* s, T* d) const						//! performs matrix-vector multiplication.
		{ for allRows(r) { T v = 0; for allCols(t) v += x[t]*s[ci[t]]; d[r] = v; } }
	void	mvm(const T* s, T* d) const						//! performs matrix-vector multiplication.
		{ for allRows(r) { T v = 0; for allCols(t) v += x[t]*s[l2g(t)]; d[r] = v; } }
	T	elem(const unsigned int r, const unsigned int c) const			//! returns element (r,c) in local coordinates.
		{ for (unsigned int t = ri[r]; t < ri[r+1]; t++)
			if (ci[t] == c) return x[t];
		  throw rtException("No element %u %u\n", r, c); }
	T&	elem(const unsigned int r, const unsigned int c)			//! returns reference to element (r,c) in local coordinates.
		{ for (unsigned int t = ri[r]; t < ri[r+1]; t++)
			if (ci[t] == c) return x[t];
		  throw rtException("No element %u %u\n", r, c); }
	void	multABt(std::vector<rcv<T>>& elC, const std::vector<rcv<T>>& elBt) const;
public:
	pmatS<T>()									//! allocates an empty partitioned sparse matrix.
		 : M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr), p(), assembled(false)
		{ }
	pmatS<T>(const unsigned int _M, const unsigned int _N, const unsigned int _nz)	//! allocates a partitioned sparse matrix.
		 : M(_M), N(_N), nz(_nz), ci(nullptr), ri(nullptr), x(nullptr), p(M), assembled(false)
		{ alloc( _M, _N, _nz); }
	pmatS<T>(const unsigned int _M, const unsigned int _N, const unsigned int _nz,
			const partition& pt)						//! allocates a sparse matrix for partition pt.
		 : M(_M), N(_N), nz(_nz), ci(nullptr), ri(nullptr), x(nullptr), p(pt), assembled(false)
		{ alloc( _M, _N, _nz); p.initialized = false; }
	pmatS<T>(const pmatS<T>& B) 							//! copies from partitioned sparse matrix B.
		 : M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr), p(), assembled(false)
		{ assign(B); }
	pmatS<T>(pmatS<T>&& B) 								//! moves from partitioned sparse matrix B.
		 : M(B.M), N(B.N), nz(B.nz), ci(B.ci), ri(B.ri), x(B.x), p(B.p), assembled(B.assembled)
		{ B.ci = nullptr; B.ri = nullptr; B.x = nullptr; }
	~pmatS<T>()									//! destructs this partitioned matrix.
		{ clear(); }
	pmatS<T>& operator=(const pmatS<T>& b)						//! assigns from sparse matrix b.
		{ if (this != &b) assign(b); return *this; }
	pmatS<T>& operator=(pmatS<T>&& b)						//! assigns from sparse matrix b.
		{ assert(this != &b);
		  M = b.M; N = b.N; nz = b.nz; p = std::move(b.p); assembled = b.assembled;
		  delete [] ci; ci = b.ci; b.ci = nullptr;
		  delete [] ri; ri = b.ci; b.ri = nullptr;
		  delete [] x; x = b.x; b.x = nullptr; return *this; }
	unsigned int m() const								//! returns rows of full matrix.
		{ return M; }
	unsigned int n() const								//! returns cols of full matrix.
		{ return N; }
        T	operator()(const unsigned int i, const unsigned int j) const		//! returns element (i,j) in global coordinates.
		{ return elem(i, j); }
	T&	operator()(const unsigned int i, const unsigned int j)			//! returns references to element (i,j) in global coordinates.
		{ return elem(i, j); }
	pvecD<T> operator*(const pvecD<T>& b) const					//! multiplies matrix by partitioned dense vector b.
		{ assert(assembled);
		  pvecD<T> a(p); b.sync(); mv(b.x, a.x); return a; }
	pvecD<T> operator*(const vecD<T>& b) const					//! multiplies matrix by dense vector b.
		{ assert(N == b.N); assert(assembled);
		  pvecD<T> a(p); mvm(b.x, a.x); return a; }
	pvecD<T> idia() const								//! returns inverse of diagonal.
		{ assert(assembled); pvecD<T> a(p);
		  for allRows(i) a[i] = T(1.0)/elem(i,i);
		  return a; }
	pvecD<T> d() const								//! returns diagonal.
		{ assert(assembled); pvecD<T> a(p);
		  for allRows(i) a[i] = elem(i,i);
		  return a; }
	vecS<T>	operator*(const vecS<T>& b) const					//! multiplies matrix by sparse vector b.
		{ assert(assembled); vecS<T> w;
		  for (unsigned int i = 0; i < b.N; i++) {
			const T v = b.x[i]; const unsigned int r = b.u[i];
			for allCols(t) {
				const unsigned int c = ci[t]; if (c >= p.len()) continue;
				w.update(c, x[t]*v); } };
		  return w; }
		unaryOpScalarPM(=)
  		unaryOpScalarPM(*=)
  		unaryOpScalarPM(/=)
	unsigned int l2g(const unsigned int t) const					//! maps column index at t to absolute column.
		{ const unsigned int c = ci[t], s = p.st(), l = p.len();
		  return c < l? c+s: p.imi(c-l); }
	uvecD	l2gmap() const								//! constructs a map from global to local column indices.
		{ const unsigned int s = p.st(), l = p.len(); uvecD u(l+p.imi.N); 
		  for (unsigned int i = 0; i < l; i++) u[i] = i+s;
		  for (unsigned int i = 0; i < p.imi.N; i++) u[i+l] = p.imi(i);
		  return u; }
	uvecD	g2lmap() const								//! constructs a map from global to local column indices.
		{ const unsigned int s = p.st(), l = p.len(); uvecD u(N); u = std::numeric_limits<unsigned int>::max(); 
		  for (unsigned int i = 0; i < l; i++) u[i+s] = i;
		  for (unsigned int i = 0; i < p.imi.N; i++) u[p.imi(i)] = l+i;
		  return u; }
	void	resize(const unsigned int _M, const unsigned int _N,
			const unsigned int _nz, const partition _p)			//! resize storage of this matrix for M rows, N columns and nz non-zero elements.
		{ clear(); p = _p; alloc(_M,_N,_nz); }
	void	initRow(const unsigned int r, const vecS<T>& v);
	void	sort() const;
	void	assemble();
	matS<T>	gather() const;
	pmatS<T> strength(const T th) const;
	pivecD	aggregate();
	T	spectralRadius(const unsigned int k = 10, const unsigned int l = 20) const;
	pmatS<T> operator+(const pmatS<T>& B) const;
	pmatS<T> operator-(const pmatS<T>& B) const;
	pmatS<T> operator*(const pmatS<T>& B) const;
};

//! \relates pmatS
template <typename T>
void pmatS<T>::initRow(const unsigned int r, const vecS<T>& v)
//! initializes a row - must be called sequentially.
{	assert(assembled == false); unsigned int t = ri[r-p.st()];			// ensure initRow is used before assembly
	for (unsigned int i = 0; i < v.N; i++) { 
		ci[t] = v.u[i]; x[t] = v.x[i]; t++; };
	ri[r-p.st()+1] = t; 
}

//! \relates pmatS
template <typename T>
void pmatS<T>::sort() const
//! sorts columns by row indices.
{	assert(assembled == false);							// will not work on assembled mat
	for allRows(r) { for allCols(t) {						// for all entries in this column
		unsigned int c = ci[t]; T v = x[t]; int i = int(t)-1;			// perform an insertion sort
		while (i >= (int)ri[r] && ci[i] > c) {
			ci[i+1] = ci[i]; x[i+1] = x[i]; i--; ci[i+1] = c; x[i+1] = v; } } };
}

//! \relates pmatS
template <typename T>
pmatS<T> pmatS<T>::strength(const T th) const
//! computes strength S given association limit th.
{	const vecD<T> dia = asDense(d()*th);						// get the whole diagonal
	unsigned int s = p.st(), nz = 0;
	for allRows(r) for allCols(t) {							// count non-zero elements
		const unsigned int c = l2g(t);
		if (norm(x[t]) >= std::abs(dia(c)*dia(r+s))) nz++; };
	pmatS<T> S(M,N,nz,p); unsigned int u = 0;					// allocate partitioned matrix
	for allRows(r) { S.ri[r] = u; for allCols(t) {
		const unsigned int c = l2g(t);
		if (norm(x[t]) >= std::abs(dia(c)*dia(r+s))) {
			S.ci[u] = c; S.x[u] = x[t]; u++; } } };				// set non-zero elements
	S.ri[p.len()] = u; S.assemble(); return S;
}

//! \relates pmatS
template <typename T>
pivecD pmatS<T>::aggregate()
//! computes aggregation vector agg from strength matrix S.
{	// on output, agg[i] corresponds to the aggregate # of vertex i
	matS<T> S = gather(); ivecD agg(M);
	if (p.isRoot()) { S.sort(); agg = S.aggregate(); };				// do this on root only
	pivecD pAgg(p); pAgg.scatter(agg); return pAgg;
}

//! \relates pmatS
template <typename T>
void pmatS<T>::assemble()
//! converts rows to offset 0 and initializes import index in partition.
{	if (assembled) return;								// don't do this twice
	sort(); uvecD ind(N); ind = 0; unsigned int ni = 0;				// sort rows by increasing column#
	for allRows(r) { for allCols(t) { unsigned int c = ci[t];			// find all import indices
		if (p.isLocal(c) || ind(c)) continue;
		ni++; ind[c] = 1; } };
	uvecD imp(ni); unsigned int j = 0;						// compact into import array
	for (unsigned int i = 0; i < N; i++) { if (ind(i)) imp[j++] = i; };
	for (unsigned int i = 0; i < ni; i++) { ind[imp[i]] = p.len()+i; };
	for allRows(r) { for allCols(t) { unsigned int c = ci[t];			// map column indices
		ci[t] = p.isLocal(c)? c-p.st(): ind(c); } };
	p.init(imp); assembled = true;							// build communication structures
}
	
//! \relates pmatS
template <typename T>
matS<T> pmatS<T>::gather() const
//! gathers a partitioned sparse matrix into sparse matrix at the root process.
{	assert(assembled); const unsigned int np = p.nprocs;
	if (np == 1) { matS<T> A(M,N,nz);
		for (unsigned int i = 0; i <= M; i++) A.ri[i] = ri[i];
		for (unsigned int i = 0; i < nz; i++) { A.ci[i] = ci[i]; A.x[i] = x[i]; }
		return A; }
	std::vector<rcv<T>> el[np]; 
	const unsigned int s = p.st(), pi = p.id(); el[pi].clear();
	for allRows(r) for allCols(t) { el[pi].push_back({r+s,l2g(t),x[t]}); };		// convert my partition into (r,c,v) format
	if (pi) { mpiSendBuf(&el[pi][0], el[pi].size()*sizeof(rcv<T>), 0);		// all non-root processes send their partitions to root
		return matS<T>(M,N,0); }						// and return an empty matrix
	for (unsigned int i = 1; i < np; i++) { unsigned int len = 0; 			// ...while root receives all partitions
		rcv<T>* e = reinterpret_cast<rcv<T>*>(mpiRecvBuf(len, i));		// ...while all others receive the rcv table
		el[i].assign(e,e+len/sizeof(rcv<T>)); delete [] e; };
	unsigned int n = 0; for (unsigned int i = 0; i < np; i++) n += el[i].size();	// count number of non-zero elements
	matS<T> A(M,N,n); unsigned int r = 0, t = 0; A.ri[r] = 0;			// assemble elements into new matrix
	for (unsigned int i = 0; i < np; i++) { 					// for all rcv packages
		for (unsigned int j = 0; j < el[i].size(); j++, t++) {
			rcv<T>& e = el[i][j]; A.ci[t] = e.c; A.x[t] = e.v;	 	// add element into matrix
			if (e.r != r) { r = e.r; A.ri[r] = t; } } };
	A.ri[M] = n; return A;								// root returns full matrix
}

//! \relates pmatS
template <typename T>
T pmatS<T>::spectralRadius(const unsigned int k, const unsigned int l) const
//! computes ratio of eigenvalues.
{	assert(assembled); unsigned int N = M, n = std::min(N,k); uniformDist<T> ud;
	matD<T> H_(n+1, n); H_ = T(0);  std::vector<pvecD<T>> v(n+1); pvecD<T> t(p);
	for (unsigned int i = 0; i < p.len(); i++) t[i] = ud();
	T f = T(1)/norm(t); v[0] = t*f;
	for (unsigned int j = 0; j < n; j++) { v[j+1] = *this*v[j]; 
		for (unsigned int i = 0; i <= j; i++) {
			H_(i,j) = dot(v[i],v[j+1]); v[j+1] -= H_(i,j)*v[i]; }
		H_(j+1,j) = norm(v[j+1]);  v[j+1] /= H_(j+1,j); };
	if (n == 0) return T(0);
	matD<T> H(n,n); vecD<T> x(n), y;
	for (unsigned int i = 0; i < n; i++)
		for (unsigned int j = 0; j < n; j++) H(i,j) = H_(i,j);
	for (unsigned int i = 0; i < n; i++) x[i] = ud();
	for (unsigned int i = 0; i < l; i++) { f = T(1)/x[x.imax()]; x *= f; y = x; x = H*x; };
	return l? norm(x)/norm(y): T(0); 
}

//! \relates pmatS
template <typename T>
pmatS<T> pmatS<T>::operator+(const pmatS<T>& B) const
//! adds sparse matrix B to this matrix. Assumes column numbers are sorted.
{	assert(assembled && B.assembled); assert(M == B.M && N == B.N);			// matrix dimensions must match
	assert(compatible(p,B.p));							// partitions must be compatible
	if (nz == 0) return B; else if (B.nz == 0) return *this;			// this is 0+B or A+0
	const unsigned int nr = p.len(); uvecD _ri(nr); _ri = 0;
	for (unsigned int r = 0; r < nr; r++) {
		unsigned int k1, k2, i1 = 0, i2 = 0; if (r > 0) _ri[r] = _ri[r-1]; 
		while (1) {
			if (i1 >= ri[r+1]-ri[r] && i2 >= B.ri[r+1]-B.ri[r]) break;
			if (i1 == ri[r+1]-ri[r]) k1 = N; else k1 = ci[ri[r]+i1]; 
			if (i2 == B.ri[r+1]-B.ri[r]) k2=N; else k2 = B.ci[B.ri[r]+i1];
			if (k1 < k2) { i1++; _ri[r]++; }
			else if (k1 > k2) { i2++; _ri[r]++; }
			else if (k1 == k2) { i1++; i2++; _ri[r]++; } } };
	unsigned int nnz = _ri[nr-1];
	pmatS<T> C(M, N, nnz, p); int k = 0; C.ri[0] = 0;
	for (unsigned int r = 0; r < nr; r++) C.ri[r+1] = _ri[r]; 
	for (unsigned int r = 0; r < nr; r++) {
		unsigned int k1, k2, i1 = 0, i2 = 0; 
		while (1) {
			if (i1 >= ri[r+1]-ri[r] && i2 >= B.ri[r+1]-B.ri[r]) break;
			if (i1 == ri[r+1]-ri[r]) k1 = N; else k1 = ci[ri[r]+i1];
			if (i2 == B.ri[r+1]-B.ri[r]) k2 = N; else k2 = B.ci[B.ri[r]+i1];
			if (k1 < k2) {
				C.ci[k] = ci[ri[r]+i1];
				C.x[k] = x[ri[r]+i1]; i1++; k++; }
			else if (k1 > k2) {
				C.ci[k] = B.ci[B.ri[r]+i2];
				C.x[k] = B.x[B.ri[r]+i2]; i2++; k++; }
			else if (k1 == k2) {
				C.ci[k] = ci[ri[r]+i1];
				C.x[k] = x[ri[r]+i1]+B.x[B.ri[r]+i2];
				i2++; i1++; k++; } } }
	C.assembled = assembled; return C;
}

//! \relates pmatS
template <typename T>
pmatS<T> pmatS<T>::operator-(const pmatS<T>& B) const
//! subtracts sparse matrix B from this matrix. Assumes column numbers are sorted.
{	assert(assembled && B.assembled); assert(M == B.M && N == B.N);			// matrix dimensions must match
	assert(compatible(p,B.p));							// partitions must be compatible
	if (nz == 0) { pmatS<T> C = B; C *= -1; return C; }				// this is 0-B
	else if (B.nz == 0) return *this;						// this is A-0
	const unsigned int nr = p.len(); uvecD _ri(nr); _ri = 0;
	for (unsigned int r = 0; r < nr; r++) {
		unsigned int k1, k2, i1 = 0, i2 = 0; if (r > 0) _ri[r] = _ri[r-1]; 
		while (1) {
			if (i1 >= ri[r+1]-ri[r] && i2 >= B.ri[r+1]-B.ri[r]) break;
			if (i1 == ri[r+1]-ri[r]) k1 = N; else k1 = ci[ri[r]+i1]; 
			if (i2 == B.ri[r+1]-B.ri[r]) k2=N; else k2 = B.ci[B.ri[r]+i1];
			if (k1 < k2) { i1++; _ri[r]++; }
			else if (k1 > k2) { i2++; _ri[r]++; }
			else if (k1 == k2) { i1++; i2++; _ri[r]++; } } };
	unsigned int nnz = _ri[nr-1];
	pmatS<T> C(M, N, nnz, p); int k = 0; C.ri[0] = 0;
	for (unsigned int r = 0; r < nr; r++) C.ri[r+1] = _ri[r]; 
	for (unsigned int r = 0; r < nr; r++) {
		unsigned int k1, k2, i1 = 0, i2 = 0; 
		while (1) {
			if (i1 >= ri[r+1]-ri[r] && i2 >= B.ri[r+1]-B.ri[r]) break;
			if (i1 == ri[r+1]-ri[r]) k1 = N; else k1 = ci[ri[r]+i1];
			if (i2 == B.ri[r+1]-B.ri[r]) k2 = N; else k2 = B.ci[B.ri[r]+i1];
			if (k1 < k2) {
				C.ci[k] = ci[ri[r]+i1];
				C.x[k] = x[ri[r]+i1]; i1++; k++; }
			else if (k1 > k2) {
				C.ci[k] = B.ci[B.ri[r]+i2];
				C.x[k] = -B.x[B.ri[r]+i2]; i2++; k++; }
			else if (k1 == k2) {
				C.ci[k] = ci[ri[r]+i1];
				C.x[k] = x[ri[r]+i1]-B.x[B.ri[r]+i2];
				i2++; i1++; k++; } } }
	C.assembled = assembled; return C;
}

//! \relates pmatS
template <typename T>
void pmatS<T>::multABt(std::vector<rcv<T>>& elC, const std::vector<rcv<T>>& elBt) const
//! multiplies this partition with transposed rcv table elBt and add results to rcv table elC.
{	const unsigned int s = p.st(), n = elBt.size(); uvecD map = l2gmap(); vecD<T> row(N);
	for (unsigned int i = 0; i < n; ) { row = T(0);                       		// for all elements in B
		rcv<T> e = elBt[i]; const unsigned int rb = e.r;
		while (i < n && e.r == rb) { row[e.c] = e.v; e = elBt[++i]; };		// unpack into row vector
		for allRows(r) { T v = T(0);
			for allCols(t) v += x[t]*row(map(ci[t]));			// compute dot product with row of A
                        if (isSmall<T>(v) == false) elC.push_back({r+s,rb,v}); } };	// save non-zero elements
}

//! \relates pmatS
template <typename T>
pmatS<T> pmatS<T>::operator*(const pmatS<T>& B) const
//! multiplies this matrix by the partitioned sparse matrix B.
{	assert(assembled && B.assembled); assert(N == B.M);				// matrix dimensions must match
	const unsigned int np = B.p.nprocs, pid = B.p.id();
	if (nz == 0 || B.nz == 0) { pmatS<T> C(M,B.N,0); C.assemble(); return C; };	// this is A*0 or 0*B
	std::vector<rcv<T>> elC;							// result in (r,c,v) format
	{ pmatS<T> Bt = trp(B); Bt.sort(); std::vector<rcv<T>> elBt;			// start a block here
	  const unsigned int s = Bt.p.st(), l = Bt.p.len();				// convert partition of Bt into (r,c,v) format
	  for (unsigned int r = 0; r < l; r++) {
		for (unsigned int t = Bt.ri[r]; t < Bt.ri[r+1]; t++) {			// for all elements of this partition
			elBt.push_back({r+s,Bt.l2g(t),Bt.x[t]}); } };
	  multABt(elC,elBt); 
	  if (np > 1) {									// now exchange matrix blocks among all partitions
		vecD<MPI_Request> req(np-1); unsigned int j = 0;
		for (unsigned int i = 0; i < np; i++) {	if (pid == i) continue;		// To all processes except myself...
			unsigned int len = elBt.size()*sizeof(rcv<T>);			// pack my Bt partition
			char* bp = reinterpret_cast<char*>(&elBt[0]);
			mpiSend(bp, len, i, &req[j++]); };				// and send it out
		for (unsigned int i = 0; i < np; i++) {	if (pid == i) continue;		// From all processes except myself...
			unsigned int len = 0; std::vector<rcv<T>> elBn;
			rcv<T>* e = reinterpret_cast<rcv<T>*>(mpiRecvBuf(len, i));	// ...I wait to receive their Bt partition
			elBn.assign(e,e+len/sizeof(rcv<T>)); 				// save as rcv table elBn
			delete [] e; multABt(elC,elBn); }				// multiply this matrix with elBn, results in elC
		mpiWait(j, req.x); } };							// wait for send operations to complete
	std::sort(elC.begin(),elC.end()); const unsigned int n = elC.size();		// sort rcv table
	pmatS<T> C(M,B.N,n,p); const unsigned int s = p.st(), l = p.len();		// assemble elements into new matrix
	unsigned int t = 0; C.ri[0] = 0;
	for (unsigned int r = 0; r < l; r++) {  					// for all rows
		while (t < n) { rcv<T>& e = elC[t]; if (e.r != r+s) break;		// while elC has elements for this row
			C.ci[t] = e.c; C.x[t] = e.v; t++; };				// add element to C
		C.ri[r+1] = t; };							// at end: set row
	C.assemble(); return C;								// build up new communication structure
}

//! \relates pmatS
template <typename T>
pmatS<T> trp(const pmatS<T>& A)
//! transposes a partitioned sparse matrix.
{	assert(A.assembled); partition p(A.N);
	const unsigned int np = A.p.nprocs, nc = 2*np, pid = A.p.id();
	if (A.nz == 0) { pmatS<T> B(A.N,B.M,0); B.assemble(); return B; };
	std::vector<rcv<T>> el[nc]; unsigned int s = A.p.st(), l = A.p.len();		// convert partition into transposed (r,c,v) format
	for (unsigned int r = 0; r < A.p.len(); r++) {
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {			// for all elements of this partition
			const unsigned int c = A.l2g(t);				// re-construct column index
			const unsigned int u = p.partitionOf(c);			// find destination partition
			el[u].push_back({c,r+s,A.x[t]}); } };				// note that we transpose indices here
	if (np > 1) { for (unsigned int i = 0; i < np; i++) {				// now exchange matrix blocks among all partitions
		if (pid == i) { 							// I send...
			for (unsigned int j = 0; j < np; j++) { if (j == pid) continue;	// to all processes except myself...
				unsigned int len = el[j].size()*sizeof(rcv<T>);
				mpiSendBuf(&el[j][0], len, j); } }			// a message containing the rcv table
		else {	unsigned int len = 0;						// ...while all other processes...
			rcv<T>* e = reinterpret_cast<rcv<T>*>(mpiRecvBuf(len, i));	// receive the rcv table
			el[np+i].assign(e,e+len/sizeof(rcv<T>)); delete [] e; } };	// and save table for compiling the transpose
		for (unsigned int i = 0; i < np; i++) {	if (i != pid) el[i].clear(); } }; // clear tables sent
	for (unsigned int i = 0; i < nc; i++) std::sort(el[i].begin(),el[i].end());	// sort all matrix entries by increasing row #
	unsigned int n = 0; for (unsigned int i = 0; i < nc; i++) n += el[i].size();	// count number of non-zero elements
	pmatS<T> B(A.N,A.M,n); uvecD rp(nc); rp = 0; 					// allocate transpose
	s = B.p.st(), l = B.p.len(); unsigned int t = 0; B.ri[0] = 0;			// assemble elements into new matrix
	for (unsigned int r = s; r < s+l; r++) {					// for all rows of the new partition
		for (unsigned int i = 0; i < nc; i++) { unsigned int j = rp[i];		// for all rcv packages
			for ( ; j < el[i].size(); j++) { rcv<T>& e = el[i][j]; 
				assert(e.r >= r); if (e.r > r) break; 			// check row number
				B.ci[t] = e.c; B.x[t] = e.v; t++; }; rp[i] = j; }; 	// add element into matrix
		B.ri[r-s+1] = t; };
	B.assemble(); return B;								// build up new communication structure
}

//! \relates pmatS
template <typename T>
void print(const pmatS<T>& A)
//! prints real partitioned matrix A to stdout.
{	assert(A.assembled); const unsigned int s = A.p.st(), l = A.p.len();
	for (unsigned int r = 0; r < l; r++)
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
			printf("%5d %5d ", r+s, A.l2g(t)); printT(A.x[t]); printf("\n"); }
}

//! \relates pmatS
template <typename T>
matD<T> asDense(const pmatS<T>& A)
//! converts partitioned sparse matrix A to dense matrix at root.
{	matS<T> B = A.gather();								// gather partitioned to root
	if (A.p.isRoot() == false) return matD<T>();					// all non-root processes return here
	matD<T> C(B.M,B.N); C = T(0); for (unsigned int r = 0; r < B.M; r++)
		for (unsigned int t = B.ri[r]; t < B.ri[r+1]; t++) C(r,B.ci[t]) = B.x[t];
	return C;									// root returns dense matrix here
}

#undef allRows
#undef allCols
#undef unaryOpScalarPM

using pimatS = pmatS<int>;
using pumatS = pmatS<unsigned int>;
using pfmatS = pmatS<float>;
using pdmatS = pmatS<double>;

#endif
