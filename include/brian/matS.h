#ifndef MATS_H
#define MATS_H

/*
 *
 * matS.h: templated sparse matrix class 
 * BRIAN Software Package Version 3.0
 *
 * $Id: matS.h 508 2017-03-26 20:13:21Z frithjof $
 *
 * 0.10 (06/10/09): for BRIAN2 by FK
 * 0.21 (15/05/11): matrix solvers revisited, mv, mvp fixed
 * 0.22 (26/05/11): sparse matrices split off
 * 0.23 (15/06/11): bug fix & test of multiply, add, subtract
 * 0.24 (22/09/11): check init for allocation problems, return bool
 * 0.25 (21/11/11): set row & col revisited
 * 0.26 (03/12/11): switched to CSR format
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.41 (11/08/15): joinBelow, joinRight, nnz, and compact added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements a templated sparse matrix class.
*/

#define unaryOpScalarS(sym) \
	const matS<T>& operator sym (const T b)  { for (unsigned int i = 0; i < nz; i++) x[i] sym b; return *this; }
#define binaryOpScalarS(sym) \
	const matS<T> operator sym (const T b) const \
	{ matS<T> a = *this; for (unsigned int i = 0; i < nz; i++) a.x[i] = x[i] sym b; return a; }

#define allRows(r)	(unsigned int r = 0; r < M; r++)
#define allCols(t)	(unsigned int t = ri[r]; t < ri[r+1]; t++)

//! Helper structure for a matrix element used to create sparse matrices.

template<typename T> class rcv {
public:
	unsigned int r;				//!< row number (absolute)
	unsigned int c;				//!< column number (absolute)
	T	v;				//!< value
	rcv(unsigned int _r = 0, unsigned int _c = 0, T _v = T(0))			//! allocates a matrix element.
		: r(_r), c(_c), v(_v) { };
        bool operator< (const rcv& b) const						//! sort by increasing row.
                { return r < b.r; }
};

//! Implements a sparse matrix.
/*! 
Parameters T and U refer to the data type of the matrix. The default is to have a
single data type, e.g., float or complex<float>. If T != U, a mixture
of real and complex types can be used. Note that not all combinations
of mixed mode operators may be implemented.
*/

template<typename T> class matS {
public:
	unsigned int M;				//!< global # of rows
	unsigned int N;				//!< global # of columns
	unsigned int nz;			//!< local # of non-zero entries
	unsigned int* ci;			//!< col indices, size nz
	unsigned int* ri;			//!< row indices, size M+1
	T* 	x;				//!< array of values, size nz
private:
	void	alloc(const unsigned int _M, const unsigned int _N, const unsigned int _nz)
		//! allocates storage for a sparse matrix with M rows, N columns and nz non-zero elements.
		{ M = _M; N = _N; nz = _nz; ri = new unsigned int [M+1];  
		  ci = new unsigned int [nz]; x = new T [nz]; ri[0] = 0; }
	void	mv(const T *s, T *d) const						//! multiplies matrix by vector.
		{ for allRows(r) { T v{0}; for allCols(t) v += x[t]*s[ci[t]];
			d[r] = v; } }
	void	tmv(const T *s, T *d) const						//! multiplies transpose matrix by vector.
		{ for (unsigned int i = 0; i < N; i++) d[i] = T(0);
		  for allRows(r) { T v = s[r]; for allCols(t) d[ci[t]] += x[t]*v; } }
	void	clear()									//! frees storage and set size to zero.
		{ M = 0; N = 0; nz = 0; delete [] ci; ci = nullptr;
		  delete [] ri; ri = nullptr; delete [] x; x = nullptr; }
	void	assign(const matS<T> &b)						//! (re)allocates and copies from sparse matrix b.
		{ if (M != b.M || N != b.N || nz != b.nz) {
			clear(); alloc(b.M,b.N,b.nz); };				// prune any excess elements
		  if (M) { for (unsigned int i = 0; i <= M; i++) ri[i] = b.ri[i]; }
		  if (nz) { for (unsigned int i = 0; i < nz; i++) {
			ci[i] = b.ci[i]; x[i] = b.x[i]; } } }
	T	elem(const unsigned int r, const unsigned int c) const			//! returns element at coordinates (r,c).
		{ for allCols(t) { if (ci[t] == c) return x[t]; }
//		  throw rtException("invalid element %u %u", r, c); }
		  return T(0); }
	T&	elem(const unsigned int r, const unsigned int c)			//! references element at coordinates (r,c).
		{ for allCols(t) { if (ci[t] == c) return x[t]; }
		  throw rtException("invalid element %u %u", r, c); }
public:
	matS<T>()									//! allocates empty sparse matrix.
		: M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr)
		{ }
	matS<T>(const unsigned int _M, const unsigned int _N, const unsigned int _nz)	//! allocates sparse matrix with M rows, N columns and nz non-zero elements.
		: M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr)
		{ alloc(_M,_N,_nz); }
	matS<T>(const matS<T>& b)							//! constructs sparse matrix from matrix b.
		: M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr)
		{ assign(b); }
	matS<T>(matS<T>&& b)								//! moves from sparse matrix b.
		: M(b.M), N(b.N), nz(b.nz), ci(b.ci), ri(b.ri), x(b.x)
		{ b.ci = nullptr; b.ri = nullptr; b.x = nullptr; }
	matS<T>(std::vector<rcv<T>>& v)							//! constructs sparse matrix from vector of rcv elements.
		: M(0), N(0), nz(0), ci(nullptr), ri(nullptr), x(nullptr)
		{ std::sort(v.begin(),v.end());	N = v.back().r+1; alloc(N,N,v.size());	// sort rcv table and allocate matrix
		  unsigned int r = 0; for (unsigned int t = 0; t < v.size(); t++) {	// for all elements
			if (v[t].r == r) { ci[t] = v[t].c; x[t] = v[t].v; t++; }	// add element to current row
			else { for (unsigned int i = r+1; i <= v[t].r; i++) ri[i] = t;	// else advance to new row
				r = v[t].r; t--; } };
		  ri[N+1] = v.size(); }
	virtual ~matS<T>() { clear(); }							// destructor
	matS<T>& operator=(const matS<T>& b)						//! assigns from matrix b.
		{ if (this != &b) assign(b); return *this; }
	matS<T>& operator=(matS<T>&& b)							//! move assigns from matrix b.
		{ assert(this != &b); M = b.M; N = b.N; nz = b.nz; 
		  delete [] ci; ci = b.ci; b.ci = nullptr;
		  delete [] ri; ri = b.ri; b.ri = nullptr;
		  delete [] x; x = b.x; b.x = nullptr; return *this; }
	unsigned int m() const								//! returns number of rows in this matrix.
		{ return M; }
	unsigned int n() const								//! returns number of columns in this matrix.
		{ return N; }
        T	operator()(const unsigned int i, const unsigned int j) const		//! returns element at coordinates (i,j).
		{ return elem(i, j); }
	T&	operator()(const unsigned int i, const unsigned int j)			//! returns reference to element at coordinates (i,j).
		{ return elem(i, j); }
		unaryOpScalarS(=)
		unaryOpScalarS(+=)
		unaryOpScalarS(-=)
		unaryOpScalarS(*=)
		unaryOpScalarS(/=)
		binaryOpScalarS(+)
		binaryOpScalarS(-)
		binaryOpScalarS(*)
		binaryOpScalarS(/)
	vecD<T>	operator*(const vecD<T>& b) const					//! multiplies this matrix by dense vector b.
		{ assert(N == b.N); vecD<T> a(M); mv(b.x,a.x); return a; }
	matD<T>	operator*(const matD<T>& B) const					//! multiplies this matrix by dense matrix b.
		{ assert(N == B.M); matD<T> A(M,B.N); 
		  for (unsigned int i = 0; i < B.N; i++) mv(B.x+i*N, A.x+i*M);
		  return A; }
	vecS<T>	operator*(const vecS<T>& b) const					//! multiplies this matrix by sparse vector b.
		{ vecD<T> s = asDense(b,N), d(N); tmv(s.x,d.x); return asSparse(d); }
	vecD<T>	tm(const vecD<T>& b) const						//! multiplies transpose matrix by dense vector b.
		{ assert(M == b.N); vecD<T> a(N); tmv(b.x,a.x); return a; }
	vecS<T>	tm(const vecS<T>& b) const						//! multiplies transpose matrix by sparse vector b.
		{ vecD<T> s = asDense(b,N), d(N); mv(s.x,d.x); return asSparse(d); }
	vecD<T>	idia() const								//! returns inverse of diagonal as a dense vector.
		{ vecD<T> a(M); for allRows(i) a[i] = T(1)/elem(i,i); return a; }
	vecD<T>	d() const								//! returns diagonal as a dense vector.
		{ vecD<T> a(M); for allRows(i) a[i] = elem(i,i); return a; }
	matS<T>& id()									//! sets sparse matrix to identity.
		{ unsigned int n = N; clear(); alloc(n, n, n);
		  for (unsigned int i = 0; i < n; i++) {
			ci[i] = i; ri[i] = i; x[i] = T(1); };
		  ri[n] = n; return *this; }
	matS<T>& addToD(const vecD<T>& b)						//! adds dense vector b to diagonal.
		{ assert(M == N && N == b.N);
		  for allRows(r) for allCols(t) {
			if (ci[t] == r) { x[t] += b(r); break; } };
		  return *this; }
	matS<T>& addToD(const T b)							//! adds scalar b to diagonal.
		{ assert(M == N); for allRows(r) for allCols(t) { 
			if (ci[t] == r) { x[t] += b; break; } };
		  return *this; }
	matS<T>& multiplyRows(const vecD<T>& b)						//! multiplies this matrix by diagonal matrix represented as vector b.
		{ assert(M == b.N); for allRows(r) {
			const T f = b(r); for allCols(t) x[t] *= f; };
		  return *this; }
	matS<T>& divideRows(const vecD<T>& b)						//! divides this matrix by diagonal matrix represented as vector b.
		{ assert(M == b.N); for allRows(r) {
			const T f = b(r); for allCols(t) x[t] /= f; };
		  return *this; }
	matS<T>& resize(const unsigned int M, const unsigned int N, const unsigned int nz)
		//! resizes storage of this matrix for M rows, N columns and nz non-zero elements.
		{ clear(); alloc(M,N,nz); return *this; }
	bool	finite() const								//! returns true if all elements are finite.
		{ for (unsigned int i = 0; i < nz; i++) {
			if (std::isfinite(x[i]) == false) return false; };
		  return true; }
	matS<T>	operator*(const matS<T>& B) const; 					// binary matS<T>-matS<T> multiplication
	matS<T> operator+(const matS<T>& B) const; 					// binary matS<T>-matS<T> addition
	matS<T> operator-(const matS<T>& B) const; 					// binary matS<T>-matS<T> subtraction
	matS<T>& setRow(const unsigned int r, const vecS<T>& v);			// sets values in a row
	matS<T>& setCol(const unsigned int c, const vecS<T>& v);			// sets values in a column
	matS<T>& initRow(const unsigned int r, const vecS<T>& v);			// initializes a row - must be called sequentially
	matS<T>& sort();
	matS<T>	strength(const T th) const;
	ivecD	aggregate() const;
	T	spectralRadius(const unsigned int k = 10, const unsigned int l = 20) const;
	matS<T>& joinBelow(const matD<T>& A);
	matS<T>& joinRight(const matD<T>& A);
	unsigned int nnz(const T lim = T(0)) const					//! returns number of elements greater than lim.
		{ unsigned int n = 0;
		  for (unsigned int i = 0; i < nz; i++) n += std::abs(x[i]) > lim;
		  return n; }
	matS<T>	prune(const T lim = T(0)) const;
	T	sum() const
		{ T s(0); for (unsigned int i = 0; i < nz; i++) s += x[i]; return s; }
};

//! \relates matS
template <typename T>
matS<T> matS<T>::operator*(const matS<T>& B) const
//! multiplies this matrix by sparse matrix B.
{	assert(N == B.M); ivecD nx(B.N); nx = -1; unsigned int nnz = 0;
	for allRows(r) { for allCols(t) { const unsigned int c = ci[t];
			for (unsigned int u = B.ri[c]; u < B.ri[c+1]; u++) {
				const unsigned int k = B.ci[u];
				if (nx[k] != int(r)) { nx[k] = int(r); nnz++; } } } };
	matS<T> C(M, B.N, nnz); vecD<T> s(B.N); s = T(0); nx = -1; nnz = 0;
	for allRows(r) { int h = -2; unsigned int n = 0; 
		for allCols(t) { unsigned int c = ci[t]; T v = x[t];
			for (unsigned int u = B.ri[c]; u < B.ri[c+1]; u++) {
				unsigned int k = B.ci[u]; s[k] += v * B.x[u];
				if (nx[k] == -1) { nx[k] = h; h = int(k); n++; } } };
		for (unsigned int i = 0; i < n; i++) {
			const unsigned int k = ITOU(h);
			if (std::abs(s[k]) > 0) { C.ci[nnz] = k; C.x[nnz] = s[k]; nnz++; };
			h = nx[k]; nx[k] = -1; s[k] = T(0); };                             
		C.ri[r+1] = nnz; };
	C.nz = nnz; return C;								// excess space will be removed during assignment.
}

//! \relates matS
template <typename T>
matS<T> matS<T>::operator+(const matS<T>& B) const
//! adds sparse matrix B to this matrix.
{	assert(M == B.M && N == B.N); matS<T> C(M,N,nz+B.nz);
	ivecD nx(M); vecD<T> a(N), b(N); a = T(0); b = T(0); nx = -1;
	unsigned int u = 0;
	for allRows(r) { int h = -2; unsigned int n = 0;
		for allCols(t) { unsigned int c = ci[t]; a[c] = x[t];
			if (nx[c] == -1) { nx[c] = h; h = int(c); n++; } };
		for (unsigned int t = B.ri[r]; t < B.ri[r+1]; t++) {
			const unsigned int c = B.ci[t]; b[c] = B.x[t];
			if (nx[c] == -1) { nx[c] = h; h = int(c); n++; } };
		for (unsigned int t = 0; t < n && h >= 0; t++) {
			const T v = a[h]+b[h]; C.ci[u] = h; C.x[u] = v; u++;
			const unsigned int k = ITOU(h);
			h = nx[k]; nx[k] = -1; a[k] = 0; b[k] = 0; };
		C.ri[r+1] = u; }
	C.sort(); return C;								// excess space will be removed during assignment.
}

//! \relates matS
template <typename T>
matS<T> matS<T>::operator-(const matS<T>& B) const
//! subtracts sparse matrix B from this matrix.
{	assert(M == B.M && N == B.N); matS<T> C(M,N,nz+B.nz);
	ivecD nx(M); vecD<T> a(N), b(N); a = T(0); b = T(0); nx = -1;
	unsigned int u = 0;
	for allRows(r) { int h = -2; unsigned int n = 0;
		for allCols(t) { unsigned int c = ci[t]; a[c] = x[t];
			if (nx[c] == -1) { nx[c] = h; h = int(c); n++; } };
		for (unsigned int t = B.ri[r]; t < B.ri[r+1]; t++) {
			const unsigned int c = B.ci[t]; b[c] = B.x[t];
			if (nx[c] == -1) { nx[c] = h; h = int(c); n++; } };
		for (unsigned int t = 0; t < n && h >= 0; t++) {
			const T v = a[h]-b[h]; C.ci[u] = h; C.x[u] = v; u++;
			const unsigned int k = ITOU(h);
			h = nx[k]; nx[k] = -1; a[k] = 0; b[k] = 0; };
		C.ri[r+1] = u; }
	C.sort(); return C;								// excess space will be removed during assignment.
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::initRow(const unsigned int r, const vecS<T>& v)
//! inits row r from sparse vector v - must be called sequentially by increasing r.
{	unsigned int t = ri[r];
	for (unsigned int i = 0; i < v.N; i++) { ci[t] = v.u[i]; x[t] = v.x[i]; t++; };
	ri[r+1] = t; return *this;
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::setRow(const unsigned int r, const vecS<T>& v)
//! sets row r from sparse vector v without updating column index.
{	assert(r < M); unsigned int i = 0;
	for allCols(t) { assert(ci[t] == v.u[i]); x[t] = v.x[i]; i++; };
	return *this;
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::setCol(const unsigned int c, const vecS<T>& v)
//! sets column c from sparse vector v without updating column index.
{	assert(c < N); for (unsigned int i = 0; i < v.N; i++) {
		const unsigned int r = v.u[i];
		for allCols(t) { if (ci[t] == c) { x[t] = v.x[i]; break; } } };
	return *this;
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::sort()
//! sorts columns by row indices.
{	for allRows(r) { for allCols(t) {						// for all entries in this column
		unsigned int c = ci[t]; T v = x[t]; int i = int(t-1);			// perform an insertion sort
		while (i >= int(ri[r]) && ci[i] > c) {
			ci[i+1] = ci[i]; x[i+1] = x[i]; i--;
			ci[i+1] = c; x[i+1] = v; } } };
	return *this;
}

//! \relates matS
template <typename T>
matS<T> matS<T>::strength(const T th) const
//! computes strength S given association limit th.
{	const vecD<T> dia = d()*th; unsigned int nz = 0;
	for allRows(r) for allCols(t) {
		if (norm(x[t]) >= std::abs(dia(ci[t])*dia(r))) nz++; };
	matS<T> S(M,N,nz); unsigned int u = 0;
	for allRows(r) { S.ri[r] = u; for allCols(t) {
		if (norm(x[t]) >= std::abs(dia(ci[t])*dia(r))) {
			S.ci[u] = ci[t]; S.x[u] = x[t]; u++; } } };
	S.ri[M] = u; return S;
}

template <typename T>
ivecD matS<T>::aggregate() const
//! computes aggregation vector agg from a strength matrix.
{	ivecD agg(M); agg = 0; int na = 1;
	for allRows(r) { if (agg[r]) continue; bool ha = false, hn = false;		// already marked
		for allCols(t) { const unsigned int c = ci[t]; if (r == c) continue;	// check if all neighbors are free
			hn = true; if (agg[c]) { ha = true; break; } };
		if (!hn) agg[r] = -int(M);						// isolated node, do not aggregate
		else if (!ha) {	agg[r] = na; for allCols(t) agg[ci[t]] = na; na++; } };	// make an aggregate of this node and its neighbors
	for allRows(r) { if (agg[r]) continue;						// add unaggregated nodes to any neighboring aggregate
		for allCols(t) { const int a = agg[ci[t]];				// check if all neighbors are free
            		if (a > 0) { agg[r] = -a; break; } } }
	na--; 
	for allRows(r) { const int a = agg[r];
		if (a != 0) {								// node r has been aggregated
			if (a > 0) agg[r] = a-1;
			else if (a == -int(M)) agg[r] = -1;
			else agg[r] = -a-1;
			continue; }
		agg[r] = na;								// else add to next aggregate
		for allCols(t) { if (agg[ci[t]] == 0) agg[ci[t]] = na; }
		na++; }
	return agg;
}

template <typename T>
T matS<T>::spectralRadius(const unsigned int k, const unsigned int l) const
//! returns largest eigenvalue of matrix A.
{	const unsigned int n = std::min(M,k); matD<T> H_(n+1,n); H_ = T(0);
	std::vector<vecD<T>> v(n+1); vecD<T> t(M); random(t); v[0] = t/T(norm(t));
	for (unsigned int j = 0; j < n; j++) { v[j+1] = *this*v[j];
		for (unsigned int i = 0; i <= j; i++)
			{ H_(i,j) = dot(v[i],v[j+1]); v[j+1] -= H_(i,j)*v[i]; }
		H_(j+1,j) = norm(v[j+1]); v[j+1] /= H_(j+1,j); };
	if (n == 0) return 0;
	matD<T> H(n,n); vecD<T> x(n), y; random(x);
	for (unsigned int i = 0; i < n; i++)
		for (unsigned int j = 0; j < n; j++) H(i,j) = H_(i,j);
	for (unsigned int i = 0; i < l; i++) {
		T m{0}; unsigned int u = 0; 
		for (unsigned int j = 0; j < n; j++) {
			if (std::abs(x[j]) > std::abs(m)) { m = x[j]; u = j; } };
		T f = T(1)/x[u]; x *= f; y = x; x = H*x; };
	return l? norm(x)/norm(y): T(0); 
}

//! \relates matS
template <typename T, typename U>
vecD<U>	operator*(const matS<T>& A, const vecD<U>& b)
//! multiplies real sparse matrix A with dense complex vector b.
{      assert(A.N == b.N); vecD<U> a(A.M); for (unsigned int r = 0; r < A.M; r++) {
		U s = U(0); for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) s += U(A.x[t],0)*b(A.ci[t]);
                a[r] = s; }; return a;
}

//! \relates matS
template <typename T, typename U>
matD<U> operator*(const matS<T>& A, const matD<U>& B)
//! multiplies real sparse matrix A with dense complex matrix B.
{	assert(A.N == B.M); matD<U> C(A.M,B.N); 
	for (unsigned int i = 0; i < B.N; i++) {
		for (unsigned int r = 0; r < A.M; r++) { U s = U(0); U *v = B.x+i*A.N;
                for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) s += A.x[t]*v[A.ci[t]];
                C(r,i) = s; } }; return C;
}

//! \relates matS
template <typename T>
void print(const matS<T>& A)
//! prints real sparse matrix A on stdout.
{	for (unsigned int r = 0; r < A.M; r++) {
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
			printf("%5d %5d ", r, A.ci[t]); printT(A.x[t]); printf("\n"); } };
}

//! \relates matS
template <typename T> 
matS<std::complex<T>> asComplex(const matS<T>& A)
//! converts real sparse matrix A to complex sparse matrix.
{	matS<std::complex<T>> B(A.M,A.N,A.nz); for (unsigned int r = 0; r < A.M; r++)
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
			B.ci[t] = A.ci[t]; B.x[t] = std::complex<T>(A.x[t],0.0); }; 
	for (unsigned int r = 0; r <= A.M; r++) B.ri[r] = A.ri[r];
	return B;
}

//! \relates matS
template <typename T>
matD<T> asDense(const matS<T>& A)
//! converts sparse matrix A to dense matrix.
{	matD<T> D(A.M,A.N); D = T(0); for (unsigned int r = 0; r < A.M; r++)
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) D.x[A.ci[t]*A.M+r] = A.x[t];
	return D;
}

//! \relates matS
template <typename T>
matS<T> asSparse(const matD<T>& A)
//! converts dense matrix A to sparse matrix.
{	matS<T> S(A.M,A.N,A.nnz()); S = T(0); unsigned int t = 0; S.ri[0] = 0;
	for (unsigned int r = 0; r < A.M; r++) {
		for (unsigned int c = 0; c < A.N; c++) { if (std::abs(A(r,c)) == T(0)) continue;
			S.x[t] = A(r,c); S.ci[t] = c; t++; }; S.ri[r+1] = t; };
	return S;
}

//! \relates matS
template<typename T>
matS<T> trp(const matS<T>& A)
//! returns transpose of a real sparse matrix A.
{	uvecD cnt(A.N); cnt = 0; for (unsigned int t = 0; t < A.nz; t++) cnt[A.ci[t]]++;
	matS<T> B(A.N, A.M, A.nz);
	for (unsigned int r = 0; r < A.N; r++) { B.ri[r+1] = cnt[r]+B.ri[r]; cnt[r] = B.ri[r]; };
	for (unsigned int r = 0; r < A.M; r++) {
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) { 
			unsigned int k = cnt[A.ci[t]]++; B.ci[k] = r; B.x[k] = A.x[t]; } };
	return B; }

//! \relates matS
template <typename T> 
matS<T> adj(const matS<T> &A)
//! returns adjoint of a real sparse matrix A.
{	return trp(A); }
			
//! \relates matS
template <typename T> 
matS<std::complex<T>> adj(const matS<std::complex<T>>& A)
//! returns adjoint of a complex sparse matrix A.
{ 	uvecD cnt(A.N); cnt = 0; for (unsigned int t = 0; t < A.nz; t++) cnt[A.ci[t]]++;
	matS<std::complex<T>> B(A.N, A.M, A.nz);
	for (unsigned int r = 0; r < A.N; r++) { B.ri[r+1] = cnt[r]+B.ri[r]; cnt[r] = B.ri[r]; };
	for (unsigned int r = 0; r < A.M; r++) {
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) { 
			unsigned int k = cnt[A.ci[t]]++; B.ci[k] = r; B.x[k] = conj(A.x[t]); } };
	return B;
}

//! \relates matS
template <typename T>
void gs(matD<T>& W, const matS<T>& M, const matD<T>& Q)
//! performs Gram-Schmidt orthonormalization of W against basis Q.
{	for (unsigned int i = 0; i < W.N; i++) { vecD<T> w = W.getCol(i);
		T nm = norm(w), ni = nm == 0? 0: T(1.0)/nm; w *= ni; vecD<T> mw = M*w;
		for (unsigned int j = 0; j < Q.N; j++) {
			const vecD<T> v = Q.getCol(j); w -= v*dot(v,mw); };
		nm *= norm(w); W.setCol(i,w*nm); };
}

//! \relates matS
template <typename T>
void ortho(matD<T>& W, const matS<T>& M, const matD<T>& Q)
//! performs SQVB orthonormalization of W against basis Q.
{	for (unsigned int it = 0; it < 100; it++) {
		gs(W,M,Q); gs(W,M,Q);							// perform two rounds of GS orthogonalization
		int c = svqb(W,M*W); if (c == 0) return;				// orthonormalize W
		if (c == -1) random(W); }						// breakdown in svqb, reinit W
}

//! \relates matS
template <typename T>
void read(const char* fname, matS<T>& A)
//! reads a real sparse matrix from file in matrix market format.
{	FILE *fp = openFile(fname,"r");
	unsigned int state = 0, M, N, nz, r, c, t = 0; T v; char line[LINELEN];
	while (fgets(line,LINELEN,fp)) {						// first pass
		if (line[0] == '%') continue;						// comment
		switch (state) {
		case 0: { const int n = sscanf(line, "%u %u %u", &M, &N, &nz);		// dimensions
			if (n != 3 || M == 0 || N == 0 || nz < M || nz < N)
				throw optException("Illegal size line: %s",line);
			A.resize(M,N,nz); for (r = 0; r <= M; r++) A.ri[r] = 0; 
			state = 1; }; break;
		case 1: { const int n = sscanf(line, "%u %u %lg\n", &r, &c, &v);	// element entry
			if (n != 3 || r > M || c > N)
				throw optException("Illegal element line: %s", line);
			A.ri[r]++; t++; }; break;
		default: break; } };
	if (t != nz) throw optException("Element lines (%u) to not match size (%u).\n", int(t), int(nz)); 
	for (r = 1; r <= M; r++) A.ri[r] += A.ri[r-1];					// build row indices 
	assert(A.ri[M] == nz);
	state = 0; rewind(fp); ivecD rc(M); rc = 0;
	while (fgets(line,LINELEN,fp)) {						// second pass
		if (line[0] == '%') continue;						// comment
		switch (state) {
		case 0: { const int n = sscanf(line, "%u %u %u", &M, &N, &nz);		// dimensions
			if (n != 3 || M == 0 || N == 0 || nz < M || nz < N)
				throw optException("Illegal size line: %s.\n",line);
			state = 1; }; break;
		case 1: { const int n = sscanf(line, "%u %u %lg\n", &r, &c, &v);	// element entry
			if (n != 3 || r > M || c > N)
				throw optException("Illegal element line: %s.\n", line);
			r--; c--; t = A.ri[r]+rc[r]; assert(t < nz);			// correct for 0-based
			A.ci[t] = c; A.x[t] = v; rc[r]++; }; break;
		default: break; } };
	closeFile(fp);
}

//! \relates matS
template <typename T>
void save(const char* fname, const matS<T>& A)
//! saves a real sparse matrix into a file in matrix market format.
{	FILE *fp = openFile(fname,"w");
	unsigned int n = 0;
	for (unsigned int t = 0; t < A.nz; t++) if (std::abs(A.x[t]) > 1e-12) n++;
	fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
	fprintf(fp, "%6d %6d %6d\n", A.M, A.N, n);
	for (unsigned int r = 0; r < A.M; r++)
		for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
			T x = A.x[t]; if (std::abs(x) < 1e-12) continue;
			fprintf(fp, "%6d %6d %e\n", r+1, A.ci[t]+1, x); }
	closeFile(fp);
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::joinBelow(const matD<T>& A)
//! joins a dense matrix A below this matrix.
{	assert(A.N == N); unsigned int n = A.nnz(), t = nz;				// find number of nonzero elements in A
	unsigned int* nri = new unsigned int [M+A.M+1];					// reallocate this matrix for A.M additional rows
	memcpy(nri,ri,(M+1)*sizeof(unsigned int)); delete [] ri; ri = nri;		// copy old matrix
	unsigned int* nci = new unsigned int [n+t];
	memcpy(nci,ci,t*sizeof(unsigned int)); delete [] ci; ci = nci;
	T* nx = new T [n+t];
	memcpy(nx,x,t*sizeof(T)); delete [] x; x = nx;
	for (unsigned int r = 0; r < A.M; r++) {					// now fill in A
		for (unsigned int c = 0; c < A.N; c++) {
			T v = A(r,c); if (v == T(0)) continue;				// skip zero elements
			x[t] = v; ci[t] = c; t++; }					// add element to *this
		ri[M+r+1] = t; };
	assert(t == n+nz); M += A.M; return *this;					// set new row size
}

//! \relates matS
template <typename T>
matS<T>& matS<T>::joinRight(const matD<T>& A)
//! joins a dense matrix A at the right side of this matrix.
{	assert(A.M == M); unsigned int n = A.nnz(), u = 0;				// find number of nonzero elements in A
	unsigned int* nri = new unsigned int [M+1];					// reallocate this matrix for A.N new columns
	unsigned int* nci = new unsigned int [n+nz];
	T* nx = new T [n+nz]; nri[0] = 0;
	for (unsigned int r = 0; r < M; r++) {						// for all rows
		for (unsigned int t = ri[r]; t < ri[r+1]; t++, u++) {			// copy current sparse row
			nx[u] = x[t]; nci[u] = ci[t]; }
		for (unsigned int c = 0; c < A.N; c++) {				// add new elements from A
			T v = A(r,c); if (v == T(0)) continue;				// skip zero elements
			nx[u] = v; nci[u] = c+N; u++; }					// add element to current row
		nri[r+1] = u; };
	assert(u == n+nz);								// check fill complete
	delete [] ri; ri = nri; delete [] ci; ci = nci; N += A.N;			// deallocate old matrix
	delete [] x; x = nx; nz = n+nz;	return *this;					// and update dimensions
}

//! \relates matS
template <typename T>
matS<T> matS<T>::prune(const T lim) const
//! removes entries less or equal than lim from matrix.
{	const unsigned int n = nnz(lim); matS<T> A(M,N,n);
	unsigned int u = 0; A.ri[0] = u;
	for allRows(r) { for allCols(t) { if (std::abs(x[t]) <= lim) continue;		// skip small values
		A.ci[u] = ci[t]; A.x[u] = x[t]; u++; }; A.ri[r+1] = u; };
	return A;
}

//! \relates matS
template <typename T>
T largestEV(const matS<T>& A, const T eps = 1e-6)
//! determines largest eigenvalue of A using Rayleigh iterations.
{	T lambda{0}; vecD<T> x(A.M); random(x);
	for (unsigned int it = 0; it < 10000; it++) {
		const vecD<T> u = A*x; T un = norm(u);
		const T l = dot(u,x)/dot(x,x), d = l/lambda-T(1);
		x = u/un; lambda = l;  if (d < eps) break; };
	return lambda;
}

//! \relates matS
template <typename T>
T cond(const matS<T>& A, const T eps = 1e-6, const bool verbose = true)
//! computes the condition number of matrix A.
// see: Avorn H, Durinsky A, Toledo S (2015).
// Spectral-Condition-Number Estimation of Large Sparse Matrices. arXiv:1301-1107v3.
{	const T smax = largestEV(A); T smin = smax;					// step 2
	vecD<T> x(A.M); random(x); const T xn = norm(x); const vecD<T> xs = x/xn;	// step 6,8
	const vecD<T> b = A*xs; T beta = norm(b); const T den = xn*smax+beta;		// step 9
	vecD<T> u = b/beta, v = A*u; T alpha = norm(v); v /= alpha;			// step 10,11
	vecD<T> w = v; x = T(0); T phib = beta, rhob = alpha;				// steps 12-14
	for (unsigned int it = 0; it < 10000; it++) {					// step 16
		u = A*v-alpha*u; beta = norm(u); u /= beta;				// steps 17-19
		v = A.tm(u)-beta*v; alpha = norm(v); v /= alpha;			// steps 20-22
		const T rho = std::sqrt(SQR(rhob)+SQR(beta));				// step 23
		const T c = rhob/rho, s = beta/rho;					// step 24
		const T theta = s*alpha; rhob = -c*alpha;				// steps 25-26
		const T phi = c*phib; phib *= s;
		x += (phi/rho)*w; w = v-(theta/rho)*w; const vecD<T> d = xs-x;		// steps 27-28
		const T dn = norm(d); if (isSmall<T>(dn)) return T(1);			// steps 31-32
		const T Adn = norm(A*d), sn = Adn/dn, cv = Adn/den;			// steps 33-35
		if (sn <= smin) smin = sn;
		if (verbose) { printf("[%3d] cond %e err %e\r", it, smax/smin, cv); fflush(stdout); };
		if (cv < eps || smax/smin >= T(64)/eps) break; };
	if (verbose) { printf("cond %-40e\n", smax/smin); };
	return smax/smin;	
}

//! \relates matS
template <typename T>
matS<T>	symmetrize(const matS<T>& A)
//! returns 0.5*(A^T+A).
{	return (A+trp(A))*T(0.5);
}

using fmatS = matS<float>;
using dmatS = matS<double>;
using cmatS = matS<fcomplex>;
using zmatS = matS<dcomplex>;

#define allACols(t,r) (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)
#define allAtCols(t,r) (unsigned int t = At.ri[r]; t < At.ri[r+1]; t++)

//! Helper class for expression A d A<sup>T</sup> for a sparse symmetric matrix A and a dense vector d.

template <typename T>
class AdAt : public matS<T> {
public:
	AdAt<T>(const matS<T>& A, const matS<T>& At);
	void	updateD(const matS<T>& A, const matS<T>& At, const vecD<T>& D);
};

template <typename T>
AdAt<T>::AdAt(const matS<T>& A, const matS<T>& At)
//! allocates structures for matrix A*d*A<SUP>T</SUP>. Uses lower triangle only.
{	this->M = A.m(); this->N = this->M; this->nz = 0; ivecD done(this->M); done = -1; 
	for (unsigned int r = 0; r < this->M; r++) for allACols(t,r) { 
		for allAtCols(u,A.ci[t]) { unsigned int c = At.ci[u]; if (c > r) continue;
			if (done[c] != int(r)) { done[c] = int(r); this->nz++; } } };
	this->resize(this->M, this->N, this->nz); this->nz = 0;
	for (unsigned int r = 0; r < this->M; r++) { this->ri[r] = this->nz; done[r] = -1; 
		for allACols(t,r) for allAtCols(u,A.ci[t]) {
			unsigned int c = At.ci[u]; if (c > r) continue;
			if (done[c] != int(r)) { done[c] = int(r); this->ci[this->nz++] = c; } } };
	this->ri[this->N] = this->nz;
}

template <typename T>
void AdAt<T>::updateD(const matS<T>& A, const matS<T>& At, const vecD<T>& D)
//! updates diagonal in A*d*A<sup>T</sup>. Assumes A,At constant between calls. Uses lower triangle only.
{	vecD<T> w(A.N), v(At.N);
	for (unsigned int r = 0; r < this->M; r++) { v = 0;
		for allACols(t,r) { unsigned int c = A.ci[t]; w[c] = A.x[t]*D(c); }
		for allACols(t,r) { unsigned int c = A.ci[t]; 
			for allAtCols(u,c) v[At.ci[u]] += w[c]*At.x[u]; };
		for (unsigned int t = this->ri[r]; t < this->ri[r+1]; t++) {
			unsigned int c = this->ci[t]; if (c <= r) this->x[t] = v[c]; } };
}

#undef allRows
#undef allCols
#undef allARows
#undef allAtRows
#undef unaryOpScalarS
#undef binaryOpScalarS

using fAdAt = AdAt<float>;
using dAdAt = AdAt<double>;

#endif

