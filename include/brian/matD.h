#ifndef MATD_H
#define MATD_H

/*
 *
 * matD.h: templated classes for 3d and dense matrices
 * BRIAN Software Package Version 3.0
 *
 * $Id: matD.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (06/10/09): for BRIAN2 by FK
 * 0.20 (15/05/11): matrix solvers revisited, mv, mvp fixed
 * 0.21 (26/05/11): sparse matrices split off
 * 0.22 (10/06/12): reworked interface&  added complex classes
 * 0.23 (23/06/12): matD tm and mcc corrected
 * 0.30 (31/12/12): released version 2.4
 * 0.31 (18/05/13): matrix solve added
 * 0.32 (24/05/13): matrix norm2 bugfix
 * 0.40 (19/12/13): documented
 * 0.41 (04/03/14): pinv bugfix for non-symmetric matrices
 * 0.42 (06/09/14): removeColumn() added
 * 0.43 (25/09/14): pinv bounds corrected for non-symmetric matrices
 * 0.50 (06/12/14): Cholesky decomposition & solve, generalized symmetric eigenvalue problem added
 * 0.51 (01/01/15): infinite norm corrected
 * 0.52 (10/01/15): methods symmetrize, inner, svqb, gs added
 * 0.53 (10/07/15): maxD added
 * 0.54 (15/10/15): pinv bug due to svd change corrected.
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Templated classes for 3D and dense matrices.
*/

#define unaryOpScalar3(sym) \
	const mat3<T>& operator sym (const T b) \
		{ for (unsigned int i = 0; i < 16; i++) x[i] sym b; return *this; }

//! Implements a templated 4 x 4 matrix in 3D

template <typename T> class mat3 {
public:
	T	x[16];

	mat3<T>()									//! allocates a 4x4 matrix.
		{ for (unsigned int i = 0; i < 16; i++) x[i] = T(0); }
        T&	operator[](const unsigned int i)					//! returns reference to element i.
		{ return x[i]; }
        T	operator[](const unsigned int i) const					//! returns element i.
		{ return x[i]; }
	T&	operator()(const unsigned int m, const unsigned int n)			//! access reference to element in mat(m,n) notation.
		{ return x[n*4+m]; }
	T	operator()(const unsigned int m, const unsigned int n) const		//! access element in mat(m,n) notation.
		{ return x[n*4+m]; }
	mat3<T> operator*(const mat3<T>& b) const					//! matrix-matrix multiplication.
		{ mat3<T> r;
		  for (unsigned int i = 0; i < 4; i++)
		  	for (unsigned int j = 0; j < 4; j++)
	    			r.x[i*4+j] = x[i*4+0]*b.x[0*4+j]+x[i*4+1]*b.x[1*4+j]
					+x[i*4+2]*b.x[2*4+j]+x[i*4+3]*b.x[3*4+j];
		  return r; }
	vec3<T>	operator*(const vec3<T>& b) const					//! matrix-vector multiplication.
		{ T out[3]; for (unsigned int i = 0; i < 3; i++)
			out[i] = b.x*x[0+i]+b.y*x[4+i]+b.z*x[8+i]+x[12+i];
		  return vec3<T>(out[0],out[1],out[2]); }
	vec3<T>	ndc(const vec3<T>& b) const						//! matrix-vector projection into normalized device coordinates.
		{ T out[4]; for (unsigned int i = 0; i < 4; i++)
			out[i] = b.x*x[0+i]+b.y*x[4+i]+b.z*x[8+i]+x[12+i];
		  return vec3<T>(out[0],out[1],out[2])/out[3]; }
	mat3<T>& translate(const vec3<T>& v)						//! translates by vector v.
		{ x[12] = x[0]*v.x+x[4]*v.y+x[8]*v.z+x[12]; x[13] = x[1]*v.x+x[5]*v.y+x[9]*v.z+x[13];
		  x[14] = x[2]*v.x+x[6]*v.y+x[10]*v.z+x[14]; x[15] = x[3]*v.x+x[7]*v.y+x[11]*v.z+x[15];
		  return *this; }
	mat3<T>& scale(const vec3<T>& v)							//! scales by vector v.
		{ mat3 m; m(0,0) = v.x; m(1,1) = v.y; m(2,2) = v.z; m(3,3) = 1.0;
		  *this = m* *this; return *this; }
	mat3<T>& id()									//! sets matrix to identity.
		{ for (unsigned int i = 0; i < 16; i++) x[i] = 0;
		  for (unsigned int i = 0; i < 4; i++) x[i*5] = 1.0;
		  return *this; }
	mat3<T>& rotate(const T a, const vec3<T>& v)					//! rotates by angle a around axis v.
		{ mat3<T> A = *this, t;	
		  const T ca = std::cos(a), sa = std::sin(a);
		  const T ca1 = 1.0f-ca, aca = v.x*ca1, cca = v.z*ca1;
		  t[0] = v.x*aca+ca; t[1] = v.y*aca+v.z*sa; t[2] = v.z*aca-v.y*sa;
		  t[4] = v.y*aca-v.z*sa; t[5] = v.y*v.y*ca1+ca; t[6] = v.y*cca+v.x*sa;
		  t[8] = v.x*cca+v.y*sa; t[9] = v.y*cca-v.x*sa; t[10] = v.z*cca+ca;
		  x[0] = A[0]*t[0]+A[4]*t[1]+A[8]*t[2];    x[1] = A[1]*t[0]+A[5]*t[1]+A[9]*t[2];
		  x[2] = A[2]*t[0]+A[6]*t[1]+A[10]*t[2];   x[3] = A[3]*t[0]+A[7]*t[1]+A[11]*t[2];
		  x[4] = A[0]*t[4]+A[4]*t[5]+A[8]*t[6];    x[5] = A[1]*t[4]+A[5]*t[5]+A[9]*t[6];
		  x[6] = A[2]*t[4]+A[6]*t[5]+A[10]*t[6];   x[7] = A[3]*t[4]+A[7]*t[5]+A[11]*t[6];
		  x[8] = A[0]*t[8]+A[4]*t[9]+A[8]*t[10];   x[9] = A[1]*t[8]+A[5]*t[9]+A[9]*t[10];
		  x[10] = A[2]*t[8]+A[6]*t[9]+A[10]*t[10]; x[11] = A[3]*t[8]+A[7]*t[9]+A[11]*t[10];
		  return *this; }
	mat3<T>& rotate(const vec3<T>& r)						//! rotates coordinates axes by angles (r.x,r.y,r.z).
		{ rotate(r.x,vec3<T>(1,0,0)); rotate(r.y,vec3<T>(0,1,0));
		  rotate(r.z,vec3<T>(0,0,1)); return *this; }
	mat3<T>& frustum(const T left, const T right, const T bottom, const T top,
			const T znear, const T zfar)					//! sets frustum from specified values.
		{ const T t1 = T(2.0)*znear, t2 = right-left, t3 = top-bottom, t4 = zfar-znear;
		  x[0] = t1/t2; x[1] = 0; x[2] = 0; x[3] = 0;
		  x[4] = 0; x[5] = t1/t3; x[6] = 0; x[7] = 0; 
		  x[8] = (right+left)/t2; x[9] = (top+bottom)/t3; x[10] = (-zfar-znear)/t4; x[11] = -1.0; 
		  x[12] = 0; x[13] = 0; x[14] = (-t1*zfar)/t4; x[15] = 1.0;
		  return *this; }
	mat3<T>& perspective(const T fovy, const T aspectRatio, const T znear, const T zfar)	//! initializes perspective projection matrix.
		{ const T ymax = znear*T(std::tan(fovy*M_PI/360.0)), xmax = ymax*aspectRatio;
		  frustum(-xmax, xmax, -ymax, ymax, znear, zfar); return *this; }
	mat3<T>& lookAt(const fvec3& eye, const fvec3& ce, const fvec3& up)		//! initializes the model view matrix.
		{ const fvec3 f = (ce-eye).normalize(), s = cross(f,up).normalize(), u = cross(s,f);
		  x[0] = s.x; x[1] = u.x; x[2] = -f.x; x[3] = 0; 
		  x[4] = s.y; x[5] = u.y; x[6] = -f.y; x[7] = 0; 
		  x[8] = s.z; x[9] = u.z; x[10] = -f.z; x[11] = 0; 
		  x[12] = 0; x[13] = 0; x[14] = 0; x[15] = 1.0;
		  return translate(eye*-1.0f); }
	mat3<T>& ortho(const T left, const T right, const T bottom, const T top,
			const T znear, const T zfar)					//! initializes the orthographic projection matrix.
		{ const T t2 = right-left, t3 = top-bottom, t4 = zfar-znear;
		  x[0] = T(2.0)/t2; x[1] = 0; x[2] = 0; x[3] = 0;
		  x[4] = 0; x[5] = T(2.0)/t3; x[6] = 0; x[7] = 0; 
		  x[8] = 0; x[9] = 0; x[10] = T(-2.0)/t4; x[11] = 0;
		  x[12] = (-right-left)/t2; x[13] = (-top-bottom)/t3;
		  x[14] = (-zfar-znear)/t4; x[15] = 1.0; return *this; }
	mat3<T>& scaleTranslation(const T b)						//! applies scaling b to translation vector of the matrix.
		{ x[12] *= b; x[13] *= b; x[14] *= b; return *this; }
		unaryOpScalar3(=)
		unaryOpScalar3(+=)
		unaryOpScalar3(-=)
		unaryOpScalar3(*=)
		unaryOpScalar3(/=)
	vec3<T>	getRotation() const							//! returns rotation angles from matrix.
		{ const T cy = std::sqrt(SQR(x[0])+SQR(x[4])); vec3<T> r = T(0);
		  if (cy < T(1e-4)) {
			r.x = std::atan2(-x[6],x[5]), r.y = std::atan2(-x[8],T(0)); }
		  else { r.x = std::atan2(x[8]/cy,x[10]/cy);
			r.y = std::atan2(-x[8],cy);
			r.z = std::atan2(x[4]/cy,x[0]/cy); };
		  return r; }
};

//! \relates mat3
template <typename T>
void print(const mat3<T>& A)
//! prints matrix A.
{ 	for (unsigned int r = 0; r < 4; r++) {
		for (unsigned int c = 0; c < 4; c++) printT(A(r,c));
		printf("\n"); };
}

//! \relates mat3
template <typename T>
mat3<T> trp(const mat3<T>& A)
//! returns transpose of matrix A.
{	mat3<T> r;
	for (unsigned int i = 0; i < 4; i++)
		  for (unsigned int j = 0; j < 4; j++) r[j+i*4] = A[i+j*4];
	return r;
}

//! \relates mat3
template <typename T>
mat3<T> inv(const mat3<T>& A)
//! returns inverse of matrix A.
{	mat3<T> inv;
	inv[0]  =  A[5]*A[10]*A[15]-A[5]*A[11]*A[14]-A[9]*A[6]*A[15]+A[9]*A[7]*A[14]+A[13]*A[6]*A[11]-A[13]*A[7]*A[10];
	inv[4]  = -A[4]*A[10]*A[15]+A[4]*A[11]*A[14]+A[8]*A[6]*A[15]-A[8]*A[7]*A[14]-A[12]*A[6]*A[11]+A[12]*A[7]*A[10];
	inv[8]  =  A[4]*A[9]*A[15] -A[4]*A[11]*A[13]-A[8]*A[5]*A[15]+A[8]*A[7]*A[13]+A[12]*A[5]*A[11]-A[12]*A[7]*A[9];
	inv[12] = -A[4]*A[9]*A[14] +A[4]*A[10]*A[13]+A[8]*A[5]*A[14]-A[8]*A[6]*A[13]-A[12]*A[5]*A[10]+A[12]*A[6]*A[9];
	inv[1]  = -A[1]*A[10]*A[15]+A[1]*A[11]*A[14]+A[9]*A[2]*A[15]-A[9]*A[3]*A[14]-A[13]*A[2]*A[11]+A[13]*A[3]*A[10];
	inv[5]  =  A[0]*A[10]*A[15]-A[0]*A[11]*A[14]-A[8]*A[2]*A[15]+A[8]*A[3]*A[14]+A[12]*A[2]*A[11]-A[12]*A[3]*A[10];
	inv[9]  = -A[0]*A[9]*A[15] +A[0]*A[11]*A[13]+A[8]*A[1]*A[15]-A[8]*A[3]*A[13]-A[12]*A[1]*A[11]+A[12]*A[3]*A[9];
	inv[13] =  A[0]*A[9]*A[14] -A[0]*A[10]*A[13]-A[8]*A[1]*A[14]+A[8]*A[2]*A[13]+A[12]*A[1]*A[10]-A[12]*A[2]*A[9];
	inv[2]  =  A[1]*A[6]*A[15] -A[1]*A[7]*A[14] -A[5]*A[2]*A[15]+A[5]*A[3]*A[14]+A[13]*A[2]*A[7] -A[13]*A[3]*A[6];
	inv[6]  = -A[0]*A[6]*A[15] +A[0]*A[7]*A[14] +A[4]*A[2]*A[15]-A[4]*A[3]*A[14]-A[12]*A[2]*A[7] +A[12]*A[3]*A[6];
	inv[10] =  A[0]*A[5]*A[15] -A[0]*A[7]*A[13] -A[4]*A[1]*A[15]+A[4]*A[3]*A[13]+A[12]*A[1]*A[7] -A[12]*A[3]*A[5];
	inv[14] = -A[0]*A[5]*A[14] +A[0]*A[6]*A[13] +A[4]*A[1]*A[14]-A[4]*A[2]*A[13]-A[12]*A[1]*A[6] +A[12]*A[2]*A[5];
	inv[3]  = -A[1]*A[6]*A[11] +A[1]*A[7]*A[10] +A[5]*A[2]*A[11]-A[5]*A[3]*A[10]-A[9]*A[2]*A[7]  +A[9]*A[3]*A[6];
	inv[7]  =  A[0]*A[6]*A[11] -A[0]*A[7]*A[10] -A[4]*A[2]*A[11]+A[4]*A[3]*A[10]+A[8]*A[2]*A[7]  -A[8]*A[3]*A[6];
	inv[11] = -A[0]*A[5]*A[11] +A[0]*A[7]*A[9]  +A[4]*A[1]*A[11]-A[4]*A[3]*A[9] -A[8]*A[1]*A[7]  +A[8]*A[3]*A[5];
	inv[15] =  A[0]*A[5]*A[10] -A[0]*A[6]*A[9]  -A[4]*A[1]*A[10]+A[4]*A[2]*A[9] +A[8]*A[1]*A[6]  -A[8]*A[2]*A[5];
	T det = A[0]*inv[0] + A[1]*inv[4] + A[2]*inv[8] + A[3]*inv[12];
	if (det == T(0)) printf("Matrix is singular\n");
	det = 1.0f/det; for (unsigned i = 0; i < 16; i++) inv[i] *= det;
	return inv;
}

//! \relates mat3
template <typename T>
T det(const mat3<T>& A)
//! returns determinant of matrix A.
{	float a =  A[5]*A[10]*A[15]-A[5]*A[11]*A[14]-A[9]*A[6]*A[15]+A[9]*A[7]*A[14]+A[13]*A[6]*A[11]-A[13]*A[7]*A[10];
	float b = -A[4]*A[10]*A[15]+A[4]*A[11]*A[14]+A[8]*A[6]*A[15]-A[8]*A[7]*A[14]-A[12]*A[6]*A[11]+A[12]*A[7]*A[10];
	float c =  A[4]*A[9]*A[15] -A[4]*A[11]*A[13]-A[8]*A[5]*A[15]+A[8]*A[7]*A[13]+A[12]*A[5]*A[11]-A[12]*A[7]*A[9];
	float d = -A[4]*A[9]*A[14] +A[4]*A[10]*A[13]+A[8]*A[5]*A[14]-A[8]*A[6]*A[13]-A[12]*A[5]*A[10]+A[12]*A[6]*A[9];
	return A[0]*a+A[1]*b+A[2]*c+A[3]*d;
}

//! \relates mat3
template <typename T>
void toString(char* buf, const mat3<T>& A)
//! converts matrix A into string in *buf.
{	char str[32]; *buf = 0;
	for (unsigned int i = 0; i < 16; i++) {
		sprintT(str, A[i]); strcat(buf, str); };
}

//! \relates mat3
template <typename T>
void toMatrix(mat3<T>& A, const char* buf)
//! converts string in *buf into matrix A.
{	const char* p = buf; double v = 0; int n = 0;
	for (unsigned int i = 0; i < 16; i++) {
		sscanf(p, "%lf%n", &v, &n); A[i] = v; p += n; };
}

//! \cond INTERNAL
using fmat3 = mat3<float>;
using dmat3 = mat3<double>;
using cmat3 = mat3<fcomplex>;
using zmat3 = mat3<dcomplex>;
//! \endcond

#define allElem(i) (unsigned int i = 0; i < M*N; i++) 

#define unaryOpScalar(sym) \
	const matD<T>& operator sym (const T b) \
	{ for allElem(i) x[i] sym b; return *this; }
#define binaryOpScalar(sym) \
	matD<T> operator sym (const T b) const \
	{ matD<T> a(M,N); for allElem(i) a.x[i] = x[i] sym b; return a; }
#define unaryOpMatrix(sym) \
	const matD<T>& operator sym (const matD<T>& b) \
	{ assert(N == b.N && M == b.M); for allElem(i) x[i] sym b.x[i]; return *this; }
#define binaryOpMatrix(sym) \
	matD<T>	operator sym (const matD<T>& b) const \
	{ assert(N == b.N && M == b.M); matD<T> a(M,N); for allElem(i) a.x[i] = x[i] sym b.x[i]; return a; }

//! Implements a templated dense matrix

template<typename T> class matD {
	void	checksize(const unsigned int _M, const unsigned int _N) const 		//! checks that computations on indices fit into an unsigned int
		{ const size_t s = static_cast<size_t>(_M)*static_cast<size_t>(_N);
		  assert(s <= static_cast<size_t>(~0u)); }
	void	alloc(const unsigned int _M, const unsigned int _N)			//! allocates space for (M,N) elements.
		{ M = _M; N = _N; x = new T [size_t(_M)*size_t(_N)]; }
	void	clear()									//! frees element array.
		{ delete [] x; x = nullptr; M = 0; N = 0; }
	void	assign(const matD<T>& b)						//! assigns from matrix b.
		{ if (M != b.M || N != b.N) { clear(); alloc(b.M,b.N); };
		  for allElem(i) x[i] = b.x[i]; }
	void	assign(const vecD<T>& b)						//! assigns a (n,1) matrix from vector b.					// (re)allocate&  copy
		{ if (M != b.N || N != 1) { clear(); alloc(b.N,1); }
		  for allElem(i) x[i] = b(i); }
	T	mrc(const T* ax, const T* bx) const					//! multiplies row *ax by column *bx.
		{ T sum{0}; for (unsigned int i = 0; i < N; i++) sum += ax[i*M]*bx[i];
		  return sum; }
	T	mcc(const T* ax, const T* bx) const					//! multiplies column *ax by column *bx.
		{ T sum{0}; for (unsigned int i = 0; i < M; i++) sum += ax[i]*bx[i];
		  return sum; }
	void	mv(T *y, T *a, T *x) const						//! multiplies matrix *xi by vector *bx, returns *ax.
		{ T *d = y, *aj = a, *xj = x;						// optimized for vector operation.
		  for (unsigned int i = 0; i < M; i++) *d++ = 0.0;
		  for (unsigned int j = 0; j < N; j++) { d = y; 
			for (unsigned int i = 0; i < M; i++) { *d++ += (*aj++) * (*xj); }
			xj++; } }
	matD<T>	getmat(const unsigned int rs, const unsigned int re,
			const unsigned int cs, const unsigned int ce) const		//! access sub-matrix (rs:re,cs:ce).
		{ assert(rs <= M && re <= M && rs <= re); const unsigned int m = re-rs; 
		  assert(cs <= N && ce <= N && cs <= ce); const unsigned int n = ce-cs;
		  matD<T> A(m,n);
		  for (unsigned int i = 0; i < m; i++)
			for (unsigned int j = 0; j < n; j++)
				A(i,j) = (*this)(i+rs,j+cs);
		  return A; }
	vecD<T>	getrow(const unsigned int r, const unsigned int cs,
			const unsigned int ce) const					//! access row (r,cs:ce).
		{ assert(r <= M && cs <= N && ce <= N && cs <= ce);
		  const unsigned int n = ce-cs; vecD<T> V(n);
		  for (unsigned int i = 0; i < n; i++) V(i) = (*this)(r,i+cs);
		  return V; }
	vecD<T>	getcol(const unsigned int rs, const unsigned int re,
			const unsigned int c) const					//! access column (rs:re,c).
		{ assert(rs <= M && re <= M && rs <= re && c <= N);
		  const unsigned int n = re-rs; vecD<T> V(n);
		  for (unsigned int i = 0; i < n; i++) V(i) = (*this)(rs+i,c);
		  return V; }
public:
	unsigned int M;				//!< number of rows
	unsigned int N;				//!< number of columns
	T*	x;				//!< pointer to data array 

	matD<T>()									//! allocates an empty matrix.
		 : M(0), N(0), x(nullptr) { }
	matD<T>(const unsigned int _M, const unsigned int _N)				//! allocates a (M,N) matrix.
		 : M(_M), N(_N), x(nullptr) { alloc(M, N); }
	matD(const matD<T>& b)								//! copies from matrix b.
		 : M(0), N(0), x(nullptr) { assign(b); }
	matD(const vecD<T>& b)								//! copies from vector b.
		 : M(0), N(0), x(nullptr) { assign(b); }
	matD(matD<T>&& b)								//! moves from matrix b.
		 : M(b.M), N(b.N), x(b.x) { b.x = nullptr; }
	~matD<T>()	{ clear(); }
	matD<T>& 	operator=(const matD<T>& b)					//! assigns from matrix b.
			{ if (this != &b) assign(b); return *this; }
	matD<T>& 	operator=(const vecD<T>& b)					//! assigns from vector b.
			{ if (this != &b) assign(b); return *this; }
	matD<T>& 	operator=(matD<T>&& b)						//! move assigns from matrix b.
			{ assert(this != &b); delete [] x; x = b.x; b.x = nullptr;
			  M = b.M; N = b.N; return *this; }
	unsigned int	m() const							//! returns number of rows.
			{ return M; }
	unsigned int	n() const							//! returns number of columns.
			{ return N; }
	size_t		nel() const
			{ return static_cast<size_t>(M)*static_cast<size_t>(N); }
	T&		operator()(unsigned int m, unsigned int n)			//! returns reference to element (m,n).
			{ return x[n*M+m]; }
	T		operator()(unsigned int m, unsigned int n) const		//! returns element (m,n).
			{ return x[n*M+m]; }
			unaryOpScalar(=)
			unaryOpScalar(+=)
			unaryOpScalar(-=)
			unaryOpScalar(*=)
			unaryOpScalar(/=)
			binaryOpScalar(+)
			binaryOpScalar(-)
			binaryOpScalar(*)
			binaryOpScalar(/)
			unaryOpMatrix(+=)
			unaryOpMatrix(-=)
			binaryOpMatrix(+)
			binaryOpMatrix(-)
        const matD<T>	operator()(const underscore& , const range& c) const		//! returns sub-matrix (_,_(cs,ce)).
			{ return getmat(0,M,c.st,c.end); }
        const matD<T>	operator()(const range& r, const underscore& ) const		//! returns sub-matrix (_(rs,re),_).
			{ return getmat(r.st,r.end,0,N); }
        const matD<T>	operator()(const range& r, const range& c) const		//! returns sub-matrix (_(rs,re),_(cs,ce)).
			{ return getmat(r.st,r.end,c.st,c.end); }
        const vecD<T>	operator()(unsigned int r, const underscore& ) const		//! returns row (r,_).
			{ return getrow(r,0,N); }
	const vecD<T>	operator()(unsigned int r, const range& c) const		//! returns sub-row (r,_(cs,ce)).
			{ return getrow(r,c.st,c.end); }
	const vecD<T>	operator()(const underscore& , unsigned int c) const		//! returns column (_,r).
			{ return getcol(0,M,c); }
	const vecD<T>	operator()(const range& r, unsigned int c) const		//! returns sub-column (_(rs,re),c).
			{ return getcol(r.st,r.end,c); }
	matD<T>&	operator*=(const matD<T>& b)					//! unary matrix-matrix multiplication.
			{ assert(N == b.M); matD<T> a(M,b.N);
			  for (unsigned int n = 0; n < b.N; n++) mv(a.x+n*M, x, b.x+n*b.M);
			  *this = a; return *this; }
	matD<T>		operator*(const matD<T>& b) const 				//! binary matrix-matrix multiplication.
			{ assert(N == b.M); matD<T> a(M,b.N);
			  for (unsigned int n = 0; n < b.N; n++) mv(a.x+n*M, x, b.x+n*b.M);
			  return a; }
	vecD<T>		operator*(const vecD<T>& b) const				//! binary matrix-vector multiplication.
			{ assert(N == b.N); vecD<T> a(M); mv(a.x,x,b.x); return a; }
	matD<T>&	set(const matD<T>& A, unsigned int rs = 0, unsigned int cs = 0)	//! inserts sub-matrix A at (rs,cs).
			{ assert(rs+A.M <= M && cs+A.N <= N);
			  for (unsigned int i = 0; i < A.M; i++)
			  	for (unsigned int j = 0; j < A.N; j++)
					(*this)(i+rs,j+cs) = A(i,j);
			  return *this; }
	matD<T>&	set(const vecD<T>& V, unsigned int c)				//! inserts vector V at column c.
			{ assert(V.N <= M && c < N);
			  for (unsigned int i = 0; i < V.N; i++) (*this)(i,c) = V(i);
			  return *this; }
	matD<T>&	set(const unsigned int r, const vecD<T>& V)			//! inserts vector V at row r.
			{ assert(V.N <= N && r < M);
			  for (unsigned int i = 0; i < V.N; i++) (*this)(r,i) = V(i);
			  return *this; }
	matD<T>&	resize(const unsigned int _M, const unsigned int _N)		//! resize matrix for (M,N) elements.
			{ checksize(_M, _N); clear(); alloc(_M, _N); return *this; }
	T		sum() const							//! returns sum over all elements.
			{ T s{0}; for allElem(i) { s += x[i]; }; return s; }
        vecD<T>         tm(const vecD<T> &b) const                      		//! multiplies vecD<T> with matD<T> transpose
                        { assert(M == b.N); vecD<T> a(N); T* ax = a.x;
                          for (unsigned int i = 0; i < N; i++) *ax++ = mcc(b.x,x+i*M);
			  return a; };
	T		chisq(const matD<T>& b)						//! returns sum of squared differences between this matrix and b.
			{ assert(N == b.N && M == b.M);
			  T s{0}; for allElem(i) { s += SQR(x[i]-b.x[i]); }; return s; }
	T		trace() const							//! returns trace of this matrix.
			{ assert(N == M); T s{0};
			  for (unsigned int i = 0; i < N; i++) s += x[i*M+i];
			  return s; }
	vecD<T>		getCol(const unsigned int j) const				//! returns column j as vector.
			{ assert(j < N); vecD<T> a(M); T *ai = a.x, *xi = x+j*M;
			  for (unsigned int i = 0; i < M; i++) *ai++ = *xi++;
			  return a; }
	matD<T>&	setCol(const unsigned int j, const vecD<T>& b)			//! sets column j from vector b.
			{ assert(j < N && M == b.N); T *bi = b.x, *xi = x+j*M;
			  for (unsigned int i = 0; i < M; i++) *xi++ = *bi++;
			  return *this; }
	matD<T>&	addCol(const vecD<T>& b) 					//! join vector b as column to this matrix.
			{ if (M) assert(M == b.N); else { M = b.N; N = 0; };
			  matD<T> A(M,N+1);
			  for allElem(i) A.x[i] = x[i];
			  for (unsigned int i = 0; i < M; i++) A(i,N) = b(i);
			  *this = A; return *this; }
	matD<T>&	delCol(const unsigned int c)
			{ assert(c < N); T* bx = new T [M*(N-1)];
			  for (unsigned int i = 0; i < c*M; i++) bx[i] = x[i];
			  for (unsigned int i = c*M; i < M*(N-1); i++) bx[i] = x[i+M];
			  delete [] x; x = bx; N--; return *this; }
	matD<T>&	swapCol(unsigned int a, unsigned int b)				//! swap columns (a,b).
			{ assert(a < N && b < N); if (a == b) return *this;
			  for (unsigned int i = 0; i < M; i++) {
				T t = x[a*M+i]; x[a*M+i] = x[b*M+i]; x[b*M+i] = t; };
			  return *this; }
	vecD<T>		getRow(const unsigned int j) const				//! returns row j as vector.
			{ assert(j < M); vecD<T> a(N); T *ai = a.x, *xi = x+j;
			  for (unsigned int i = 0; i < N; i++) *ai++ = xi[i*M];
			  return a; }
	matD<T>&	setRow(const unsigned int j, const vecD<T>& b)			//! sets row j from vector b.
			{ assert(j < M && N == b.N); T *bi = b.x, *xi = x+j;
			  for (unsigned int i = 0; i < N; i++) xi[i*M] = *bi++;
			  return *this; }
	matD<T>&	addRow(const vecD<T>& b) 					//! join vector b as row to this matrix.
			{ assert(N == b.N); matD<T> A(M+1,N); unsigned int k = 0;
			  for (unsigned int j = 0; j < N; j++)
			  	for (unsigned int i = 0; i < M; i++) A(i,j) = x[k++];
			  for (unsigned int i = 0; i < N; i++) A(M,i) = b(i);
			  *this = A; return *this; }
	matD<T>&	delRow(const unsigned int r)
			{ assert(r < M); T* bx = new T [(M-1)*N]; unsigned int k = 0, l = 0;
			  for (unsigned int j = 0; j < N; j++) {
			  	for (unsigned int i = 0; i < M; i++, l++)
					if (i == r) k++; else bx[l-k] = x[l]; }
			  assert(k == N); delete [] x; x = bx; M--; return *this; }
	matD<T>&	swapRow(unsigned int a, unsigned int b)				//! swap rows (a,b).
			{ assert(a < M && b < M); if (a == b) return *this;
			  for (unsigned int i = 0; i < N; i++) {
				T t = x[i*M+a]; x[i*M+a] = x[i*M+b]; x[i*M+b] = t; };
			  return *this; }
	matD<T>&	left()								//! performs a circular shift left of this matrix.
			{ for (unsigned int i = 0; i < M; i++) { T *xi = x+i; T t = *xi;
				for (unsigned int j = 0; j < N-1; j++) xi[j*M] = xi[(j+1)*M];
				xi[(N-1)*M] = t; };
			  return *this; }
	matD<T>&	right()								//! performs a circular shift right of this matrix.
			{ for (unsigned int i = 0; i < M; i++) { T *xi = x+i; T t = xi[(N-1)*M];
				for (unsigned int j = 0; j < N-1; j++) xi[(j+1)*M] = xi[j*M];
				*xi = t; };
			  return *this; }
	matD<T>&	id()								//! sets matrix to identity.
			{ *this = T(0); addToD(T(1)); return *this; }
	matD<T>&	addToD(const vecD<T>& b)					//! adds vector b to diagonal.
			{ const unsigned int K = std::min(M,N);
			  for (unsigned int i = 0; i < K; i++) x[i*M+i] += b(i);
			  return *this; }
	matD<T>&	addToD(const T b)						//! add scalar b to diagonal.
			{ const unsigned int K = std::min(M,N);
			  for (unsigned int i = 0; i < K; i++) x[i*M+i] += b;
			  return *this; }
	T		maxD() const							//! returns maximum element on diagonal.
			{ const unsigned int K = std::min(M,N);
			  T vm = std::numeric_limits<T>::lowest();
			  for (unsigned int i = 0; i < K; i++) {
				const T v = x[i*M+i]; vm = std::max(vm,v); };
			  return vm; }
	unsigned int	imax() const							//! returns linear index of maximum element.
			{ unsigned int imax = 0; T vmax = std::numeric_limits<T>::lowest();
			  for (unsigned int i = 0; i < M*N; i++) {
				if (x[i] > vmax) { vmax = x[i]; imax = i; } };
			  return imax; }
	unsigned int	amax() const							//! returns linear index of absolute maximum element.
			{ unsigned int imax = 0; T vmax{0};
			  for (unsigned int i = 0; i < M*N; i++) {
				if (std::abs(x[i]) > vmax) { vmax = std::abs(x[i]); imax = i; } };
			  return imax; }
	matD<T>&	scaleEigenvectors(const vecD<T>& wi)				//! scales complex conjugate eigenvectors to unit length.
	  		{ for (unsigned int i = 0; i < N; i++) {
				vecD<T> e = getCol(i); T en = norm(e);
				if (wi(i) == 0) { e /= en; setCol(i,e); }		// scale real eigenvector
				else if (wi(i) > 0) {					// scale conjugate pair
					vecD<T> ei = getCol(i+1); T ein = norm(ei);
					T sc = T(1.0/std::sqrt(SQR(en)+SQR(ein)));
					e *= sc; ei *= sc; setCol(i,e); setCol(i+1,ei); } };
			  return *this; }
	matD<T>&	scaleEigenvectors()						//! scales real eigenvectors to unit length.
	  		{ for (unsigned int i = 0; i < N; i++) {
				vecD<T> e = getCol(i); T en = norm(e); e /= en; setCol(i,e); };
			  return *this; }
	matD<T>&	splitConjugatePairs(matD<std::complex<T>>& Q, vecD<T>& w)		//! modifies eigenvectors of complex conjugate pairs.
			{ Q.resize(N,N); for (unsigned int i = 0; i < N; i++) {
				if (w[i] == T(0)) {
					for (unsigned int j = 0; j < N; j++)
						Q(j,i) = std::complex<T>(x[i*N+j],T(0)); }
	  			else if (w[i] > T(0)) {
					for (unsigned int j = 0; j < N; j++)
						Q(j,i) = std::complex<T>(x[i*N+j],x[(i+1)*N+j]); }
	  			else {	for (unsigned int j = 0; j < N; j++)
						Q(j,i) = std::complex<T>(x[(i-1)*N+j],-x[i*N+j]); } };
			 return *this; }
	vecD<T>		d() const							//! returns diagonal.
			{ unsigned int n = std::min(N,M); vecD<T> a(n); 
			  for (unsigned int i = 0; i < n; i++) a[i] = x[i*M+i];
			  return a; }
	unsigned int	nnz(const T lim = T(0)) const					//! returns number of non-zero elements.
			{ unsigned int n = 0; for allElem(i) n += std::abs(x[i]) > lim;
			  return n; }
	bool		finite() const							//! checks that all elements are finite.
			{ for allElem(i) if (std::isfinite(x[i]) == false) return false;
			  return true; }
	matD<T>&	joinBelow(const matD<T>& A);
	matD<T>&	joinRight(const matD<T>& A);
	vecD<T>		rowsum() const							//! returns a vector of row-wise sums.
			{ vecD<T> a(M); a = T(0); for (unsigned int i = 0; i < M; i++) {
				for (unsigned int j = 0; j < N; j++) a[i] += x[j*M+i]; };
			  return a; }
	vecD<T>		colsum() const							//! returns a vector of column-wise sums.
			{ vecD<T> a(N); a = T(0); for (unsigned int j = 0; j < N; j++) {
				for (unsigned int i = 0; i < M; i++) a[j] += x[j*M+i]; };
			  return a; }
	vecD<T>		rowdeg() const							//! returns a vector of off-diagonal row-wise sums.
			{ vecD<T> a(M); a = T(0); for (unsigned int i = 0; i < M; i++) {
				for (unsigned int j = 0; j < N; j++)
					if (i != j) a[i] += x[j*M+i]; };
			  return a; }
	vecD<T>		coldeg() const							//! returns a vector of off-diagonal column-wise sums.
			{ vecD<T> a(N); a = T(0); for (unsigned int j = 0; j < N; j++) {
				for (unsigned int i = 0; i < M; i++)
					if (i != j) a[j] += x[j*M+i]; };
			  return a; }
	uvecD		rowmax() const							//! returns vector of row-wise maxima.
		  	{ uvecD a(M); for (unsigned int i = 0; i < M; i++) { 
				T m = std::numeric_limits<T>::lowest(); unsigned int ix = ~0u;
			  	for (unsigned int j = 0; j < N; j++) { const T v = x[j*M+i]; 
				 	if (v > m) { m = v; ix = j; } };
				a[i] = ix; };
			  return a; }
	uvecD		colmax() const							//! returns vector of column-wise maxima.
		  	{ uvecD a(N); for (unsigned int j = 0; j < N; j++) {
				T m = std::numeric_limits<T>::lowest(); unsigned int ix = ~0u;
				for (unsigned int i = 0; i < M; i++) { const T v = x[j*M+i]; 
				 	if (v > m) { m = v; ix = i; } };
				a[j] = ix; };
			  return a; }
	T		nmi() const							//! computes normalized mutual information.
			{ const T t = sum(); if (t == T(0)) throw rtException("nmi: zero matrix");
			  vecD<T> hi(M), hj(N); hi = T(0); hj = T(0); T e = T(0);	// normalize joint histogram
			  for (unsigned int i = 0; i < M; i++)				// marginalize histograms
				for (unsigned int j = 0; j < N; j++) {
					T h = x[j*M+i]/t; if (h == T(0)) continue;
					hi[i] += h; hj[j] += h; };
			  for (unsigned int i = 0; i < M; i++)
				for (unsigned int j = 0; j < N; j++) {
					T h = x[j*M+i]/t, d = hi[i]*hj[j];
					if (h == T(0) || d == T(0)) continue;
					e -= h*std::log(h/d); };
			  const T em = T(0.5)*(hi.entropy()+hj.entropy()); 		// compute entropies
			  return isSmall<T>(em)? T(0): e/em; }
};

#undef unaryOpScalar
#undef binaryOpScalar
#undef unaryOpMatrix
#undef binaryOpMatrix
#undef allElem

//! \cond INTERNAL
using fmatD = matD<float>;
using dmatD = matD<double>;
using cmatD = matD<fcomplex>;
using zmatD = matD<dcomplex>;
using vfmat3 = std::vector<fmat3>;
using vdmat3 = std::vector<dmat3>;
using vfmatD = std::vector<fmatD>;
using vdmatD = std::vector<dmatD>;
using imatD = matD<int>;
//! \endcond

// overloaded lapack interface functions
extern void getri(fmatD& A, ivecD& ipiv);
extern void getri(dmatD& A, ivecD& ipiv);
extern void getri(cmatD& A, ivecD& ipiv);
extern void getri(zmatD& A, ivecD& ipiv);
extern void getrf(fmatD& A, ivecD& ipiv);
extern void getrf(dmatD& A, ivecD& ipiv);
extern void getrf(cmatD& A, ivecD& ipiv);
extern void getrf(zmatD& A, ivecD& ipiv);
extern float gecon(fmatD& A, float nm);
extern double gecon(dmatD& A, double nm);
extern float gecon(cmatD& A, float nm);
extern double gecon(zmatD& A, double nm);
extern void gesv(fmatD& A, fvecD& b);
extern void gesv(dmatD& A, dvecD& b);
extern void gesv(cmatD& A, cvecD& b);
extern void gesv(zmatD& A, zvecD& b);
extern void gesv(fmatD& A, fmatD& B);
extern void gesv(dmatD& A, dmatD& B);
extern void gesv(cmatD& A, cmatD& B);
extern void gesv(zmatD& A, zmatD& B);
extern void gels(fmatD& A, fvecD& b);
extern void gels(dmatD& A, dvecD& b);
extern void gels(cmatD& A, cvecD& b);
extern void gels(zmatD& A, zvecD& b);
extern void gesvd(fmatD& A, fvecD& d, fmatD& U, fmatD& V, const bool econ);
extern void gesvd(dmatD& A, dvecD& d, dmatD& U, dmatD& V, const bool econ);
extern void gesvd(cmatD& A, fvecD& d, cmatD& U, cmatD& V, const bool econ);
extern void gesvd(zmatD& A, dvecD& d, zmatD& U, zmatD& V, const bool econ);
extern void gesvd(fmatD& A, fvecD& d);
extern void gesvd(dmatD& A, dvecD& d);
extern void gesvd(cmatD& A, fvecD& d);
extern void gesvd(zmatD& A, dvecD& d);
extern void syev(fmatD& A, fvecD& d);
extern void syev(dmatD& A, dvecD& d);
extern void heevd(cmatD& A, fvecD& d);
extern void heevd(zmatD& A, dvecD& d);
extern bool potrf(fmatD& A);
extern bool potrf(dmatD& A);
extern bool potrf(cmatD& A);
extern bool potrf(zmatD& A);
extern bool trtri(fmatD& A);
extern bool trtri(dmatD& A);
extern bool trtri(cmatD& A);
extern bool trtri(zmatD& A);
extern bool potrs(fmatD& A, fvecD& b);
extern bool potrs(dmatD& A, dvecD& b);
extern bool potrs(cmatD& A, cvecD& b);
extern bool potrs(zmatD& A, zvecD& b);
extern void geev(fmatD& A, fvecD& wr, fvecD& wi, fmatD& V);
extern void geev(dmatD& A, dvecD& wr, dvecD& wi, dmatD& V);
extern void geev(cmatD& A, cvecD& w, cmatD& V);
extern void geev(zmatD& A, zvecD& w, zmatD& V);
extern void ggev(fmatD& A, fmatD& B, fvecD& ar, fvecD& ai, fvecD& bt, fmatD& V);
extern void ggev(dmatD& A, dmatD& B, dvecD& ar, dvecD& ai, dvecD& bt, dmatD& V);
extern void ggev(cmatD& A, cmatD& B, cvecD& al, cvecD& bt, cmatD& V);
extern void ggev(zmatD& A, zmatD& B, zvecD& al, zvecD& bt, zmatD& V);
extern void geqrf(fmatD& A, fvecD& d);
extern void geqrf(dmatD& A, dvecD& d);
extern void geqrf(cmatD& A, cvecD& d);
extern void geqrf(zmatD& A, zvecD& d);
extern void orgqr(fmatD& A, fvecD& d);
extern void orgqr(dmatD& A, dvecD& d);
extern void orgqr(cmatD& A, cvecD& d);
extern void orgqr(zmatD& A, zvecD& d);
extern void gees(fmatD& A, fvecD& wr, fvecD& wi, fmatD& V);
extern void gees(dmatD& A, dvecD& wr, dvecD& wi, dmatD& V);
extern void gees(cmatD& A, cvecD& w, cmatD& V);
extern void gees(zmatD& A, zvecD& w, zmatD& V);
extern void gges(fmatD& A, fmatD& B, fmatD& Vl, fmatD& Vr);
extern void gges(dmatD& A, dmatD& B, dmatD& Vl, dmatD& Vr);
extern void gges(cmatD& A, cmatD& B, cmatD& Vl, cmatD& Vr);
extern void gges(zmatD& A, zmatD& B, zmatD& Vl, zmatD& Vr);
extern void gehrd(fmatD& A, fvecD& d);
extern void gehrd(dmatD& A, dvecD& d);
extern void gehrd(cmatD& A, cvecD& d);
extern void gehrd(zmatD& A, zvecD& d);
extern void hseqr(fmatD& A, fvecD& wr, fvecD& wi, fmatD& V);
extern void hseqr(dmatD& A, dvecD& wr, dvecD& wi, dmatD& V);
extern void hseqr(cmatD& A, cvecD& w, cmatD& V);
extern void hseqr(zmatD& A, zvecD& w, zmatD& V);
extern void trevc(fmatD& A, fmatD& V);
extern void trevc(dmatD& A, dmatD& V);
extern void trevc(cmatD& A, cmatD& V);
extern void trevc(zmatD& A, zmatD& V);
extern void orghr(fmatD& A, fvecD& d);
extern void orghr(dmatD& A, dvecD& d);
extern void orghr(cmatD& A, cvecD& d);
extern void orghr(zmatD& A, zvecD& d);
extern void trexc(fmatD& A, fmatD& V, int j1, int j2);
extern void trexc(dmatD& A, dmatD& V, int j1, int j2);
extern void trexc(cmatD& A, cmatD& V, int j1, int j2);
extern void trexc(zmatD& A, zmatD& V, int j1, int j2);
extern void tgexc(fmatD& A, fmatD& B, fmatD& Q, fmatD& Z, int si, int di);
extern void tgexc(dmatD& A, dmatD& B, dmatD& Q, dmatD& Z, int si, int di);
extern void tgexc(cmatD& A, cmatD& B, cmatD& Q, cmatD& Z, int si, int di);
extern void tgexc(zmatD& A, zmatD& B, zmatD& Q, zmatD& Z, int si, int di);
extern void gebal(fmatD& A, int& ilo, int& ihi, uvecD& pm, fvecD& sc);
extern void gebal(dmatD& A, int& ilo, int& ihi, uvecD& pm, dvecD& sc);
extern void gebal(cmatD& A, int& ilo, int& ihi, uvecD& pm, fvecD& sc);
extern void gebal(zmatD& A, int& ilo, int& ihi, uvecD& pm, dvecD& sc);
extern void sygv(fmatD& A, fmatD& B, fvecD& v);
extern void sygv(dmatD& A, dmatD& B, dvecD& v);

// non-member dense matrix functions

#define allElem(i) (unsigned int i = 0; i < A.N*A.M; i++)
#define allCR(A) (unsigned int c = 0; c < A.N; c++) for (unsigned int r = 0; r < A.M; r++)

//! \relates matD
template <typename T>
matD<T> operator*(const T a, const matD<T>& A)
//! commutative matrix-scalar multiplication.
{ 	return A*a;
}

//! \relates matD
template <typename T>
void print(const matD<T>& A)
//! prints real matrix A to stdout.
{ 	for (unsigned int r = 0; r < A.M; r++) {
		for (unsigned int c = 0; c < A.N; c++) printT(A(r,c));
		printf("\n"); };
}

//! \relates matD
template<typename T, typename U> 
vecD<U> operator*(const matD<T>& A, const vecD<U>& b)
//! multiplies a real matrix A with a complex vector b.
{	assert(A.N == b.N); vecD<U> a(A.M);
	for (unsigned int m = 0; m < A.M; m++) {
		U s = U(0); for (unsigned int n = 0; n < A.N; n++) s += A(m,n)*b(n); a[m] = s; };
	return a;
}

//! \relates matD
template <typename T>
vecD<T> operator*(const vecD<T> a, const matD<T>& A)
//! vector-matrix multiplication.
{	assert(a.N == A.M); vecD<T> b(A.N); T* x = A.x;
	for (unsigned int c = 0; c < A.N; c++) { T s{0}, *ax = a.x;
		for (unsigned int r = 0; r < A.M; r++) s += (*x++)*(*ax++);		// dot product of a and column c of A.
		b[c] = s; };
	return b;
}

//! \relates matD
template <typename T> 
matD<std::complex<T>> asComplex(const matD<T>& A)
//! converts a real matrix A into a complex matrix.
{	matD<std::complex<T>> B(A.M,A.N);
	for allElem(i) B.x[i] = std::complex<T>(A.x[i],0.0);
	return B;
}

//! \relates matD
template <typename T> 
matD<T> real(const matD<std::complex<T>>& A)
//! returns real portion of a complex matrix A.
{	matD<T> B(A.M,A.N); for allElem(i) B.x[i] = A.x[i].real(); return B;
}

//! \relates matD
template <typename T> 
matD<T> imag(const matD<std::complex<T>>& A)
//! returns imaginary portion of a complex matrix A.
{	matD<T> B(A.M,A.N); for allElem(i) B.x[i] = A.x[i].imag(); return B;
}

//! \relates matD
template <typename T> 
matD<T> asMatrix(const vecD<T>& b)
//! converts a vector b into a matrix.
{  	matD<T> A(b.N,1);
	for (unsigned int i = 0; i < b.N; i++) A(i,0) = b(i);
	return A;
}

//! \relates matD
template <typename T> 
void random(matD<T>& A)
//! generates real random matrix A in [-1,1[.
{	uniformDist<T> ud(T(-1),T(1)); for allElem(i) A.x[i] = ud();
}

//! \relates matD
template <typename T> 
void random(matD<std::complex<T>>& A)
//! generates complex random matrix A in [-1,1[.
{	uniformDist<T> ud(T(-1),T(1));
	for allElem(i) A.x[i] = std::complex<T>(ud(),ud());
}

//! \relates matD
template <typename T> 
matD<T> conj(const matD<T>& A)
//! generates complex conjugate matrix of a real matrix A (noop copy).
{	matD<T> B = A; return B;
}

//! \relates matD
template <typename T> 
matD<std::complex<T>> conj(const matD<std::complex<T>>& A)
//! generates complex conjugate matrix of a complex matrix A.
{ 	matD<std::complex<T>> B(A.M,A.N);
	for allElem(i) B.x[i] = conj(A.x[i]);
	return B;
}

//! \relates matD
template <typename T> 
matD<T> adj(const matD<T>& A)
//! generates adjoint matrix of a real matrix A.
{	return trp(A);
}
			
//! \relates matD
template <typename T> 
matD<std::complex<T>> adj(const matD<std::complex<T>>& A)
//! generates adjoint matrix of a complex matrix A.
{ 	matD<std::complex<T>> B(A.N,A.M); std::complex<T> *d = B.x, *s = A.x; 
	for allCR(A) d[r*A.N+c] = conj(*s++); return B;
}

//! \relates matD
template <typename T>
matD<T> inv(const matD<T>& A)
//! returns inverse of square matrix A.
{	assert(A.M == A.N); matD<T> S = A; ivecD ipiv(A.N);
	getrf(S, ipiv);	getri(S, ipiv); return S;
}

//! \relates matD
template <typename T>
T cond(const matD<T>& A)
//! computes the condition number of real matrix A.
{	matD<T> S = A; const unsigned int ns = std::min(A.M,A.N);
	vecD<T> d(ns); gesvd(S, d); return d(0)/d(ns-1);
}

//! \relates matD
template <typename T>
std::complex<T> cond(const matD<std::complex<T>>& A)
//! computes the condition number of complex matrix A.
{	matD<std::complex<T>> S = A; const unsigned int ns = std::min(A.M,A.N);
	vecD<T> d(ns); gesvd(S, d); return d(0)/d(ns-1);
}

//! \relates matD
template <typename T>
T rcond(const matD<T>& A)
//! returns inverse of the condition number of real matrix A.
{	matD<T> S = A; ivecD ipiv(std::min(A.M,A.N));
	getrf(S, ipiv);	return gecon(S, norm1(A));
}

//! \relates matD
template <typename T>
T rcond(const matD<std::complex<T>>& A)
//! returns inverse of the condition number of complex matrix A.
{	matD<std::complex<T>> S = A; ivecD ipiv(std::min(A.M,A.N));
	getrf(S, ipiv);	return gecon(S, norm1(A));
}

//! \relates matD
template <typename T>
vecD<T> solve(const matD<T>& A, const vecD<T>& b)
//! solves A x = b for x, given vector b.
{	assert(A.M == A.N && A.M == b.N); if (A.M == 0) return vecD<T>();
	matD<T> S = A; vecD<T> x = b;
	gesv(S, x); return x;
}

//! \relates matD
template <typename T>
matD<T> solve(const matD<T>& A, const matD<T>& B)
//! solves A x = B for x, given matrix B.
{	assert(A.M == A.N && A.M == B.N); if (A.M == 0) return matD<T>();
	matD<T> S = A; matD<T> X = B; gesv(S, X); return X;
}

//! \relates matD
template <typename T>
vecD<T> lsqr(const matD<T>& A, const vecD<T>& b)
//! solves least-squares problem A x = b for x, given non-square matrix A.
{	assert(A.M == b.N); matD<T> S = A; vecD<T> x = b;
	gels(S, x); x.N = A.N; return x;
}

template <typename T> 
struct retMVM {
	matD<T>	U;				//!< left basis
	vecD<T> d;				//!< eigenvalues
	matD<T>	V;				//!< right basis
};

//! \relates matD
template <typename T> 
retMVM<T> svd(const matD<T>& A, const bool econ = false)
//! computes a singular value decomposition of a real matrix A;
//! returns left and right matrices (U,V), and diagonal vector.
{	const unsigned int m = A.M, n = A.N, ns = std::min(m, n);
	matD<T> S = A, U, V; vecD<T> d(ns);
	if (econ) { U.resize(m,ns); V.resize(ns,n); }
	else { U.resize(m,m); V.resize(n,n); }
	gesvd(S, d, U, V, econ); return { U,d,trp(V) };
}

template <typename T> 
struct retMcVMc {
	matD<std::complex<T>> U;		//!< left basis
	vecD<T> d;				//!< eigenvalues
	matD<std::complex<T>> V;		//!< right basis
};

//! \relates matD
template <typename T>
retMcVMc<T> svd(const matD<std::complex<T>>& A, const bool econ = false)
//! computes a singular value decomposition of a complex matrix A;
//! returns left and right matrices (U,V), and diagonal vector.
{  	const unsigned int m = A.M, n = A.N, ns = std::min(m, n);
	matD<std::complex<T>> S = A, U, V; vecD<T> d(ns);
	if (econ) { U.resize(m,ns); V.resize(ns, n); }
	else { U.resize(m, m); V.resize(n, n); }
	gesvd(S, d, U, V, econ); return { U,d,adj(V) };
}

//! \relates matD
template <typename T>
matD<T> pinv(const matD<T>& A)								
//! returns pseudo-inverse of rectangular real matrix A.
{	const auto S = svd(A); matD<T> W(A.N,A.M); W = T(0);
	const unsigned int cmax = std::min(S.V.N,W.N);
	for (unsigned int c = 0; c < cmax; c++) {
		T d = S.d(c) != T(0)? T(1)/S.d(c): T(0);
		for (unsigned int r = 0; r < S.V.M; r++) W(r,c) = S.V(r,c)*d; };
	return W*trp(S.U);
}

//! \relates matD
template <typename T>
matD<std::complex<T>> pinv(const matD<std::complex<T>>& A)								
//! returns pseudo-inverse of rectangular complex matrix A.
{  	const auto S = svd(A); matD<std::complex<T>> W(A.N,A.M); W = std::complex<T>(0);
	const unsigned int cmax = std::min(S.V.N,W.N);
	for (unsigned int c = 0; c < cmax; c++) {
		T d = S.d(c) != T(0)? T(1)/S.d(c): T(0);
		for (unsigned int r = 0; r < S.V.M; r++) W(r,c) = S.V(r,c)*d; };
	return W*adj(S.U);
}

template <typename T> 
struct retMV {
	matD<T>	U;				//!< basis
	vecD<T> d;				//!< eigenvalues
};

//! \relates matD
template <typename T>
retMV<T> sev(const matD<T>& A)
//! computes an eigen-decomposition of real symmetric matrix A;
//! returns eigenvectors in U and eigenvalues in d.
{	assert(A.M == A.N); matD<T> U = A; vecD<T> d(A.N);
	syev(U,d); return { U,d };
}

template <typename T> 
struct retMcV {
	matD<std::complex<T>> U;		//!< basis
	vecD<T> d;				//!< eigenvalues
};

//! \relates matD
template <typename T>
retMcV<T> sev(const matD<std::complex<T>>& A)
//! computes an eigen-decomposition of complex symmetric matrix A;
//! returns eigenvectors in U and eigenvalues in d.
{	assert(A.M == A.N); matD<T> U = A; vecD<T> d(A.N); 
	heevd(U,d); return { U,d };
}

//! \relates matD
template <typename T>
T norm1(const matD<T>& A)
//! computes the 1-norm of matrix A (max abs column sum).
{	T m{0}; for (unsigned int j = 0; j < A.N; j++) { T s{0};
		for (unsigned int i = 0; i < A.M; i++) s += std::abs(A(i,j));
		m = std::max(m,s); };
	return m;
}

//! \relates matD
template <typename T>
T norm2(const matD<T>& A)
//! computes the 2-norm of real matrix A.
{	const unsigned int ns = std::min(A.M,A.N);
	vecD<T> d(ns); matD<T> S = A; 
	gesvd(S, d); return d(0);
}

//! \relates matD
template <typename T>
T norm2(const matD<std::complex<T>>& A)
//! computes the 2-norm of complex matrix A.
{	const unsigned int ns = std::min(A.M,A.N);
	vecD<T> d(ns); matD<std::complex<T>> S = A; 
	gesvd(S, d); return d(0);
}

//! \relates matD
template <typename T>
T normI(const matD<T>& A)
//! computes the infinite norm of matrix A (max abs row sum).
{	T m{0}; for (unsigned int i = 0; i < A.M; i++) { T s{0};
		for (unsigned int j = 0; j < A.N; j++) s += std::abs(A(i,j));
		m = std::max(m,s); };
	return m;
}

//! \relates matD
template <typename T>
T normF(const matD<T>& A)
//! computes the Frobenius norm of matrix A.
{	T m{0}; for (unsigned int i = 0; i < A.M; i++)
		for (unsigned int j = 0; j < A.N; j++) m += SQR(A(i,j));
	return std::sqrt(m);
}

//! \relates matD
template <typename T>
bool isPSD(const matD<T>& A)
//! checks if matrix A is positive semi-definite.
{ 	assert(A.M == A.N); matD<T> S = A; return potrf(S);
}

//! \relates matD
template <typename T>
retMV<T> gsev(const matD<T>& A, const matD<T>& B)
//! computes a general eigen-decomposition of real symmetric matrix pair (A,B);
//! returns eigenvectors in Q and eigenvalues.
{ 	assert(A.M == A.N); matD<T> SA = A, SB = B; vecD<T> s(A.N);
	sygv(SA, SB, s); return { SA,s };
}

template <typename T> 
struct retMcVc {
	matD<std::complex<T>> U;		//!< basis
	vecD<std::complex<T>> d;		//!< eigenvalues
};

//! \relates matD
template <typename T>
retMcVc<T> nev(const matD<T>& A)
//! computes an eigen-decomposition of real non-symmetric matrix A;
//! returns eigenvectors in U and eigenvalues in d.
{ 	assert(A.M == A.N); const unsigned int n = A.N;
	matD<T> S = A, V(n,n); vecD<T> wr(n), wi(n); geev(S, wr, wi, V); 
	matD<std::complex<T>> C; V.splitConjugatePairs(C, wi);
	return { C,asComplex(wr,wi) };
}

//! \relates matD
template <typename T>
retMcVc<T> nev(const matD<std::complex<T>>& A)
//! computes an eigen-decomposition of complex non-symmetric matrix A;
//! returns eigenvectors in U and eigenvalues in d.
{	assert(A.M == A.N); const unsigned int n = A.N;
	matD<std::complex<T>> S = A, U(n,n); vecD<std::complex<T>> w(n); 
	geev(S, w, U); return { U,w };
}

template <typename T> 
struct retMcVcV {
	matD<std::complex<T>> U;		//!< basis
	vecD<std::complex<T>> alpha;		//!< left eigenvalues
	vecD<T> beta;				//!< right eigenvalues
};

//! \relates matD
template <typename T>
retMcVcV<T> gev(const matD<T>& A, const matD<T>& B)
//! computes a general eigen-decomposition of real non-symmetric matrix pair (A,B);
//! returns eigenvectors in U, left and right eigenvalues in (alpha,beta).
{ 	assert(A.M == A.N); const unsigned int n = A.N;
	matD<T> SA = A, SB = B, V(n,n), Vr; vecD<T> ar(n), ai(n), beta(n);
	ggev(SA, SB, ar, ai, beta, V); V.splitConjugatePairs(Vr, ai);
	return { Vr, asComplex(ar,ai), beta };
}

template <typename T> 
struct retMcVcVc {
	matD<std::complex<T>> U;		//!< basis
	vecD<std::complex<T>> alpha;		//!< left eigenvalues
	vecD<std::complex<T>> beta;		//!< right eigenvalues
};

//! \relates matD
template <typename T>
retMcVcVc<T> gev(const matD<std::complex<T>>& A, const matD<std::complex<T>>& B)
//! computes a general eigen-decomposition of complex non-symmetric matrix pair (A,B);
//! returns eigenvectors in Vr, left and right eigenvalues as vectors (alpha,beta).
{	assert(A.M == A.N); const unsigned int n = A.N;
	matD<std::complex<T>> SA = A, SB = B, U(n,n);
	vecD<std::complex<T>> alpha(n), beta(n);
	ggev(SA, SB, alpha, beta, U);
	return { U, alpha, beta };
}

//! \relates matD
template <typename T>
T det(const matD<T>& A)
//! returns determinant of a square matrix A.
{	assert(A.M == A.N); matD<T> S = A; ivecD ipiv(A.N);
	getrf(S, ipiv); T det = 1.0; bool neg = false; 
	for (unsigned int i = 0; i < A.N; i++) { 
		det *= S.x[i*A.N+i]; if (ipiv[i] != int(i+1)) neg = !neg; };
	return neg? -det: det;
}

//! \relates matD
template <typename T>
matD<T> log(const matD<T>& A)
//! returns logarithm of matrix A.
{	assert(A.M == A.N); matD<T> Ul(A.N,A.N); auto S = sev(A);
	for (unsigned int i = 0; i < S.d.N; i++) {
		if (S.d(i) > std::numeric_limits<T>::epsilon()) continue;
		S.d(i) = std::numeric_limits<T>::epsilon(); };
	for allCR(S.U) Ul(c,r) = std::log(S.d(c))*S.U(r,c); return S.U*Ul;
}

//! \relates matD
template <typename T>
matD<T> exp(const matD<T>& A)
//! returns exponential of matrix A.
{	assert(A.M == A.N); matD<T> Ul(A.N,A.N); const auto S = sev(A);
	for allCR(S.U) Ul(c,r) = T(std::exp(S.d(c))*S.U(r,c)); return S.U*Ul;
}

//! \relates matD
template <typename T>
matD<T> trp(const matD<T>& A)
//! returns transpose of matrix A.
{ 	matD<T> B(A.N,A.M); T *d = B.x, *s = A.x;
	for allCR(A) d[r*A.N+c] = *s++;
	return B;
}

template <typename T> 
struct retMM {
	matD<T> Q;			//!< basis
	matD<T> R;			//!< basis
};

//! \relates matD
template <typename T>
retMM<T> qr(const matD<T>& A)
//! returns a QR decomposition of matrix A in R, returns Q and R.
{	const unsigned int n = A.N, d = std::min(A.M,A.N);
	matD<T> Q = A, R(d,n); R = T(0); vecD<T> tau(d);
	geqrf(Q, tau); 
	for (unsigned int i = 0; i < R.M; i++)
		for (unsigned int j = 0; j <= i; j++) R(j,i) = Q(j,i);	  
	orgqr(Q, tau); return { Q, R };
}

//! \relates matD
template <typename T>
retMM<T> Schur(const matD<T>& A)
//! returns a Schur decomposition ZSZ of a real matrix A in Z, returns Z and S.
{	assert(A.M == A.N); const unsigned int n = A.N;
	matD<T> S = A, Z(n,n); vecD<T> wr(n), wi(n);
	gees(S, wr, wi, Z); return { Z, S };
}

template <typename T> 
struct retMcMc {
	matD<std::complex<T>> Q;	//!< basis
	matD<std::complex<T>> R;	//!< basis
};

//! \relates matD
template <typename T>
retMcMc<T> Schur(const matD<std::complex<T>>& A)
//! returns a Schur decomposition of a complex matrix A.
{	assert(A.M == A.N); const unsigned int n = A.N;
	matD<std::complex<T>> S = A, Z(n,n); vecD<std::complex<T>> w(n);
	gees(S, w, Z); return { Z,S };
}

template <typename T> 
struct retMcMcMcMc {
	matD<T> Sa;			//!< left matrix pair
	matD<T> Vl;			//!<
	matD<T> Sb;			//!< eight matrix pair
	matD<T> Vr;			//!<
};

//! \relates matD
template <typename T>
retMcMcMcMc<T> qz(const matD<T>& A, const matD<T>& B)
//! returns a generalized Schur decomposition of matrix pair (A,B)
//! as left and right matrix pairs (Sa,Vl), (Sb,Vr).
{	assert(A.M == A.N && B.M == B.N && A.M == B.N); const unsigned int n = A.N;
	matD<T> Sa = A, Sb = B, Vl(n,n), Vr(n,n);
	gges(Sa, Sb, Vl, Vr); return { Sa,Vl,Sb,Vr };
}

//! \relates matD
template <typename T>
void exchSchur(matD<T>& A, matD<T>& Q, unsigned int si, unsigned int di)
//! exchanges blocks (si,di) in a Schur decomposition (A,Q).
{	assert(A.M == A.N);
	trexc(A, Q, si+1, di+1);							// adjust for Fortran indices
}

//! \relates matD
template <typename T>
void exchGeneralSchur(matD<T>& Sa, matD<T>& Sb, matD<T>& Vl, matD<T>& Vr,
	const unsigned int si, const unsigned int di)
//! exchanges blocks (si,di) in a generalized Schur decomposition (Sa,Vl), (Sb,Vr).
{	assert(Sa.M == Sa.N && Sb.M == Sb.N && Sa.M == Sb.N);
	tgexc(Sa, Sb, Vl, Vr, si+1, di+1);						// adjust for Fortran indices
}

//! \relates matD
template <typename T>
void transformEigenvectors(matD<T>& S, matD<T>& V)
//! computes eigenvectors V of a Schur decomposition S.
{	trevc(S, V);
}

//! \relates matD
template <typename T>
retMM<T> Hessenberg(const matD<T>& A)	
//! computes a Hessenberg decomposition of a real symmetric matrix A in Q, returns H.
{ 	assert(A.M == A.N); const unsigned int n = A.N; vecD<T> tau(n-1);
	matD<T> H = A, Q(n,n); Q = T(0); gehrd(H, tau);
	for (unsigned int j = 1; j < n; j++) {						// disentangle H and Q
		for (unsigned int i = 0; i < j-1; i++) {
			Q(j,i) = H(j,i); H(j,i) = T(0); } };	  
	orghr(Q, tau); return { Q,H };
}

//! \relates matD
template <typename T>
retMcVc<T> hev(const matD<T>& A)
//! computes an eigen-decomposition of a real Hessenberg matrix A in U, returns eigenvectors.
{	assert(A.M == A.N); const unsigned int n = A.N;
	matD<T> S = A, U, V(n,n); vecD<T> dr(n), di(n);
	hseqr(S, dr, di, V); trevc(S, V);
	V.scaleEigenvectors(di).splitConjugatePairs(U,di);
	return { U,asComplex(dr,di) };
}

//! \relates matD
template <typename T>
retMcVc<T> hev(const matD<std::complex<T>>& A)
//! computes an eigen-decomposition of a complex Hessenberg matrix A in U, returns eigenvectors.
{	assert(A.M == A.N); const unsigned int n = A.N; vecD<std::complex<T>> d(n);
	matD<std::complex<T>> S = A, U(n,n); U = std::complex<T>(0);
	hseqr(S, d, U); trevc(S, U); 
	U.scaleEigenvectors(); return { U,d };
}

//! \relates matD
template <typename T>
matD<T> outer(const vecD<T>& a, const vecD<T>& b)
//! returns outer product of vectors (a,b) as a matrix.
{ 	matD<T> A(a.N,b.N);
	for (unsigned int i = 0; i < A.M; i++)
		for (unsigned int j = 0; j < A.N; j++) A(i,j) = a(i)*b(j);
	return A;
}

//! \relates matD
template <typename T>
matD<T> outer(const vec3<T>& a, const vec3<T>& b)
//! returns outer product of 3D vectors (a,b) as a matrix.
{ 	matD<T> A(3,3);
	A(0,0) = a.x*b.x; A(0,1) = a.x*b.y; A(0,2) = a.x*b.z; 
	A(1,0) = a.y*b.x; A(1,1) = a.y*b.y; A(1,2) = a.y*b.z; 
	A(2,0) = a.z*b.x; A(2,1) = a.z*b.y; A(2,2) = a.z*b.z; 
	return A;
}

//! \relates matD
template <typename T>
matD<T> CholDec(const matD<T>& A)
//! returns upper triangle of a Cholesky decomposition of the PSD matrix A.
{	assert(A.M == A.N); matD<T> U = A; potrf(U); return U;
}

//! \relates matD
template <typename T>
matD<T> CholInv(const matD<T>& A)
//! computes inverse of a Cholesky decomposition, returns upper triangle.
{	assert(A.M == A.N); matD<T> Ai = A; trtri(Ai); 
	for (unsigned int j = 0; j < Ai.N; j++)
		for (unsigned int i = j+1; i < Ai.M; i++) Ai(i,j) = T(0);
	return Ai;
}

//! \relates matD
template <typename T>
vecD<T> CholSolve(const matD<T>& A, const vecD<T>& b)
//! returns upper triangle of a Cholesky decomposition of the PSD matrix A.
{	assert(A.M == A.N); matD<T> U = A; vecD<T> x = b;
	potrs(U,x); return x;
}

//! \relates matD
template <typename T>
void balance(matD<T>& A, int& ilo, int& ihi, uvecD& pm, vecD<T>& sc)
//! balances real matrix A.
{	assert(A.M == A.N); pm.resize(A.N); sc.resize(A.N);
	gebal(A, ilo, ihi, pm, sc);
}

//! \relates matD
template <typename T, typename U>
void balance(matD<U>& A, int& ilo, int& ihi, uvecD& pm, vecD<T>& sc)
//! balances complex matrix A.
{	assert(A.M == A.N); pm.resize(A.N); sc.resize(A.N);
	gebal(A, ilo, ihi, pm, sc);
}

//! \relates matD
template <typename T>
matD<T>& matD<T>::joinBelow(const matD<T>& A)
//! joins matrix A below this matrix.
{	assert(A.N == N); const unsigned int nm = M+A.M; T* nx = new T [nm*N];		// allocate storage for new number of rows
	for (unsigned int c = 0; c < N; c++) memcpy(nx+nm*c,x+M*c,M*sizeof(T));		// copy columns of this matrix
	for (unsigned int c = 0; c < N; c++) memcpy(nx+nm*c+M,A.x+A.M*c,A.M*sizeof(T));	// copy columns of A below
	delete [] x; x = nx; M = nm; return *this;					// update this matrix
}

//! \relates matD
template <typename T>
matD<T>& matD<T>::joinRight(const matD<T>& A)
//! joins matrix A besides this matrix.
{	assert(A.M == M); const unsigned int nn = N+A.N; T* nx = new T [nn*M];		// allocate storage for new number of columns
	memcpy(nx,x,M*N*sizeof(T)); memcpy(nx+M*N,A.x,A.M*A.N*sizeof(T));		// copy this matrix and join A besides
	delete [] x; x = nx; N = nn; return *this;					// update this matrix
}

//! \relates matD
template <typename T>
matD<T>	symmetrize(const matD<T>& A)
//! returns 0.5*(A^T+A).
{	return (A+trp(A))*T(0.5);
}

//! \relates matD
template <typename T>
matD<T>	inner(const matD<T>& A, const matD<T>& B)
//! computes inner product of matrices (A,B).
{	matD<T> M(A.N,B.N); assert(A.M == B.M); const unsigned int n = A.M;
	for (unsigned int j = 0; j < B.N; j++) { const T* b = B.x+j*B.M;
		for (unsigned int i = 0; i < A.N; i++) { const T* a = A.x+i*A.M;
			T s{0}; for (unsigned int k = 0; k < n; k++) s += a[k]*b[k];
			M(i,j) = s; } };
	return M;
}

//! \relates matD
template <typename T>
void gs(matD<T>& W, const matD<T>& M, const matD<T>& Q)
//! performs Gram-Schmidt orthonormalization of W against basis Q.
{	for (unsigned int i = 0; i < W.N; i++) { vecD<T> w = W.getCol(i);
		T nm = norm(w), ni = nm == T(0)? T(0): T(1)/nm;
		w *= ni; vecD<T> mw = M*w;
		for (unsigned int j = 0; j < Q.N; j++) {
			const vecD<T> v = Q.getCol(j); w -= v*dot(v,mw); };
		nm *= norm(w); W.setCol(i,w*nm); };
}

//! \relates matD
template <typename T>
int svqb(matD<T>& W, const matD<T>& MW)
//! applies and SVQB decomposition to matrix W and image M*W. Returns 0 (ok), -1 (breakdown), 1 (no convergence).
{	matD<T> A = symmetrize(inner(W,MW)); const unsigned int n = A.M; vecD<T> d(n);	// compute S = W^T M W
	for (unsigned int i = 0; i < n; i++) {						// safely invert diagonal
		T v = A(i,i); d[i] = (v == T(0))? T(0): T(1)/std::sqrt(v); };
	for (unsigned int i = 0; i < n; i++)						// scale S
		for (unsigned int j = 0; j < n; j++) A(i,j) *= d(i)*d(j);
	auto S = sev(A);								// decompose S
	const T tau = T(10)*std::numeric_limits<T>::epsilon()*S.d(n-1);
	if (S.d(0) < tau) return -1;							// check for breakdown
	const T c = (S.d(n-1)/S.d(0))-T(1);
	for (unsigned int i = 0; i < n; i++) {						// safely invert eigenvalues
		T v = S.d(i); S.d[i] = (v == T(0))? T(0): T(1)/std::sqrt(v); };
	for (unsigned int i = 0; i < n; i++)						// scale eigenvectors
		for (unsigned int j = 0; j < n; j++) S.U(i,j) *= d(i)*S.d(j);
	W = W*S.U; return c < T(2)*tau? 0: 1;						// apply decomposition and check for convergence
}

//! \relates matD
template <typename T>
void ortho(matD<T>& W, const matD<T>& M, const matD<T>& Q)
//! performs SQVB orthonormalization of W against basis Q.
{	for (unsigned int it = 0; it < 100; it++) {
		gs(W,M,Q); gs(W,M,Q);							// perform two rounds of GS orthogonalization
		const int c = svqb(W,M*W); if (c == 0) return;				// orthonormalize W
		if (c == -1) random(W);	}						// breakdown in svqb, reinit W
}

//! \relates matD
template <typename T>
matD<T> centerCols(const matD<T>& A)
//! subtracts column-wise means from A.
{	const vecD<T> mn = A.colsum()/T(A.M); matD<T> B = A;
	for (unsigned int j = 0; j < B.N; j++) { T mi = mn(j), *x = B.x+j*B.M;
		for (unsigned int i = 0; i < B.M; i++) x[i] -= mi; };
	return B;
}

//! \relates matD
template <typename T>
matD<T>	normalizeCols(const matD<T>& A)
//! normalizes columns to zero mean and unit variance.
{	matD<T> B = A; 
	for (unsigned int j = 0; j < B.N; j++) { T m{0}, v{0}, *x = B.x+j*B.M;
		for (unsigned int i = 0; i < B.M; i++) m += x[i];
		m = isSmall<T>(m)? 0: m/B.M;
		for (unsigned int i = 0; i < B.M; i++) v += SQR(x[i]-m);
		v = B.M > 1? v/T(B.M-1): v; v = T(1)/std::sqrt(v);
		for (unsigned int i = 0; i < B.M; i++) x[i] = (x[i]-m)*v; };
	return B;
}

//! \relates matD
template <typename T>
matD<T>	pca(const matD<T>& A, const unsigned int l)
//! computes the PCA score matrix from a normalized matrix A.
{	matD<T> U,V; vecD<T> s = svd(A,U,V);
	const unsigned int n = std::min(s.N,l); const underscore _;
	for (unsigned int j = 0; j < n; j++) { T* x = U.x+j*U.M;
		for (unsigned int i = 0; i < U.M; i++) x[i] *= s[j]; };
	return U(_,_(0,n));
}

//! \relates matD
template <typename T>
T absdiff(const matD<T>& A, const matD<T>& B)
//! returns absolute difference of A and B.
{	T e{0}, *ax = A.x, *bx = B.x;
	for (unsigned int i = 0; i < A.M*A.N; i++) e += std::abs(*ax++-*bx++);
	return e;
}

//! \relates matD
template <typename T>
matD<T>	bistNormalize(const matD<T>& A, const unsigned int maxit = 1000, const T eps = T(1e-3))
//! performs Sinkhorn's bi-stochastic normalization (Section 2.4).
{	matD<T> B; const underscore _;
	if (A.M == A.N) B = A;								// symmetric problem
	else { const unsigned int n = std::max(A.M,A.N);				// else add slack variables
		B.resize(n,n); B = T(1); B.set(A); };
	for (unsigned int it = 0; it < maxit; it++) { matD<T> C = B;	 		// for nit iterations
		for (unsigned int i = 0; i < C.M; i++) { T s{0};			// normalize rows
			for (unsigned int j = 0; j < C.N; j++) s += C(i,j);
			assert(isSmall<T>(s) == false);
			for (unsigned int j = 0; j < C.N; j++) C(i,j) /= s; };
		for (unsigned int j = 0; j < C.N; j++) { T s{0};			// normalize columns
			for (unsigned int i = 0; i < C.M; i++) s += C(i,j);
			assert(isSmall<T>(s) == false);
			for (unsigned int i = 0; i < C.M; i++) C(i,j) /= s; };
		const T e = absdiff(B,C); B = C; if (e < eps) break; };
	return A.M == A.N? B: B(_(0,A.M),_(0,A.N));
}

//! \relates matD
template <typename T>
uvecD greedyMax(const matD<T>& A)
//! returns vector of column-wise maxima - greedy version.
{	matD<T> B = A; uvecD a(B.N); a = ~0u;						// copy matrix and init result vector
	while (true) { const unsigned int imax = B.imax();				// repeat: get index of maximum value in matrix
		if (B.x[imax] == T(0)) break;						// until matrix is empty
		const unsigned int r = imax%B.M, c = imax/B.M; a[c] = r;		// convert row & column, save result
		for (unsigned int j = 0; j < B.N; j++) B(r,j) = T(0);			// clear this column
		for (unsigned int i = 0; i < B.M; i++) B(i,c) = T(0); };		// clear this row
	return a;									// return column-wise maxima
}

//! \relates matD
template <typename T>
vecD<T> asVector(const matD<T>& A)
//! returns vector of columns of A.
{	const unsigned int n = A.M*A.N; vecD<T> a(n);
	for (unsigned int i = 0; i < n; i++) a.x[i] = A.x[i];
	return a;
}

//! \relates matD
template <typename T>
matD<T> asMatrix(const vecD<T>& a, const unsigned int m, const unsigned int n)
//! returns matrix (m,n) from vector a.
{	assert(m*n == a.N); matD<T> A(m,n);
	for (unsigned int i = 0; i < n; i++) A.x[i] = a.x[i];
	return A;									// return column-wise maxima
}

//! \relates matD
template <typename T>
matD<T> sqrt(const matD<T>& A)
//! returns square root of matrix A. NB: use real if A is psd, else complex.
{	assert(A.M == A.N); const unsigned int n = A.M;
	const auto S = Schur(A); matD<T> R(n,n); R = T(0);
	for (unsigned int j = 0; j < n; j++) { const T s = S.R(j,j);
		if (!isSmall<float,T>(s)) R(j,j) = std::sqrt(s);			// be generous here, use float.
		for (unsigned int i = j-1; i != ~0u; i--) { T s{0};
			for (unsigned int k = i+1; k < j; k++) s += R(i,k)*R(k,j);
			const T r = R(i,i)+R(j,j);
			if (std::abs(r) > T(0)) R(i,j) = (S.R(i,j)-s)/r; } };
	return S.Q*R*adj(S.Q);
}

template <typename T>
matD<T> invSqrt(const matD<T>& A, const T eps = std::numeric_limits<T>::epsilon())
//! returns A^{-0.5}. Works for A less than full rank.
{	assert(A.M == A.N); const unsigned int n = A.N;
	const auto S = sev(A);
	matD<T> D(n,n); D = T(0); const T lim = S.d(S.d.N-1)*eps;
	for (unsigned int i = 0; i < n; i++) {
		if (S.d(i) > lim) D(i,i) = T(1)/std::sqrt(S.d(i)); };
	return S.U*D;
}

//! \relates matD
template<typename T>
matD<T> minimumSpanningTree(const matD<T>& A)
//! returns minimum spanning tree using Kruskal's algorithm - undirected version.
{	class djSets {
		const unsigned int n;
		uvecD	par;
		uvecD	rnk;
	public:
		djSets(const unsigned int _n)
			: n(_n), par(n+1), rnk(n+1)
			{ rnk = 0; for (unsigned int i = 0; i <= n; i++) par[i] = i; }
		unsigned int find(const unsigned int u)					//! returns the parent of a node u.
			{ if (u != par[u]) { par[u] = find(par[u]); }; return par[u]; }
		void	merge(unsigned int x, unsigned int y)				//! performs union by rank
			{ x = find(x); y = find(y);
			  if (rnk[x] > rnk[y]) par[y] = x; else par[x] = y;		// make tree with smaller height a subtree of the other tree
			  if (rnk[x] == rnk[y]) rnk[y]++; }
	};

	assert(A.M == A.N); djSets ds(A.M);						// create disjoint sets
	std::vector<std::pair<T,std::pair<unsigned int,unsigned int>>> ed;		// represents an edge
	for (unsigned int i = 0; i < A.M; i++)						// create edges from dense matrix
		for (unsigned int j = i+1; j < A.M; j++)				// assume symmetric matrix
			if (A(i,j) != T(0)) ed.push_back({-A(i,j),{i,j}});		// use negative sign for highest weight
	std::sort(ed.begin(),ed.end()); matD<T> B(A.M,A.N); B.id();			// allocate output matrix
	for (const auto& e : ed) {
		const unsigned int u = e.second.first, v = e.second.second;		// get nodes
		const unsigned int su = ds.find(u), sv = ds.find(v);			// get sets
		if (su == sv) continue;        						// skip if this edge is creating a cycle
		B(u,v) = -e.first; B(v,u) = -e.first; ds.merge(su,sv); };		// else to the MST and merge two sets
	return B;
}

#undef allElem
#endif
