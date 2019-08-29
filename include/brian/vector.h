#ifndef VECTOR_H
#define VECTOR_H

/*
 *
 * vector.h: templated vector classes
 * BRIAN Software Package Version 3.0
 *
 * $Id: vector.h 438 2016-11-10 00:58:28Z frithjof $
 *
 * 0.10 (17/04/09): initial version
 * 0.11 (28/09/09): vecD and vecS added
 * 0.12 (01/06/10): ssq, imin, imax added
 * 0.20 (10/12/11): let vecS grow
 * 0.21 (08/06/12): changed for non-member functions, check complex
 * 0.22 (21/06/12): vecS dot product revised
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (20/12/13): documented
 * 0.50 (11/11/14): vecS pack, unpack & size added
 * 0.51 (16/10/15): vecD rank added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for dense and sparse vectors.
*/

template <class T> class uniformDist;		// forward declaration for random()

#define unaryOpScalar3(sym) \
	vec3<T>& operator sym (const T b) { x sym b; y sym b; z sym b; return *this; }
#define binaryOpScalar3(sym) \
	vec3<T> operator sym (const T b) const { return vec3<T>(x sym b, y sym b, z sym b); }
#define unaryOpVec3(sym) \
	vec3<T>& operator sym (const vec3<T>& b) { x sym b.x; y sym b.y; z sym b.z; return *this; }
#define binaryOpVec3(sym) \
	vec3<T> operator sym (const vec3<T>& b) const { return vec3<T>(x sym b.x, y sym b.y, z sym b.z); }

//! Implements a templated vector in 3D

template<typename T> class vec3 {
public:
	T	x;				//!< component x
	T	y;				//!< component y
	T	z;				//!< component z

	vec3<T>()									//! allocates an empty 3D vector.
		 : x(T(0)), y(T(0)), z(T(0)) { }
	vec3<T>(const T _x, const T _y, const T _z)					//! allocates a 3D vector at (x,y,z).
		 : x(_x), y(_y), z(_z) { }
	vec3<T>(const T b)								//! allocates a 3D vector at (b,b,b).
		 : x(b), y(b), z(b) { }
	template<typename U>
	vec3<T>(const vec3<U>& b)							//! copies from vector b - note type change.
		 : x(T(b.x)), y(T(b.y)), z(T(b.z)) { }
	template<typename U>
	vec3<T>(vec3<U>&& b)								//! moves from vector b - note type change.
		 : x(T(b.x)), y(T(b.y)), z(T(b.z)) { }
	template<typename U>
	vec3<T>& operator=(const vec3<U>& b)						//! assigns from vector b - note type change.
		{ x = T(b.x); y = T(b.y); z = T(b.z); return *this; }
	template<typename U>
	vec3<T>& operator=(vec3<U>&& b)							//! move assigns from vector b - note type change.
		{ x = T(b.x); y = T(b.y); z = T(b.z); return *this; }
	bool	operator==(const vec3<T>& b) const					//! equal if all components are equal.
		{ T d = SQR(x-b.x)+SQR(y-b.y)+SQR(z-b.z); return isSmall<T>(d); }
	bool	operator!=(const vec3<T>& b) const					//! unequal if any component is unequal.
		{ T d = SQR(x-b.x)+SQR(y-b.y)+SQR(z-b.z); return isSmall<T>(d) == false; }
	bool	operator<(const vec3<T>& b) const					//! unequal if any component is unequal.
		{ return norm2(*this) < norm2(b); }
		unaryOpScalar3(=)
		unaryOpScalar3(+=)
		unaryOpScalar3(-=)
		unaryOpScalar3(*=)
		unaryOpScalar3(/=)
		binaryOpScalar3(+)
		binaryOpScalar3(-)
		binaryOpScalar3(*)
		binaryOpScalar3(/)
		unaryOpVec3(+=)
		unaryOpVec3(-=)
		unaryOpVec3(*=)								// component-wise multiplication
		unaryOpVec3(/=)								// component-wise division
		binaryOpVec3(+)
		binaryOpVec3(-)
		binaryOpVec3(*)								// component-wise multiplication
		binaryOpVec3(/)								// component-wise division
        vec3<T>	max(const vec3<T>& b) const						//! returns vector of largest components.
		{ return vec3<T>((x > b.x? x: b.x), (y > b.y? y: b.y), (z > b.z? z: b.z)); }
        vec3<T>	min(const vec3<T>& b) const						//! returns vector of smallest components.
		{ return vec3<T>((x < b.x? x: b.x), (y < b.y? y: b.y), (z < b.z? z: b.z)); }
        T	max() const								//! returns largest component.
		{ return x > y? (x > z? x: (y > z? y: z)): (y > z? y: z); }
        T	min() const								//! returns smallest component.
		{ return x < y? (x < z? x: (y < z? y: z)): (y < z? y: z); }
        T	sum() const								//! returns sum of components.
		{ return x+y+z; }
	vec3<T>	normalize() const							//! returns normalized vector.
		{ T s = norm(*this); return s > T(0.0)? vec3<T>(x/s, y/s, z/s): *this; }
	bool	finite() const								//! returns true if all elements are finite.
		{ return std::isfinite(x) && std::isfinite(y) && std::isfinite(z); }
	void	read(is& in)
		{ in.read(x); in.read(y); in.read(z); }
	void	save(os& out) const
		{ out.save(x); out.save(y); out.save(z); }
};

//! \relates vec3
template <>
inline bool vec3<unsigned int>::operator== (const vec3<unsigned int>& b) const
//! returns true if two 3D unsigned vectors are equal.
{	return x == b.x && y == b.y && z == b.z; }

//! \relates vec3
template <>
inline bool vec3<int>::operator== (const vec3<int>& b) const
//! returns true if two 3D signed vectors are equal.
{	return x == b.x && y == b.y && z == b.z; }

//! \relates vec3
template <>
inline bool vec3<unsigned int>::operator!= (const vec3<unsigned int>& b) const
//! returns true if two 3D unsigned vectors are not equal.
{	return x != b.x || y != b.y || z != b.z; }

//! \relates vec3
template <>
inline bool vec3<unsigned int>::operator< (const vec3<unsigned int>& b) const
//! returns true if this vector is closer to the origin than b.
{	return (x > b.x || y > b.y)? false: z < b.z; }

//! \relates vec3
template <>
inline bool vec3<int>::operator!= (const vec3<int>& b) const
//! tests if two 3D signed vectors are not equal.
{	return x != b.x || y != b.y || z != b.z; }

//! \relates vec3
template <typename T>
vec3<T> operator*(const T a, const vec3<T>& b) 
//! multiplies scalar a and vector b.
{	return b * a;
}

template <typename T>
void print(const vec3<T>& b)
//! prints vector b on stdout.
{	printT(b.x); printT(b.y); printT(b.z); printf("\n");
}

//! \relates vec3
template <typename T>
vec3<T> conj(const vec3<T>& b)
//! returns real conjugate of b - noop.
{	return b;
}

//! \relates vec3
template <typename T>
vec3<std::complex<T>> conj(const vec3<std::complex<T>>& b)
//! returns complex conjugate of b.
{	return vec3<std::complex<T>>(conj(b.x), conj(b.y), conj(b.z));
}

//! \relates vec3
template <typename T>
T dot(const vec3<T>& a, const vec3<T>& b) 
//! returns dot product of real vectors (a,b).
{	return (a.x*b.x+a.y*b.y+a.z*b.z);
}

//! \relates vec3
template <typename T>
T dot(const vec3<T>& a, const vec3<std::complex<T>>& b)
//! returns dot product of real vectors a and complex vector b.
{	return (a.x*b.x.conj()+a.y*b.y.conj()+a.z*b.z.conj());
}

//! \relates vec3
template <typename T>
T dot(const vec3<std::complex<T>>& a, const vec3<T>& b)
//! returns dot product of complex vectors a and real vector b.
{	return (a.x.conj()*b.x+a.y.conj()*b.y+a.z.conj()*b.z);
}

//! \relates vec3
template <typename T>
T norm2(const vec3<T>& a)
//! returns infinite norm.
{	return SQR(a.x)+SQR(a.y)+SQR(a.z);
}

//! \relates vec3
template <typename T>
T norm(const vec3<T> a)
//! returns 2-norm.
{	return T(std::sqrt(norm2(a)));
}

//! \relates vec3
template <typename T>
vec3<T> cross(const vec3<T>& a, const vec3<T>& b)
//! returns cross product of real vectors (a,b).
{	return vec3<T>(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

//! \relates vec3
template <typename T> 
void random(vec3<T>& b)
//! generates real random vector -1...1.
{	uniformDist<T> ud(T(-1),T(1));
	b.x = ud(); b.y = ud(); b.z = ud();
}

//! \relates vec3
template<typename T>
T det3(const vec3<T>& a, const vec3<T>& b, const vec3<T>& c)
//! returns the determinant of a 3x3 matrix represented as column vectors.
{	return a.x*(b.y*c.z-b.z*c.y)-b.x*(a.y*c.z-a.z*c.y)+c.x*(a.y*b.z-a.z*b.y); }

#undef unaryOpScalar3
#undef binaryOpScalar3
#undef unaryOpVec3
#undef binaryOpVec3

using ivec3 = vec3<int>;
using uvec3 = vec3<unsigned int>;
using fvec3 = vec3<float>;
using dvec3 = vec3<double>;
using cvec3 = vec3<fcomplex>;
using zvec3 = vec3<dcomplex>;
using vuvec3 = std::vector<uvec3>;
using vfvec3 = std::vector<fvec3>;
using vdvec3 = std::vector<dvec3>;

// class vecD implements a templated vector of arbitrary length

#define allElem(i) (unsigned int i = 0; i < N; i++) 

#define unaryOpScalarD(sym) \
	vecD<T>& operator sym (const T b) \
	{ for allElem(i) x[i] sym b; return *this; }
#define binaryOpScalarD(sym) \
	vecD<T>	operator sym (const T b) const \
	{ vecD<T> a(N); for allElem(i) a[i] = x[i] sym b; return a; }
#define unaryOpVecD(sym) \
	vecD<T>& operator sym (const vecD<T>& b) \
	{ assert(N == b.N); for allElem(i) x[i] sym b(i); return *this; }
#define binaryOpVecD(sym) \
	vecD<T>	operator sym (const vecD<T>& b) const \
	{ assert(N == b.N); vecD<T> a(N); for allElem(i) a[i] = x[i] sym b(i); return a; }

//! Implements a templated dense vector

template<typename T> class vecD {
protected:
	void		alloc(const unsigned int _N)					//! allocates space for N elements.
			{ N = _N; x = new T [N]; }
	void		clear()								//! frees up space.
			{ delete [] x; x = nullptr; N = 0; }
	void		assign(const vecD<T>& b)					// assigns from vector b.
			{ if (N != b.N) { clear(); alloc(b.N); };
			  for allElem(i) x[i] = b(i); }
	vecD<T>		getRange(const unsigned int s, const unsigned int e) const	//! access range (s:e).
			{ assert(s <= N && e <= N && s <= e );
			  const unsigned int n = e-s; vecD<T> v(n);
			  for (unsigned int i = 0; i < n; i++) v(i) = (*this)(s+i);
			  return v; }
public:
	unsigned int N;				//!< dimension of vector
	T*	x;				//!< array of elements

	vecD<T>(const unsigned int _N = 0)						//! allocates a vector of N elements.
			: N(0), x(nullptr) { alloc(_N); }
	vecD(const vecD<T>& b)								//! copies from vector b.
			: N(0), x(nullptr) { assign(b); }
	vecD(vecD<T>&& b)								//! moves from vector b.
			: N(b.N), x(b.x) { b.x = nullptr; }
	~vecD<T>()
			{ clear(); }
	vecD<T>&	operator=(const vecD<T>& b)					//! assigns from vector b.
			{ if (this != &b) assign(b); return *this; }
	vecD<T>&	operator=(vecD<T>&& b)						//! move assigns from vector b.
			{ assert(this != &b); delete [] x; x = b.x; b.x = nullptr; N = b.N;
		  	  return *this; }
        T&		operator[](const unsigned int i)				//! returns reference to element i.
			{ return x[i]; }
        T		operator()(const unsigned int i) const				//! returns element i.
			{ return x[i]; }
        T&		operator()(const unsigned int i)				//! returns reference to element i.
			{ return x[i]; }
        vecD<T>		operator()(const range& r) const				//! returns row (r,_).
			{ return getRange(r.st,r.end); }
			unaryOpScalarD(=)
			unaryOpScalarD(+=)
			unaryOpScalarD(-=)
			unaryOpScalarD(*=)
			unaryOpScalarD(/=)
			binaryOpScalarD(+)
			binaryOpScalarD(-)
			binaryOpScalarD(*)
			binaryOpScalarD(/)
			unaryOpVecD(+=)
			unaryOpVecD(-=)
			unaryOpVecD(*=)
			unaryOpVecD(/=)
			binaryOpVecD(+)
			binaryOpVecD(-)
			binaryOpVecD(*)
			binaryOpVecD(/)
	T		sum() const							//! returns sum over all elements.
			{ T s{0}; for allElem(i) s += x[i]; return s; }
	vecD<T>&	normalize()							//! normalizes vector to sum of 1.
			{ T s = sum(); if (s == T(0)) return *this;
			  s = T(1)/s; for allElem(i) { x[i] *= s; }; return *this; }
	T		entropy() const							//! returns the entropy of this vector.
			{ T e{0}; for allElem(i) {
				if (x[i] > T(0)) e += x[i]*T(std::log(x[i])); };
			  return e; }
	T		chisq(const vecD<T>& b) const					//! computes sum-of-squared differences to vector b.
			{ assert(N == b.N); T s{0}; for allElem(i) s += SQR(x[i]-b.x[i]);
			  return s; }
	unsigned int	imax() const							//! returns position of maximum element.
			{ T m = std::numeric_limits<T>::lowest(); unsigned int k = 0;
			  for allElem(i) { if (x[i] > m) { m = x[i]; k = i; } };
			  return k; }
	unsigned int	amax() const							//! returns position of absolute maximum element.
			{ T m{0}; unsigned int k = 0; 
			  for allElem(i) { if (std::abs(x[i]) > m) { m = std::abs(x[i]); k = i; } };
			  return k; }
	unsigned int	imin() const							//! returns position of minimum element.
			{ T m = std::numeric_limits<T>::max(); unsigned int k = 0;
			  for allElem(i) { if (x[i] < m) { m = x[i]; k = i; } };
			  return k; }
	bool		resize(const unsigned int n)					//! resizes vector for n elements.
			{ clear(); alloc(n); return x != nullptr; }
	vecD<T>&	swap(vecD<T>& b)						//! swaps with vector b.
			{ std::swap(N,b.N); std::swap(x,b.x); return *this; }
	unsigned int	nnz() const							//! returns number of non-zero elements.
			{ unsigned int n = 0; for allElem(i) {
				if (std::abs(x[i]) != T(0)) n++; };
			  return n; }
	void		range(T& a, T& b) const						//! returns min & max in this vector.
			{ a = std::numeric_limits<T>::max(); b = -a;
			  for allElem(i) { b = std::max(b,x[i]); a = std::min(a,x[i]); } }
	vecD<T>		histogram(const float vmin, const float vmax,
				const unsigned int nbins = 256) const 			//! returns histogram of samples in this vector.
			{ vecD<T> hs(nbins); hs = T(0); for allElem(i) {
				const T a = round((nbins-1)*((x[i]-vmin)/(vmax-vmin)));
				if (a < T(0)) continue;
				const unsigned int b = FTOU(a); if (b < nbins) hs[b] += T(1); };
			  T s = T(1.0/hs.sum()); hs *= s; return hs; }
	vecD<T>		histogram(const unsigned int nbins = 256) const			//! computes normalized histogram of this vector.
			{ T vmin, vmax; range(vmin, vmax);
			  return histogram(vmin,vmax,nbins); }
	vecD<T>		quantiles() const						//! returns quantiles for the samples in this vector.
			{ std::vector<T> v(N); for allElem(i) v[i] = x[i];
			  std::sort(v.begin(), v.end()); vecD<T> q(7);
			  q[0] = v[0]; q[1] = v[N/40]; q[2] = v[N/4]; q[3] = v[N/2];
			  q[4] = v[3*N/4]; q[5] = v[39*N/40]; q[6] = v[N-1]; return q; }
	T		mean() const							//! returns sample mean.
			{ return ::mean(x,N); }
	T		var() const							//! returns sample variance.
			{ return ::var(x,N); }
	bool		finite() const							//! returns true if all elements are finite.
			{ for allElem(i) {
				if (std::isfinite(x[i]) == false) return false; };
			  return true; }
	vecD<T>		sample(const unsigned int n)					//! resample vector to size n.
			{ vecD<T> a(n); uniformDist<T> ud(0,N);
			  for (unsigned int i = 0; i < n; i++) a[i] = x[FTOU(ud())];
			  return a; }
	vecD<T>&	center()							//! normalizes to zero mean
			{ const T m = ::mean(x,N); for allElem(i) { x[i] -= m; };
			  return *this; }
	vecD<T>&	zmuv()								//! normalizes to zero mean and unit variance.
			{ center(); if (N <= 1) return *this;
			  T v{0}; for allElem(i) v += SQR(x[i]);
			  if (v == T(0)) return *this;
			  v = std::sqrt((N-1)/v); for allElem(i) x[i] *= v;
			  return *this; }
#undef allElem
};

template<typename T> vecD<std::complex<T>> fft(const vecD<T>&);				// compute a forward FFT of a float vector, see fft.C
template<typename T> vecD<T> ifft(const vecD<std::complex<T>>&);			// compute an inverse FFT of a complex vector, see fft.C

#define allElem(i) (unsigned int i = 0; i < b.N; i++)

//! \relates vecD
template<typename T, typename U> 
vecD<U> operator*(const vecD<T>& v, const U b)
//! multiplies real vector with complex scalar.
{	vecD<U> a(v.N); for (unsigned int i = 0; i < v.N; i++) a[i] = v(i)*b;
	return a;
}

//! \relates vecD
template <typename T>
T vvD(const T *ax, const T *bx, unsigned int N)
//! real vector-vector multiplication.
{	T s{0}; for (unsigned int i = 0; i < N; i++) s += ax[i]*bx[i];
	return s;
}

//! \relates vecD
template <typename T>
T vcD(const T *ax, const T *bx, unsigned int N)
//! complex vector-vector multiplication.
{	T s{0}; for (unsigned int i = 0; i < N; i++) s += conj(ax[i])*bx[i];
	return s;
}

//! \relates vecD
template <typename T>
vecD<T> operator*(const T a, const vecD<T>& b)
//! multiplies scalar with vector.
{	return b * a;
} 

//! \relates vecD
template <typename T>
void print(const vecD<T>& b)
//! prints real vector b to stdout.
{	for allElem(i) { printf("%5d ", i); printT(b(i)); printf("\n"); }
}

//! \relates vecD
template <typename T>
T norm1(const vecD<T>& b)
//! returns 1-norm of real vector b.
{	T s{0}; for allElem(i) s += std::abs(b(i));
	return s;
}

//! \relates vecD
template <typename T>
T norm2(const vecD<T>& b)
//! returns infinite norm of real vector b.
{	return dot(b,b);
}

//! \relates vecD
template <typename T>
T norm2(const vecD<std::complex<T>>& b)
//! returns real part of infinite norm of complex vector b.
{	return dot(b,b).real(); 
}

//! \relates vecD
template <typename T>
T norm(const vecD<T>& b)
//! returns norm of real vector b.
{	return T(std::sqrt(norm2(b)));
}

//! \relates vecD
template <typename T>
T norm(const vecD<std::complex<T>>& b)
//! returns norm of complex vector b.
{	return T(std::sqrt(norm2(b)));
}

//! \relates vecD
template <typename T> 
vecD<std::complex<T>> asComplex(const vecD<T>& b)
//! converts real to complex vector.
{	vecD<std::complex<T>> a(b.N);
	for allElem(i) a[i] = std::complex<T>(b(i));
	return a;
}

//! \relates vecD
template <typename T> 
vecD<std::complex<T>> asComplex(const vecD<T>& b, const vecD<T>& bi)
//! converts real & imaginary vectors to complex vector.
{	assert(b.N == bi.N); vecD<std::complex<T>> a(b.N);
	for allElem(i) a[i] = std::complex<T>(b(i),bi(i));
	return a;
}

//! \relates vecD
template <typename T> 
vecD<T> real(const vecD<std::complex<T>>& b)
//! returns real portion of complex vector.
{	vecD<T> a(b.N); for allElem(i) a[i] = b(i).real();
	return a;
}

//! \relates vecD
template <typename T> 
vecD<T> imag(const vecD<std::complex<T>>& b)
//! returns imaginary portion of complex vector.
{	vecD<T> a(b.N); for allElem(i) a[i] = b(i).imag();
	return a;
}

//! \relates vecD
template <typename T> 
void random(vecD<T>& b)
//! generates real random vector -1...1.
{	uniformDist<T> ud(T(-1),T(1));
	for allElem(i) b[i] = ud();
}

//! \relates vecD
template <typename T> 
void random(vecD<std::complex<T>>& b)
//! generates complex random vector -1...1.
{	uniformDist<T> ud(T(-1),T(1));
	for allElem(i) b[i] = std::complex<T>(ud(),ud());
}

//! \relates vecD
template <typename T> 
void sort(vecD<T>& b, const bool desc = false)
//! sorts vector b.
{	std::vector<T> v; v.assign(b.x,b.x+b.N);
	if (desc) std::sort(v.begin(),v.end(),std::greater<T>());			// descending sort
	else std::sort(v.begin(),v.end());						// ascending sort
	for allElem(i) b[i] = v[i];
}

//! \relates vecD
template <typename T> 
vecD<T> unique(const vecD<T>& b)
//! returns vector of unique elements in vector b.
{	std::vector<T> v; v.assign(b.x,b.x+b.N); std::sort(v.begin(),v.end());		// ascending sort
	const auto e = std::unique(v.begin(),v.end()); 
	const unsigned int n = std::distance(v.begin(),e);
	vecD<T> a(n); for (unsigned int i = 0; i < n; i++) a[i] = v[i];
	return a;
}

//! \relates vecD
template <typename T> 
vecD<unsigned int> rank(const vecD<T>& b, const bool desc = false)
//! returns vector with rank order 
{	std::vector<std::pair<T,unsigned int>> v(b.N); vecD<unsigned int> r(b.N);
	for (unsigned int i = 0; i < b.N; i++) v[i] = std::pair<T,unsigned int>(b(i),i);
	if (desc) std::sort(v.begin(),v.end(),std::greater<std::pair<T,unsigned int>>());
	else std::sort(v.begin(),v.end(),std::less<std::pair<T,unsigned int>>());
	for (unsigned int i = 0; i < b.N; i++) r[i] = v[i].second;
        return r;
}

//! \relates vecD
template <typename T> 
vecD<T> conj(const vecD<T>& b)
//! returns conjugate of real vector b - noop.
{	return b;
}

//! \relates vecD
template <typename T> 
vecD<std::complex<T>> conj(const vecD<std::complex<T>>& b)
//! returns conjugate of conjugate of complex vector b.
{	vecD<std::complex<T>> a(b.N);
	for allElem(i) a[i] = conj(b(i));
	return a;
}

//! \relates vecD
template <typename T> 
T dot(const vecD<T>& a, const vecD<T>& b)
//! returns dot product of real vectors (a,b).	
{	assert(a.N == b.N); return vvD(a.x,b.x,a.N);
}

//! \relates vecD
template <typename T> 
std::complex<T> dot(const vecD<std::complex<T>>& a, const vecD<std::complex<T>>& b)
//! returns dot product of complex vectors (a,b).	
{	assert(a.N == b.N); return vcD(a.x,b.x,a.N);
}

//! \relates vecD
template <typename T> 
std::complex<T> dot(const vecD<std::complex<T>>& a, const vecD<T>& b)
//! returns dot product of complex vector a and real vector b.
{	assert(a.N == b.N); return vcD(a.x,b.x,a.N);
}

//! \relates vecD
template <typename T> 
std::complex<T> dot(const vecD<T>& b, const vecD<std::complex<T>>& a)
//! returns dot product of real vector a and complex vector b.
{	assert(a.N == b.N); return vcD(a.x,b.x,a.N);
}

//! \relates matD
template <typename T>
vecD<T> stack(const vecD<T>& a, const vecD<T>& b)
//! stacks vector a above b.
{	vecD<T> r(a.N+b.N);
	for (unsigned int i = 0; i < a.N; i++) r[i] = a(i);
	for (unsigned int i = 0; i < b.N; i++) r[i+a.N] = b(i);
	return r;
}

using bvecD = vecD<bool>;
using ivecD = vecD<int>;
using uvecD = vecD<unsigned int>;
using fvecD = vecD<float>;
using dvecD = vecD<double>;
using cvecD = vecD<fcomplex>;
using zvecD = vecD<dcomplex>;
using vbvecD = std::vector<bvecD>;
using vuvecD = std::vector<uvecD>;
using vfvecD = std::vector<fvecD>;
using vdvecD = std::vector<dvecD>;

#undef allElem
#undef unaryOpScalarD
#undef binaryOpScalarD
#undef unaryOpVecD
#undef binaryOpVecD

#define allElem(i) (unsigned int i = 0; i < N; i++) 

#define unaryOpScalarS(sym) \
	vecS<T>& operator sym (const T b) \
		{ for allElem(i) x[i] sym b; return *this; }
#define binaryOpScalarS(sym) \
	vecS<T> operator sym (const T b) const \
		{ vecS v(*this); for allElem(i) v.x[i] = x[i] sym b; return v; }
#define unaryOpVecS(sym) \
	vecS<T>& operator sym (const vecS<T>& b) \
		{ for allElem(i) x[i] sym b.x[i]; return *this; }
#define binaryOpVecS(sym) \
	vecS<T> operator sym (const vecS<T>& b) const \
		{ vecS v(*this); for allElem(i) v.x[i] = x[i] sym b.x[i]; return v; }

//! Implements a templated sparse vector of length N

template<typename T> class vecS {
protected:
	void		alloc(const unsigned int _N)					//! allocates space for N non-zero elements.
			{ N = _N; u = new unsigned int [N]; x = new T [N]; }
	void		clear()								//! frees space.
			{ delete [] x; x = nullptr; delete [] u; u = nullptr; N = 0; }
	void		assign(const vecS<T>& b)					//! assigns from sparse vector b.
			{ if (N != b.N) { clear(); alloc(b.N); };
			  for allElem(i) { u[i] = b.u[i]; x[i] = b.x[i]; } }
	void		append(const unsigned int ul, T xl)				//! appends entry (u,x).
			{ unsigned int *un = new unsigned int [N+1]; T *xn = new T [N+1]; 
			  for allElem(i) { un[i] = u[i]; xn[i] = x[i]; };
			  un[N] = ul; xn[N] = xl; clear(); u = un; x = xn; N = N+1; }
public:
	unsigned int N;				//!< length of vector
	unsigned int* u;			//!< array of indices
	T*	x;				//!< array of values

	vecS<T>(const unsigned int _N = 0)						//! allocates a vector with N non-zero elements.
			: N(0), u(nullptr), x(nullptr) { alloc(_N); }
	vecS(const vecS<T>& b)								//! copies from sparse vector b.
			: N(0), u(nullptr), x(nullptr) { assign(b); }
	vecS(vecS<T>&& b)								//! moves from sparse vector b.
			: N(b.N), u(b.u), x(b.x) { b.u = nullptr; b.x = nullptr; }
	~vecS<T>() 	{ clear(); }
	vecS<T>&	operator=(const vecS<T>& b)					//! assigns from sparse vector b.
			{ if (this != &b) assign(b); return *this; }
	vecS<T>&	operator=(vecS<T>&& b)						//! move assigns from sparse vector b.
			{ assert(this != &b);
		 	 delete [] u; u = b.u; b.u = nullptr;
			  delete [] x; x = b.x; b.x = nullptr;
			  N = b.N; return *this; }
        T&		operator[](const unsigned int id)				//! returns LHS reference to element i.
			{ unsigned int t = 0; assert(find(id, t)); return x[t]; }
        T		operator()(const unsigned int id) const 			//! returns element i.
			{ unsigned int t = 0; assert(find(id, t)); return x[t]; }
        T&		operator()(const unsigned int id)				//! returns reference to element i
			{ unsigned int t = 0; assert(find(id, t)); return x[t]; }
			unaryOpScalarS(=)
			unaryOpScalarS(+=)
			unaryOpScalarS(-=)
			unaryOpScalarS(*=)
			unaryOpScalarS(/=)
			binaryOpScalarS(+)
			binaryOpScalarS(-)
			binaryOpScalarS(*)
			binaryOpScalarS(/)
			unaryOpVecS(+=)
			unaryOpVecS(-=)
			unaryOpVecS(*=)
			unaryOpVecS(/=)
			binaryOpVecS(+)
			binaryOpVecS(-)
			binaryOpVecS(*)
			binaryOpVecS(/)
	vecS<T>&	set(const unsigned int i, const unsigned int ui, const T xi) 	//! sets element i to (ui,xi).
			{ assert(i < N); u[i] = ui; x[i] = xi; return *this; }
	bool		resize(const unsigned int n)					//! resize vector for N non-zero elements.
			{ clear(); alloc(n); return u != nullptr && x != nullptr; }
	vecS<T>&	update(const unsigned int i, const T v) 			//! sets value of element i to v.
			{ unsigned int t = 0; if (find(i, t)) x[t] += v; else append(i,v);
			  return *this; }
	bool		find(const unsigned int id, unsigned int &t) const 		//! returns position t of element id.
			{ for allElem(i) if (id == u[i]) { t = i; return true; };
		 	  return false; }
	unsigned int	size() const							//! returns size (in bytes) for this vector.
			{ return (N+1)*sizeof(unsigned int)+N*sizeof(T); }
	unsigned int	unpack(const unsigned char* buf)				//! unpacks a vector from buffer* buf.
			{ const unsigned int n = *reinterpret_cast<const unsigned int*>(buf);
			  resize(n); unsigned int t = sizeof(unsigned int);
			  for allElem(i) { u[i] = *reinterpret_cast<const unsigned int*>(buf+t);
				t += sizeof(unsigned int); }
			  for allElem(i) { x[i] = *reinterpret_cast<const T*>(buf+t); t += sizeof(T); }
			  return t; }
	unsigned int	pack(unsigned char* buf) const					//! packs a vector into buffer* buf.
			{ *reinterpret_cast<unsigned int*>(buf) = N;
			  unsigned int t = sizeof(unsigned int);
			  for allElem(i) { *reinterpret_cast<unsigned int*>(buf+t) = u[i];
				t += sizeof(unsigned int); }
			  for allElem(i) { *reinterpret_cast<T*>(buf+t) = x[i]; t += sizeof(T); }
			  return t; }
	bool		finite() const							//! returns true if all elements are finite.
			{ for allElem(i) { if (std::isfinite(x[i]) == false) return false; };
			  return true; }
#undef allElem
};

#define allElem(i) (unsigned int i = 0; i < b.N; i++)

//! \relates vecS
template <typename T>
T vvS(const vecS<T>& a, const vecS<T>& b)
//! multiplies real sparse vectors (a,b).
{	T s{0}; unsigned int t = 0;
	for (unsigned int i = 0; i < a.N; i++) {
		if (b.find(a.u[i],t)) s += a.x[i]*b.x[t]; };
	return s;
}

//! \relates vecS
template <typename T>
std::complex<T> vcS(const vecS<std::complex<T>>& a, const vecS<T>& b)
//! multiplies complex sparse vector a with real sparse vector b.
{	std::complex<T> s{0}; unsigned int t = 0;
	for (unsigned int i = 0; i < a.N; i++) {
		if (b.find(a.u[i],t)) s += conj(a.x[i])*b.x[t]; };
	return s;
}

//! \relates vecS
template <typename T>
std::complex<T> vcS(const vecS<std::complex<T>>& a, const vecS<std::complex<T>>& b)
//! multiplies complex sparse vectors (a,b).
{	std::complex<T> s{0}; unsigned int t = 0;
	for (unsigned int i = 0; i < a.N; i++) {
		if (b.find(a.u[i],t)) s += conj(a.x[i])*b.x[t]; };
	return s;
}

//! \relates vecS
template <typename T>
vecS<T> operator* (const T a, const vecS<T>& b)
//! multiplies scalar a with real sparse vectors b.
{	return b*a;
}

//! \relates vecS
template <typename T>
void print(const vecS<T>& b)
//! prints vector b to stdout.
{	for allElem(i) { printf("%5d ", b.u[i]); printT(b.x[i]); printf("\n"); }; 
}

//! \relates vecS
template <typename T>
T norm2(const vecS<T>& b)
//! returns infinite norm of real vector b.
{	return dot(b,b);
}

//! \relates vecS
template <typename T>
T norm2(const vecS<std::complex<T>>& b)
//! returns infinite norm of complex vector b.
{	return real(dot(b,b));
}

//! \relates vecS
template <typename T>
T norm(const vecS<T>& b)
//! returns infinite norm of real vector b.
{	return T(std::sqrt(norm2(b)));
}

//! \relates vecS
template <typename T>
T norm(const vecS<std::complex<T>>& b)
//! returns infinite norm of complex vector b.
{	return T(std::sqrt(norm2(b)));
}

//! \relates vecS
template <typename T> 
vecS<std::complex<T>> asComplex(const vecS<T>& b)
//! converts real to complex vector.
{	vecS<std::complex<T>> a(b.N);
	for (unsigned int i = 0; i < a.N; i++) {
		a.u[i] = b.u[i]; a.x[i] = std::complex<T>(b(i)); };
	return a;
}

//! \relates vecS
template <typename T> 
vecS<T> real(const vecS<std::complex<T>>& b)
//! returns real portion of complex vector b.
{	vecS<T> a(b.N);
	for allElem(i) { a.u[i] = b.u[i]; a.x[i] = real(b(i)); };
	return a;
}

//! \relates vecS
template <typename T> 
vecS<T> imag(const vecS<std::complex<T>>& b)
//! return imaginary portion of complex vector.
{	vecS<T> a(b.N);
	for allElem(i) { a.u[i] = b.u[i]; a.x[i] = imag(b(i)); };
	return a;
}

//! \relates vecS
template <typename T> 
void random(vecS<T>& b)
//! generates real random vector -1...1.
{	uniformDist<T> ud(T(-1),T(1));
	for allElem(i) b[i] = ud.sample();
}

//! \relates vecS
template <typename T> 
void random(vecS<std::complex<T>>& b)
//! generates complex random vector -1...1.
{	uniformDist<T> ud(T(-1),T(1));
	for allElem(i) b[i] = std::complex<T>(ud.sample(),ud.sample());
}

//! \relates vecS
template <typename T> 
vecS<T> conj(const vecS<T>& b)
//! returns conjugate of real vector b - noop.
{	return b;
}
			
//! \relates vecS
template <typename T> 
vecS<std::complex<T>> conj(const vecS<std::complex<T>>& b)
//! returns conjugate of complex vector b.
{	vecS<std::complex<T>> a(b.N);
	for allElem(i) { a.u[i] = b.u[i]; a.x[i] = conj(b(i)); };
	return a;
}

//! \relates vecS
template <typename T> 
T dot(const vecS<T>& a, const vecS<T>& b)
//! returns dot product of real vectors (a,b).	
{	return vvS(a,b);
}

//! \relates vecS
template <typename T> 
std::complex<T> dot(const vecS<std::complex<T>>& a, const vecS<std::complex<T>>& b)
//! returns dot product of complex vectors (a,b).	
{	return vcS(a,b);
}

//! \relates vecS
template <typename T> 
std::complex<T> dot(const vecS<std::complex<T>>& a, const vecS<T>& b)
//! returns dot product of complex vectors a and real vector b.	
{	return vcS(a,b);
}

template <typename T> 
std::complex<T> dot(const vecS<T>& a, const vecS<std::complex<T>>& b)
//! returns dot product of real vectors a and complex vector b.		
{	return vcS(b,a);
}

//! \relates vecS
template <typename T>
vecD<T> asDense(const vecS<T>& a, const unsigned int n)
//! convert sparse to dense vector
{	vecD<T> d(n); d = T(0);
	for (unsigned int i = 0; i < a.N; i++) d[a.u[i]] = a.x[i];
	return d;
}

//! \relates vecS
template <typename T>
vecS<T> asSparse(const vecD<T>& a)
//! converts dense to sparse vector.
{	vecS<T> d(a.nnz()); unsigned int j = 0;
	for (unsigned int i = 0; i < a.N; i++) {
		if (std::abs(a.x[i]) != 0) d.set(j++, i, a.x[i]); };
	return d;
}

using fvecS = vecS<float>;
using dvecS = vecS<double>;
using cvecS = vecS<fcomplex>;
using zvecS = vecS<dcomplex>;

#undef allElem
#undef unaryOpScalarS
#undef binaryOpScalarS
#undef unaryOpVecS
#undef binaryOpVecS

#define unaryOpScalar2(sym) \
	svec<T>& operator sym (const T b) { u sym b; v sym b; return *this; }
#define binaryOpScalar2(sym) \
	svec<T> operator sym (const T& b) const { return svec<T>(u sym b, v sym b); }
#define unaryOpsvec(sym) \
	svec<T>& operator sym (const svec<T>& b) { u sym b.u; v sym b.v; return *this; }
#define binaryOpsvec(sym) \
	svec<T> operator sym (const svec<T>& b) const { return svec<T>(u sym b.u, v sym b.v); }

//! Implements a templated 2d spherical vector

template <typename T> class svec {
public:
	T	u;            			//!< u = theta, angle in x,y-plane, in [0,2pi[
	T	v;              		//!< v = phi, angle in z-direction, in [-pi/2, pi/2[
	void	toSphere(const vec3<T>& p)  						//! converts position p to spherical coordinates.
		{ toSphere(p.x,p.y,p.z); }
	void	toSphere(const T x, const T y, const T z) 				//! converts position (x,y,z) to spherical coordinates.
		{ u = std::atan2(y,x); v = std::acos(z); }
	vec3<T>	toCart() const 								//! converts spherical coordinates to position p.
		{ return vec3<T>(std::cos(u)*std::sin(v), std::sin(u)*std::sin(v), std::cos(v)); }
	svec<T>()									//! allocates an empty spherical vector.
		 : u(0), v(0) { }
	svec<T>(const T _u, const T _v)							//! allocates a spherical vector (u,v).
		 : u(_u), v(_v) { }
	template<typename U>
	svec<T>(const vec3<U>& p)							//! allocates a spherical vector from position p - allows type change.
		{ toSphere(p.x,p.y,p.z); }
	template<typename U>								//! allocates a spherical vector from position (x,y,z) - allows type change.
	svec<T>(const U x, const U y, const U z)
		{ toSphere(x,y,z); }
	vec3<T>	cart() const								//! converts to position p.
		{ return toCart(); }
		unaryOpScalar2(=)
		unaryOpScalar2(+=)
		unaryOpScalar2(-=)
		unaryOpScalar2(*=)
		unaryOpScalar2(/=)
		binaryOpScalar2(+)
		binaryOpScalar2(-)
		binaryOpScalar2(*)
		binaryOpScalar2(/)
		unaryOpsvec(+=)
		unaryOpsvec(-=)
		unaryOpsvec(*=)								// component-wise multiplication
		unaryOpsvec(/=)								// component-wise division
		binaryOpsvec(+)
		binaryOpsvec(-)
		binaryOpsvec(*)								// component-wise multiplication
		binaryOpsvec(/)								// component-wise division
	bool	operator==(const svec<T>& b) const					//! equal if both components equal.
		{ return u == b.u && v == b.v; }
	bool	operator!=(const svec<T>& b) const					//! unequal if any components unequal.
		{ return u != b.u || v != b.v; }
};

//! \relates svec
template <typename T>
svec<T> operator*(const T a, const svec<T>& b)
//! multiplies scalar a with spherical vector b.
{	return b * a;
}

//! \relates svec
template <typename T>
void print(const svec<T>& b)
//! prints spherical vector to stdout.
{	printT(b.u); printT(b.v);
}

using fsvec = svec<float>;
using dsvec = svec<double>;
using csvec = svec<fcomplex>;
using zsvec = svec<dcomplex>;
using vfsvec = std::vector<fsvec>;
using vdsvec = std::vector<dsvec>;

#undef unaryOpScalar2
#undef binaryOpScalar2
#undef unaryOpsvec
#undef binaryOpsvec

#endif

