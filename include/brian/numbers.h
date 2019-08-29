#ifndef XX_H
#define XX_H

using fcomplex = std::complex<float>;
using dcomplex = std::complex<double>;

template<typename T>
T safeMin()
{	T eps = std::numeric_limits<T>::epsilon();
	T min = std::numeric_limits<T>::min();
	T base = std::numeric_limits<T>::radix;
	return T(pow(base,1.0+std::floor(std::log(min/eps)/std::log(base))/2.0));
}

template<typename T>
bool almostEqual(const T a, const T b,
	const T dabs = std::numeric_limits<T>::epsilon(),
	const T drel = std::numeric_limits<T>::epsilon())
{	// see C. Ericson, Real-Time Collision Detection
	return std::abs(a-b) <= std::max(dabs, drel*std::max(std::abs(a),std::abs(b)));
}

template<typename T, typename U = T>
bool isSmall(const U a, const T tol = std::numeric_limits<T>::epsilon()) 
{	return T(std::abs(a)) <= tol;
}

template<typename T>
inline bool checkNZ(const T a)
//! check if a is greater than epsilon.
{ return std::abs(a) > std::numeric_limits<T>::epsilon(); }

template<typename T>
inline bool checkEPS(const T a)
//! checks if a is greater than the safe minimum.
{ return std::abs(a) > safeMin<T>(); }

template<typename T>
inline void assertNZ(const T a)
//! asserts that a is greater than epsilon.
{ assert(checkNZ(a)); }

template <typename T>
void printT(const T& v)									// prints a templated value.
{	if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
		printf("%12.5e ", double(v));
	else if (typeid(T) == typeid(int))
		printf("%4d ", int(v));
	else if (typeid(T) == typeid(unsigned int))
		printf("%4u ", (unsigned int)v);
	else if (typeid(T) == typeid(bool))
		printf("%-5s ", v? "true": "false");
	else throw rtException("Cannot print type %s", typeid(T).name());
}

template <typename T>
void printT(const std::complex<T>& v)							// prints a templated complex value.
{	if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
		printf("(%12.5e %12.5ei) ", real(v), imag(v));
	else throw rtException("Cannot print type %s", typeid(T).name());
}

template <typename T>
void sprintT(char* buf, const T& v)							// prints a templated value into buffer.
{	if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
		sprintf(buf, "%12.5e ", double(v));
	else if (typeid(T) == typeid(int))
		sprintf(buf, "%4d ", int(v));
	else if (typeid(T) == typeid(unsigned int))
		sprintf(buf, "%4u ", (unsigned int)v);
	else if (typeid(T) == typeid(bool))
		sprintf(buf, "%-5s ", v? "true": "false");
	else throw rtException("Cannot print type %s", typeid(T).name());
}

template <typename T>
void sprintT(char* buf, const std::complex<T>& v)					// prints a templated complex value into buffer.
{	if (typeid(T) == typeid(double) || typeid(T) == typeid(float))
		sprintf(buf, "(%12.5e %12.5ei) ", real(v), imag(v));
	else throw rtException("Cannot print type %s", typeid(T).name());
}

inline float conj(const float a) { return a; }		// need specialization for real types
inline double conj(const double a) { return a; }
inline float real(const float a) { return a; }		// need specialization for real types
inline double real(const double a) { return a; }
inline float norm(const float a) { return a*a; }	// need specialization for real types
inline double norm(const double a) { return a*a; }
inline float imag(const float) { return 0; }		// need specialization for real types
inline double imag(const double) { return 0; }
inline float arg(const float a) { return a; }		// need specialization for real types
inline double arg(const double a) { return a; }

inline float mean(const float* x, const unsigned int n)
{ float s = 0; for (unsigned int i = 0; i < n; i++) s += x[i]; return s/float(n); }
inline double mean(const double* x, const unsigned int n)
{ double s = 0; for (unsigned int i = 0; i < n; i++) s += x[i]; return s/float(n); }
inline float var(const float* x, const unsigned int n)
{ float s2 = 0, m = mean(x,n); for (unsigned int i = 0; i < n; i++) s2 += (x[i]-m)*(x[i]-m); return s2/float(n-1); }
inline double var(const double* x, const unsigned int n)
{ double s2 = 0, m = mean(x,n); for (unsigned int i = 0; i < n; i++) s2 += (x[i]-m)*(x[i]-m); return s2/float(n-1); }
inline float dot(const float* x, const float* y, const unsigned int n)
{ float s = 0; for (unsigned int i = 0; i < n; i++) s += x[i]*y[i]; return s; }
inline double dot(const double* x, const double* y, const unsigned int n)
{ double s = 0; for (unsigned int i = 0; i < n; i++) s += x[i]*y[i]; return s; }
inline void abs(float* x, const unsigned int n)
{ for (unsigned int i = 0; i < n; i++) x[i] = std::abs(x[i]); }

//! Implements a range record for accessing matrices

class range {
public:
	unsigned int st;				//!< start index
	unsigned int end;				//!< end index
	range()										//! allocates an empty range.
		 : st(0), end(0) { }
	range(const unsigned int _st, const unsigned int _end)				//! allocates a range from st to end.
		 : st(_st), end(_end) { }
};

//! Implements a range operator for accessing matrices

class underscore {
public:
        underscore()									//! allocates an empty range operator.
		{ }
        range operator()(const unsigned int st, const unsigned int end) const		//! allocates a range operator from st to end.
		{ return range(st, end); }
};

#endif
