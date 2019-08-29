#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

/*
 *
 * distributions.h: distributions
 * BRIAN Software Package Version 3.0
 *
 * $Id: distributions.h 493 2017-02-22 22:51:12Z kruggel $
 *
 * 0.30 (29/03/13): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.41 (09/08/14): mods to Watson distribution
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Classes for probability distributions.
*/

//! Implements a Mersenne twister pseudo-random number generator.
/*!
Implements the algorithm in:\n
Matsumoto & Nishimura (1988) ACM T Mod Comp Sim 8, pp 3-30.
*/ 

class mtrand {
	static const unsigned int n = 624;	//!< size of state array
	static const unsigned int m = 397;	//!< shift m
	static unsigned long s[n];		//!< state array
	static unsigned int p;			//!< position in s array

	unsigned long twiddle(const unsigned long u, const unsigned long v) const	//! twist u and v.
		{ return (((u & 0x80000000ul) | (v & 0x7FFFFFFFul)) >> 1) 
				^ (v & 1ul? 0x9908B0DFul: 0x0ul); }
	void	generate()								//! generate new state value.
		{ for (unsigned int i = 0; i < n-m; i++)
			s[i] = s[i+m] ^ twiddle(s[i], s[i+1]);
		  for (unsigned int i = n-m; i < n-1; i++)
			s[i] = s[i+m-n] ^ twiddle(s[i], s[i+1]);
		  s[n-1] = s[m-1] ^ twiddle(s[n-1], s[0]); p = 0; }
	void	seed()									//! initialize RNG from current time.
		{ s[0] = static_cast<unsigned long>(time(0l)) & 0xFFFFFFFFul;
		  for (unsigned int i = 1; i < n; i++) { 
			s[i] = 1812433253ul*(s[i-1]^(s[i-1] >> 30))+i; s[i] &= 0xFFFFFFFFul; }
		  p = n; }
public:
	mtrand()									//! allocates and initializes a RNG.
		{ seed(); }
	unsigned long randint()								//! draw random integer.
		{ if (p == n) generate();						// new s vector needed
		  unsigned long x = s[p++]; x ^= (x >> 11);
		  x ^= (x << 7) & 0x9D2C5680ul; x ^= (x << 15) & 0xEFC60000ul;
		  return x ^ (x >> 18); }
};

//! Implements a uniform distribution

template<typename T> class uniformDist {
	const T	min;				//!< lower bound
	const T	max;				//!< upper bound
public:
	uniformDist<T>(const T _min = T(0), const T _max = T(1))			//! constructs a uniform distribution in [min,max[.
		: min(std::min(_min,_max)), max(std::max(_min,_max))
		{ if (min == max) throw optException("uniformDist: Invalid parameter"); }
	T	operator()() const;							//! returns a sample of a uniform distribution.
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ T p = T(1)/(max-min);
		  return (x >= min && x <= max)? T(std::log(p)): std::numeric_limits<T>::lowest(); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ T p = T(1)/(max-min);
		  return (x >= min && x <= max)? p: T(0); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ if (x < min) return T(0); else if (x >= max) return T(1);
		  return (x-min)/(max-min); }
	T	mean() const								//! returns the mean of a uniform distribution.
		{ return T(0.5)*(min+max); }
	T	median() const								//! returns the median of a uniform distribution.
		{ return mean(); }
	T	var() const								//! returns the variance of a uniform distribution.
		{ return T(1.0/12.0)*SQR(max-min); }
	T	entropy() const								//! returns the entropy of a uniform distribution.
		{ return std::log(max-min); }
	T	quantile(const T p) const						//! returns the quantile at p of a uniform distribution.
		{ if (p < T(0) || p > T(1)) throw optException("uniformDist: Invalid parameter");
		  return min+p*(max-min); }
};

//! Implements a normal distribution

template<typename T> class normalDist {
	const T	mu;				//!< mean
	const T	sigma2;				//!< variance
	const T	c;				//!< normalization factor
public:
	normalDist<T>(const T _mu = T(0), const T _sigma2 = T(1))			//! constructs a normal distribution N(&mu;,&sigma;<SUP>2</SUP>).
		: mu(_mu), sigma2(_sigma2),
		c(T(std::log(std::sqrt(2.0*M_PI*sigma2))))
		{ if (sigma2 <= T(0)) throw optException("normalDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a normal distribution.
		{ return T(mu+sample()*std::sqrt(sigma2)); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ T t = T(0.5)*SQR(x-mu)/sigma2; return -t-c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return T(0.5)*(T(1)+std::erf((x-mu)/std::sqrt(T(2)*sigma2))); }
	T	mean() const								//! returns the mean of a normal distribution.
		{ return mu; }
	T	median() const								//! returns the median of a normal distribution.
		{ return mean(); }
	T	var() const								//! returns the variance of a normal distribution.
		{ return sigma2; }
	T	entropy() const								//! returns the entropy of a normal distribution.
		{ return T(0.5)*std::log(T(2)*M_PI*M_E*sigma2); }
	T	quantile(const T p) const						//! returns the quantile at p of a normal distribution.
		{ if (p < T(0) || p > T(1)) throw optException("normalDist: Invalid parameter");
		  return mu+T(std::sqrt(sigma2)*quant(p)); }
	static T sample();
	static T quant(const T p);
};

//! Implements a univariate std::exponential distribution

template<typename T> class expDist {
	const T	lambda;				//!< shape factor
	const T	c;				//!< normalization factor
public:
	expDist<T>(const T _lambda = T(1))						//! constructs an std::exponential distribution E(&lambda;).
		: lambda(_lambda), c(std::log(lambda))
		{ if (lambda <= 0) throw optException("std::exponentialDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of an std::exponential distribution.
		{ return sexponential()/lambda; }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  return -lambda*x+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x < T(0)? T(0): std::exp(lpdf(x)); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): T(1)-std::exp(-lambda*x); }
	T	mean() const								//! returns the mean of an std::exponential distribution.
		{ return T(1)/lambda; }
	T	median() const								//! returns the median of an std::exponential distribution.
		{ return T(std::log(T(2)))/lambda; }
	T	var() const								//! returns the variance of an std::exponential distribution.
		{ return T(1)/SQR(lambda); }
	T	entropy() const								//! returns the entropy of an std::exponential distribution.
		{ return T(1)-T(std::log(lambda)); }
	T	quantile(const T p) const						//! returns the quantile at p of an std::exponential distribution.
		{ if (p < T(0) || p > T(1)) throw optException("std::exponentialDist: Invalid parameter");
		  return -std::log(T(1)-p)/lambda; }
	static T sexponential();
};

//! Implements a Poisson distribution

template<typename T> class PoissonDist {
	const T	lambda;				//!< shape factor
public:
	PoissonDist<T>(const T _lambda = T(1))					//! constructs a Poisson distribution P(&lambda;).
		: lambda(_lambda)
		{ if (lambda <= T(0)) throw optException("PoissonDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a Poisson distribution.
		{ return sample(lambda); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ T xf = T(std::floor(x+T(0.5)));
		  if (xf < 0) return std::numeric_limits<T>::lowest();
		  else if (xf == 0) return -lambda;
		  else return T(-stirlerr(xf)-bd0(xf,lambda)-std::log(std::sqrt(2*M_PI*xf))); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x < T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): T(1.0-pGamma(x+1.0,lambda)); }
	T	mean() const								//! returns the mean of a Poisson distribution.
		{ return lambda; }
	T	median() const								//! returns the median of a Poisson distribution.
		{ return T(lambda+1.0/3.0-0.02/lambda); }				// this is an approximation
	T	var() const								//! returns the variance of a Poisson distribution.
		{ return lambda; }
	T	quantile(const T p) const						//! returns the quantile at p of a Poisson distribution.
		{ if (p < T(0) || p > T(1)) throw optException("PoissonDist: Invalid parameter");
		  return quant(p, lambda); }
	static T sample(const T lambda);
	static T quant(const T p, const T lambda);
};

//! Implements a gamma distribution

template<typename T> class gammaDist {
	const T	alpha;				//!< shape factor
	const T	beta;				//!< shape factor
	const T	c;				//!< normalization factor
public:
	gammaDist<T>(const T _alpha, const T _beta)					//! constructs a gamma distribution gamma(&alpha;,&beta;).
		: alpha(_alpha), beta(_beta),
		c(alpha*T(std::log(beta))-T(std::lgamma(alpha)))
		{ if (alpha < T(0) || beta <= T(0)) throw optException("gammaDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a gamma distribution.
		{ return sample(alpha)/beta; }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  T t = -x*beta+(alpha-T(1))*T(std::log(x)); return t+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x < T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): pGamma(alpha,beta*x); }
	T	mean() const								//! returns the mean of a gamma distribution.
		{ return alpha/beta; }
	T	median() const								//! returns the median of a gamma distribution.
		{ return mean()*(T(3)*alpha-T(0.8))/(T(3)*alpha+T(0.2)); }		// this is an approximation
	T	var() const								//! returns the variance of a gamma distribution.
		{ return alpha/SQR(beta); }
	T	entropy() const								//! returns the entropy of a gamma distribution.
		{ return alpha-std::log(beta)+std::lgamma(alpha)+(T(1)-alpha)*digamma(alpha); }
	T	quantile(const T p) const						//! returns the quantile at p of a gamma distribution.
		{ if (p < T(0) || p > T(1)) throw optException("gammaDist: Invalid parameter");
		  return quant(p,alpha,beta); }
	static T sample(const T alpha);
	static T quant(const T p, const T alpha, const T beta);
};

//! Implements an inverse gamma distribution

template<typename T> class invGammaDist {
	const T	alpha;				//!< shape factor
	const T	beta;				//!< shape factor
	const T	c;				//!< normalization factor
public:
	invGammaDist<T>(const T _alpha, const T _beta)					//! constructs an inverse gamma distribution invGamma(&alpha;,&beta;).
		: alpha(_alpha), beta(_beta),
		c(alpha*std::log(beta)-std::lgamma(alpha))
		{ if (alpha < T(0) || beta <= T(0)) throw optException("invGammaDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of an inverse gamma distribution.
		{ return beta/gammaDist<T>::sample(alpha); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  T t = -beta/x+(-alpha-T(1))*T(std::log(x)); return t+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x <= T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): qGamma(alpha,beta/x); }
	T	mean() const								//! returns the mean of an inverse gamma distribution.
		{ if (alpha <= T(1)) throw optException("invGammaDist: Invalid parameter");
		  return beta/(alpha-T(1)); }
	T	var() const								//! returns the variance of an inverse gamma distribution.
		{ if (alpha <= T(2)) throw optException("invGammaDist: Invalid parameter");
		  return SQR(beta)/(SQR(alpha-T(1))*(alpha-T(2))); }
	T	entropy() const								//! returns the entropy of an inverse gamma distribution.
		{ return alpha+std::log(beta)+std::lgamma(alpha)-(T(1)+alpha)*digamma(alpha); }
};

//! Implements a central &chi;<SUP>2</SUP> distribution

template<typename T> class chiSqDist : public gammaDist<T> {
public:
	chiSqDist<T>(const unsigned int k, const T _sigma2 = T(1))			//! defined from gamma distribution gamma(k/2,1/(2*&sigma;<SUP>2</SUP>)).
		 : gammaDist<T>(T(0.5)*k,T(0.5/_sigma2)) { }
};

//! Implements a non-central &chi;<SUP>2</SUP> distribution

template<typename T> class nChiSqDist {
	const unsigned int k;			//!< DOF
	const T	mu;				//!< mean
	const T	sigma2;				//!< variance
	const T	lambda;				//!< shape factor
public:
	nChiSqDist<T>(const unsigned int _k, const T _mu, const T _sigma2 = T(1))	//! constructs a non-central &chi;<SUP>2</SUP> distribution.
		: k(_k), mu(_mu), sigma2(_sigma2), lambda(k*SQR(mu)) { }
	T	operator()() const							//! returns a sample of a non-central &chi;<SUP>2</SUP> distribution.
		{ T s = 0; for (unsigned int i = 0; i < k; i++) {
			T v = mu+normalDist<T>::sample()*T(std::sqrt(sigma2)); s += SQR(v)/sigma2; }
		  return s; }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= 0) return std::numeric_limits<T>::lowest();
		  T t1 = T(0.25*k-0.5)*std::log(x/lambda)-T(0.5)*SQR(std::sqrt(x)-std::sqrt(lambda))/sigma2;
		  T t2 = std::log(besseli_sc(T(0.5*k-1.0), std::sqrt(x*lambda)/sigma2));
		  return t1+t2+std::log(T(0.5)/sigma2); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x <= T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): T(1)-MarcumQ(T(0.5*k),std::sqrt(lambda),std::sqrt(x)); }
	T	mean() const								//! returns the mean of a non-central &chi;<SUP>2</SUP> distribution.
		{ return k*sigma2+lambda; }
	T	var() const								//! returns the variance of a non-central &chi;<SUP>2</SUP> distribution.
		{ return 2*sigma2*(k+2*lambda); }
};

//! Implements a central &chi; distribution

template<typename T> class chiDist {
	const unsigned int k;			//!< DOF
	const T	sigma2;				//!< variance
	const T	c;				//!< normalization factor
public:
	chiDist<T>(const unsigned int _k, const T _sigma2 = T(1))			//! constructs a central &chi; distribution.
		: k(_k), sigma2(_sigma2), c(T((1.0-0.5*k)*M_LN2-std::lgamma(0.5*k)-0.5*std::log(sigma2)))
		{ if (k < 2) throw optException("chiDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a central &chi; distribution.
		{ T s = 0; for (unsigned int i = 0; i < k; i++) {
			T v = normalDist<T>::sample(); s += SQR(v); };
		  return std::sqrt(s*sigma2); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= 0) return std::numeric_limits<T>::lowest();
		  T x2 = SQR(x)/sigma2, t = T(-0.5*x2+0.5*(k-1)*std::log(x2)); return t+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x <= T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): pGamma(T(0.5)*k,T(0.5)*SQR(x)/sigma2); }
	T	mean() const								//! returns the mean of a central &chi; distribution.
		{ const T k2 = T(0.5*k), t = std::lgamma(k2+T(0.5))-std::lgamma(k2);	// checked 23/06/16 FK
		  return std::sqrt(T(2)*sigma2)*std::exp(t); }
	T	var() const								//! returns the variance of a central &chi; distribution.
		{ const T k2 = T(0.5*k), t = std::lgamma(k2+T(1))-std::lgamma(k2);	// checked 23/06/16 FK
		  const T m2 = T(2)*sigma2*std::exp(t), mn = mean(); return m2-SQR(mn); }
};

//! Implements a non-central &chi; distribution

template<typename T> class nChiDist {							// see Anderson (1981) Indian J Stat 48, 58-67
	const unsigned int k;			//!< DOF
	const T	mu;				//!< mean
	const T	sigma2;				//!< variance
	const T	c;				//!< normalization factor
public:
	nChiDist<T>(const unsigned int _k, const T _mu, const T _sigma2 = T(1))	//! constructs a non-central &chi; distribution.
		: k(_k), mu(_mu), sigma2(_sigma2), c(std::log(mu/sigma2))
		{ if (mu <= T(0) || sigma2 <= 0) throw optException("nChiDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a non-central &chi; distribution.
		{ T s = 0; for (unsigned int i = 0; i < k; i++) {
			T v = mu+normalDist<T>::sample()*T(std::sqrt(sigma2)); s += SQR(v); }
		  return std::sqrt(s); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= 0) return std::numeric_limits<T>::lowest();
		  T t1 = T(0.5*k*std::log(x/mu)-0.5*SQR(x-mu)/sigma2);
		  T t2 = T(std::log(besseli_sc(0.5*k-1.0, x*mu/sigma2)));
		  return t1+t2+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x <= T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ return x < T(0)? T(0): T(1.0-MarcumQ(T(0.5)*k,mu/std::sqrt(sigma2),x/std::sqrt(sigma2))); }
	T	mean() const								//! returns the mean of a non-central &chi; distribution.
		{ const T k2 = T(0.5*k), t = T(std::lgamma(k2+0.5)-std::lgamma(k2));	// checked 23/06/16 FK
		  const T s2 = T(2.0*sigma2), x2 = T(k*SQR(mu)/s2);
		  return T(std::sqrt(s2)*std::exp(t)*Hyper1F1(-x2,-0.5,k2)); }
	T	var() const								//! returns the variance of a non-central &chi; distribution.
//		{ const T m = mean(); return mu+sigma2*k; }
		{ const T k2 = T(0.5*k), t = std::lgamma(k2+T(1))-std::lgamma(k2);	// checked 23/06/16 FK
		  const T m2 = T(2)*sigma2*std::exp(t), mn = mean();
		  return m2-SQR(mn)+k*SQR(mu); }
};

//! Implements a Dirichlet distribution

template<typename T> class DirichletDist {
	const vecD<T> alpha;			//!< shape parameter vector
	const T	c;				//!< normalization factor

	vecD<T>	init(unsigned int n, const T al)					//! initializes shape vector of length n with &alpha;.
		{ vecD<T> v(n); for (unsigned int i = 0; i < n; i++) v[i] = T(al/n);
		  return v; }
	T	norm()									//! returns normalization factor.
		{ T s{0}; for (unsigned int i = 0; i < alpha.N; i++)
				s += T(std::lgamma(alpha(i)));
		  return s-T(std::lgamma(alpha.sum())); }
public:
	DirichletDist<T>(const vecD<T>& _alpha)						//! constructs a Dirichlet distribution from vector &alpha;.
		: alpha(_alpha), c(norm())
		{ if (alpha.N < 2) throw optException("DirichletDist: Invalid parameter"); }
	DirichletDist<T>(const unsigned int n, T al)					//! constructs a Dirichlet distribution of dimension n and common &alpha;.
		: alpha(init(n,al)), c(norm()) { }
	vecD<T> operator()() const							//! returns a sample of a Dirichlet distribution.
		{ unsigned int n = alpha.N; vecD<T> a(n);
		  for (unsigned int i = 0; i < n; i++) a[i] = gammaDist<T>::sample(alpha(i));
		  return a/a.sum(); }
	T	lpdf(const vecD<T>& x) const						//! returns log PDF at x.
		{ T xp = 0; assert(x.N == alpha.N);
		  for (unsigned int i = 0; i < x.N; i++) {
			if (x(i) <= T(0)) return std::numeric_limits<T>::lowest();
			xp += std::log(x(i))*(alpha(i)-T(1)); };
		  return xp-c; }
	T	pdf(const vecD<T>& x) const						//! returns PDF at x.
		{ assert(x.N == alpha.N);
		  for (unsigned int i = 0; i < x.N; i++) { if (x(i) <= T(0)) return 0; };
		  return T(std::exp(lpdf(x))); }
	vecD<T> mean() const								//! returns the mean of a Dirichlet distribution.
		{ return alpha/alpha.sum(); }
	vecD<T> var() const								//! returns the variance of a Dirichlet distribution.
		{ T s = alpha.sum(); unsigned int n = alpha.N; vecD<T> a(n);
		  for (unsigned int i = 0; i < n; i++) a[i] = alpha(i)*(s-alpha(i));
		  T d = SQR(s)*(s+T(1)); return a/d; }
};

//! Implements a beta distribution

template<typename T> class betaDist {
	const T alpha;				//!< shape factor
	const T beta;				//!< shape factor
	const T	c;				//!< normalization factor
public:
	betaDist<T>(const T _alpha, const T _beta)					//! constructs a beta distribution.
		: alpha(_alpha), beta(_beta), c(lBeta(alpha,beta)) { }
	T	operator()() const							//! returns a sample of a beta distribution.
		{ return sample(alpha,beta); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  if (x >= T(1)) return std::numeric_limits<T>::lowest();
		  return (alpha-T(1))*std::log(x)+(beta-T(1))*std::log(T(1.0-x))-c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ if (x <= T(0)) return T(0); if (x >= T(1)) return T(1);
		  return iBeta(alpha,beta,x); }
	T	mean() const								//! returns the mean of a beta distribution.
		{ return alpha/(alpha+beta); }
	T	var() const								//! returns the variance of a beta distribution.
		{ return alpha*beta/(SQR(alpha+beta)*(alpha+beta+T(1))); }
	T	entropy() const								//! returns the entropy of a beta distribution.
		{ T a = (alpha-T(1))*digamma(alpha), b = (beta-T(1))*digamma(beta);
		  T s = (alpha+beta-T(2))*digamma(alpha+beta); return c-a-b+s; }
	T	quantile(const T p) const						//! returns the quantile of a beta distribution.
		{ if (p < T(0) || p > T(1)) throw optException("tDist: Invalid parameter");
		  return quant(p, alpha, beta); }
	static T sample(const T alpha, const T beta);
	static T quant(const T p, const T alpha, const T beta);
};

//! Implements a binomial distribution

template<typename T> class binomialDist {
	const unsigned int n;			//!< number of draws
	const T p;				//!< probability

	double	search(T y, T& z, const T pr, const T d) const				//! helper for quantile function.
		{ if (z >= pr) { while (1) { T newz = cdf(y-d);
			if (y == T(0) || newz  < pr) return y;
			y = std::max(T(0), y-d); z = newz; } }
		  else { while (1) { y = std::min(y+d, n); z = cdf(y);
			if (y == n || z  >= pr) return y; } } }
public:
	binomialDist<T>(const unsigned int _n, const T _p)				//! constructs a binomial distribution.
		: n(_n), p(_p)
		{ if (p < T(0) || p > T(1)) throw optException("binomialDist: Invalid parameter"); }
	int	operator()() const							//! returns a sample of a binomial distribution.
		{ return sample(n,p); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ return lpBinomial(std::floor(x+0.5),n,p,1.0-p); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ if (x < T(0)) return T(0); if (x >= n) return T(1);
		  T xf = T(std::floor(x+1e-6)), n1 = T(std::floor(n+0.5));
		  return iBeta(n1-xf,1.0+xf,1.0-p); }
	T	mean() const								//! returns the mean of a binomial distribution.
		{ return n*p; }
	T	var() const								//! returns the variance of a binomial distribution.
		{ return n*p*(T(1)-p); }
	T	entropy() const								//! returns the entropy of a binomial distribution.
		{ T t = T(std::log(2.0*M_PI*M_E*var())); return T(0.5)*t/M_LN2; }
	T	quantile(T x) const							//! returns the quantile of a binomial distribution.
		{ if (x < T(0) || x > T(1)) throw optException("binomialDist: Invalid parameter");
		  if (p == T(0) || n == 0) return T(0);
		  T q = T(1)-p; if (q == T(0)) return n;
		  T mu = n*p, sigma = std::sqrt(n*p*q), gamma = (q-p)/sigma;
		  if (x+T(1.01)*std::numeric_limits<T>::epsilon() >= T(1)) return n;
		  T z = normalDist<T>().quantile(x);
		  T y = T(std::floor(mu+sigma*(z+gamma*(z*z-1.0)/6.0)+0.5)); if (y > n) y = n;
		  z = cdf(y); x *= T(1)-T(64)*std::numeric_limits<T>::epsilon();
		  if (n < 1e5) return search(y,z,x,1);
		  T d = T(std::floor(n*0.001)), o;
		  do { o = d; y = search(y,z,x,d); d = T(std::max(1.0,std::floor(0.01*d))); }
		  while (o > 1 && d > T(n*1e-6));
		  return y; }
	static unsigned int sample(unsigned int n, T p);
	static T lpBinomial(const T x, const T n, const T p, const T q);
};

//! Implements a multinomial distribution

template<typename T> class multinomialDist {
	const unsigned int n;			//!< number of draws
	const vecD<T>& p;			//!< probability vector
public:
	multinomialDist<T>(const unsigned int _n, const vecD<T>& _p)			//! constructs a multinomial distribution.
		: n(_n), p(_p)
		{ if (p.sum() != T(1)) throw optException("binomialDist: Invalid parameter"); }
	uvecD	operator()() const							//! returns a sample of a multinomial distribution.
		{ return sample(n,p); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  if (x >= n) return std::numeric_limits<T>::lowest();
		  return binomialDist<T>::lpBinomial(std::floor(x+0.5),n,p,1.0-p); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	vecD<T>	mean() const								//! returns the mean of a multinomial distribution.
		{ unsigned int N = p.V; vecD<T> mu(N);
		  for (unsigned int i = 0; i < N; i++) mu[i] = n*p(i);
		  return mu; }
	T	var() const								//! returns the variance of a multinomial distribution.
		{ unsigned int N = p.V; vecD<T> v(N);
		  for (unsigned int i = 0; i < N; i++) v[i] = n*p(i)*(T(1)-p(i));
		  return v; }
	static uvecD sample(const unsigned int n, const vecD<T>& p);
};

//! Implements a multivariate normal distribution

template<typename T> class mvNormalDist {
	const vecD<T> mu;			//!< mean vector
	const matD<T> Sigma;			//!< covariance matrix
	const unsigned int N;			//!< dimensionality
	const matD<T> S;			//!< for sample generation

	matD<T>	initS() const
		{ matD<T> D(N,N),V; vecD<T> v = sev(Sigma,V); D = T(0);
		  for (unsigned int i = 0; i < N; i++) D(i,i) = T(std::sqrt(v(i)));
		  return V*D; }
public:
	mvNormalDist<T>(const vecD<T>& _mu, const matD<T>& _Sigma) 			//! constructs a multivariate normal distribution.
		: mu(_mu), Sigma(_Sigma), N(mu.N), S(initS())
		{ if (Sigma.M != N || Sigma.N != N)
			throw optException("mvNormalDist: Incompatible arguments"); }
	T	operator()() const							//! returns a sample of a multivariate normal distribution.
		{ vecD<T> rn(N); for (unsigned int i = 0; i < N; i++)
			rn[i] = normalDist<T>::sample();
		  return mu+S*rn; }
	T	lpdf(const vecD<T>& x) const						//! returns log PDF at x.
		{ if (x.N != N) throw optException("mvNormalDist: Incompatible arguments");
		  const vecD<T> d = x-mu; const T t = T(0.5)*dot(d,inv(Sigma)*d);
		  const T c = T(0.5)*(N*std::log(2.0*M_PI)+std::log(det(Sigma)));
		  return -t-c; }
	T	pdf(const vecD<T>& x) const						//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	vecD<T> mean() const								//! returns the mean of a multivariate normal distribution.
		{ return mu; }
	matD<T> var() const								//! returns the variance of a multivariate normal distribution.
		{ return Sigma; }
	T	entropy() const								//! returns the entropy of a multivariate normal distribution.
		{ const T t = std::log(det(Sigma));
		  return T(0.5)*(N+N*std::log(T(2)*M_PI)+t); }
};

//! Implements a Wishart distribution

template<typename T> class WishartDist {
	const matD<T>& V;			//!< scale matrix
	const unsigned int n;			//!< DOF
public:
	WishartDist<T>(const matD<T>& _V, const unsigned int _n) 			//! constructs a Wishart distribution.
		: V(_V), n(_n) { }
	T	operator()() const							//! returns a sample of a Wishart distribution.
		{ const unsigned int N = V.N; if (n < N-1)				// see Smith & Hocking, Applied Statistics, 21, 341
			throw optException("WishartDist: Argument error");
		  matD<T> L = CholDec(V), A(N,N); A = 0;
		  for (unsigned int j = 0; j < N; j++) {
			A(j,j) = chiDist<T>(n-j+1)();
			for (unsigned int i = j+1; i < N; i++)
				A(j,i) = normalDist<T>::sample(); };
		  matD<T> LA = L*A; return LA*trp(LA); }
	T	lpdf(const matD<T>& X) const						//! returns log PDF at X(p,p).
		{ const unsigned int p = V.N;						// see Wikipedia
		  if (p != X.N)	throw optException("WishartDist: Argument error");
		  const T n1 = T(0.5)*(n-p-1)*std::log(det(X));				// nominator, part 1
		  const T n2 = -T(0.5)*(inv(V)*X).trace();				// nominator, part 2
		  const T d1 = T(0.5)*n*p*M_LN2+T(0.5)*n*std::log(det(V));		// denominator, part 1
		  const T d2 = lmvGamma(p,T(0.5)*n);					// denominator, part 2 - uses multivariate gamma function
		  return n1+n2-d1-d2; }
	T	pdf(const vecD<T>& x) const						//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	matD<T> mean() const								//! returns the mean of a Wishart distribution.
		{ return n*V; }
	matD<T> var() const								//! returns the variance of a Wishart distribution.
		{ unsigned int p = V.N; matD<T> v(p,p);
		  for (unsigned int i = 0; i < p; i++)
		  	for (unsigned int j = 0; j < p; j++) 
				v(i,j) = n*(SQR(V(i,j))+V(i,i)*V(j,j));
		  return v; }
	T	entropy() const								//! returns the entropy of a Wishart distribution.
		{ unsigned int p = V.N; T t1 = T(0.5)*(p+1)*(std::log(det(V))+p*M_LN2);
		  T t2 = lmvGamma(p,T(0.5)*n)+T(0.5)*n*p, t3 = T(0);
		  for (unsigned int i = 0; i < p; i++) t3 += digamma(T(0.5)*(n-i));
		  return t1+t2-T(0.5)*(n-p-1)*t3; }
};

//! Implements a Watson distribution on the sphere

template<typename T> class WatsonDist {
	const vec3<T>& mu;			//!< mean vector
	const T	k;				//!< concentration
public:
	WatsonDist<T>(const vec3<T>& _mu, const T _k) 					//! constructs a Watson distribution.
		: mu(_mu), k(_k) { }
	T	operator()() const							//! returns a sample of a Watson distribution.
		{ // see Best & Fisher (1986) Austral J Stat 28, 13-31.
		  if (k <= 0) throw optException("WatsonDist: Argument error");
		  uniformDist<T> ud;
		  T c = k/(std::exp(k)-T(1)), u0 = ud(), u1 = ud();
		  T y = std::log(u0*k/c + T(1))/k, z = y(u1 < std::exp(k*(y*y-y)));
		  T r = std::sqrt(1-z*z), phi = 2*M_PI*ud();
		  vecD<T> X(3); X[0] = r*cos(phi); X[1] = r*sin(phi); X[2] = z;
		  matD<T> W(3,3), R(3,3); W = 0; R.id();
		  W(0,2) = mu(0); W(1,2) = mu(1); W(2,0) = -mu(0); W(2,1) = -mu(1); 
		  W /= std::sqrt(SQR(mu(0))+SQR(mu(1)));
		  R += std::sqrt(T(1)-SQR(mu(2)))*W + (1.0-mu(2))*(W*W);
		  return R*X; }
	T	lpdf(const vec3<T>& x) const						//! returns log PDF at x.
		{ T t = dot(mu,x); return k*SQR(t)-logHyper1F1(k,0.5,1.5)-std::log(4*M_PI); }
	T	pdf(const vec3<T>& x) const						//! returns PDF at x.
		{ T t = dot(mu,x); return T(std::exp(k*SQR(t))/(4*M_PI*Hyper1F1(k,0.5,1.5))); }
};

//! Implements an F distribution

template<typename T> class FDist {
	const T d1;				//!< DOF1
	const T d2;				//!< DOF2
	const T c;				//!< normalization factor
public:
	FDist<T>(const T _d1, const T _d2) 						//! constructs an F distribution F(d1,d2).
		: d1(_d1), d2(_d2), c(T(0.5*d1*std::log(d1/d2)-lBeta(0.5*d1,0.5*d2)))
		{ if (d1 <= T(0) || d2 <= T(1))
			throw optException("FDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of an F distribution.
		{ const T s1 = gammaDist<T>::sample(T(0.5)*d1)/d1;
		  const T s2 = gammaDist<T>::sample(T(0.5)*d2)/d2; return s1/s2; }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= T(0)) return std::numeric_limits<T>::lowest();
		  T t1 = T((0.5*d1-1.0)*std::log(x));
		  T t2 = T((0.5*(d1+d2))*std::log(1.0+d1*x/d2));
		  return t1-t2+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return x <= T(0)? T(0): T(std::exp(lpdf(x))); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ if (x <= T(0)) return T(0);
		  return T(iBeta(0.5*d1,0.5*d2,d1*x/(d1*x+d2))); }
	T	mean() const								//! returns the mean of an F distribution.
		{ if (d2 <= T(2)) throw optException("FDist: Invalid parameter");
		  return d2/(d2-T(2)); }
	T	var() const								//! returns the variance of an F distribution.
		{ if (d2 <= T(4)) throw optException("FDist: Invalid parameter");
		  const T n = T(2)*SQR(d2)*(d1+d2-T(2));
		  const T d = d1*SQR(d2-T(2.0))*(d2-T(4)); return n/d; }
	T	quantile(const T p) const						//! returns the quantile of an F distribution.
		{ if (p < T(0) || p > T(1)) throw optException("FDist: Invalid parameter");
		  if (d1 <= d2 && d2 > 4e5) {
			chiSqDist<T> cd(d1); return cd.quantile(p)/d1; }
		  else if (d1 > 4e5) {
			chiSqDist<T> cd(d2); return d2/(T(1)-cd.quantile(p)); }
		  else return (T(1)/(T(1)-qBeta(p,T(0.5*d1),T(0.5*d2)))-T(1))*(d2/d1); }
};

//! Implements a Student t distribution

template<typename T> class tDist {
	const T df;				//!< DOF
	const T c;				//!< normalization factor
public:
	tDist<T>(const T _df) 								//! constructs a t distribution t(df).
		: df(_df),c(std::lgamma(T(0.5)*(df+T(1)))-std::lgamma(T(0.5)*df)
				-std::log(std::sqrt(df*M_PI)))
		{ if (df <= T(0)) throw optException("tDist: Invalid parameter"); }
	T	operator()() const							//! returns a sample of a t distribution.
		{ T z = normalDist<T>::sample(), v = T(2)*gammaDist<T>::sample(T(0.5)*df);
		  return z*std::sqrt(df/v); }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ T t = (T(0.5)*(df+T(1)))*std::log(T(1)+SQR(x)/df); return -t+c; }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return std::exp(lpdf(x)); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ const double nx = 1.0+(double(x)/double(df))*double(x);		// leave this as double
		  const double v = nx < 1e100? iBeta(0.5*df,0.5,1.0/nx):		// because iBeta has double arguments
		  	-0.5*df*(2.0*std::log(std::abs(x))-std::log(df))-lBeta(0.5*df,0.5)-std::log(0.5*df);
		  return x <= T(0)? T(0.5*v): T(1.0-0.5*v); }
	T	mean() const								//! returns the mean of a t distribution.
		{ if (df <= T(1)) throw optException("tDist: Invalid parameter");
		  return T(0); }
	T	var() const								//! returns the variance of a t distribution.
		{ if (df <= T(2)) throw optException("tDist: Invalid parameter");
		  return df/(df-T(2)); }
	T	entropy() const								//! returns the entropy of a t distribution.
		{ T v = T(0.5)*(T(1)+df), t1 = digamma(v)-digamma(0.5*df);
		  T t2 = lBeta(0.5*df,0.5)+std::log(T(0.5)*df);
		  return v*t1+t2; }
	T	quantile(const T p) const						//! returns the quantile of a t distribution.
		{ if (p < T(0) || p > T(1)) throw optException("tDist: Invalid parameter");
		  return quant(p,df); }
	static T quant(const T p, const T df);
};

//! Implements an empirical distribution

template<typename T> class empiricalDist {
protected:
	const T* v;				//!< samples
	const unsigned int N;			//!< number of samples
	T	min;				//!< lowest sample
	T	max;				//!< highest sample
	vecD<T>	hs;				//!< histogram

	T	bw() const								//! returns bandwidth of sample bin.
		{ return almostEqual(min,max)? T(0): (max-min)/(hs.N-1); }
	unsigned int bin(const T x) const						//! computes bin for observation x.
		{ return almostEqual(min,max)? 0: FTOU((x-min)/bw()); }
	void	setBounds()								//! determines bounds of sample.
		{ min = v[0]; max = v[0];
		  for (unsigned int i = 1; i < N; i++) {
			min = std::min(min,v[i]); max = std::max(max,v[i]); } }
	void	hist()									//! computes the empirical PDF.
		{ hs = 0; if (almostEqual(min,max)) { hs[0] = T(1); return; };
		  for (unsigned int i = 0; i < N; i++) {
			unsigned int a = bin(v[i]); if (a < hs.N) hs[a] += T(1); };
		  T s = hs.sum(); hs /= s; }
	void	init()									//! determines bounds and PDF.
		{ if (N < 1) throw optException("empiricalDist: Invalid parameter"); 
		  setBounds(); hist(); }
	empiricalDist<T>(const empiricalDist<T>& e) = delete;
	empiricalDist<T> operator=(const empiricalDist<T>& e) = delete;
public:
	empiricalDist<T>(const vecD<T>& x, const unsigned int B = 1000)			//! constructs an empirical distribution from a vector.
		: v(x.x), N(x.N), min(0), max(0), hs(std::min(B,N)) { init(); }
	empiricalDist<T>(const image<T>& f, const unsigned int B = 1000)		//! constructs an empirical distribution from an image.
		: v(f.x), N(f.nel()), min(0), max(0), hs(std::min(B,N)) { init(); }
	T	operator()() const							//! returns a sample of the observations.
		{ uniformDist<unsigned int> ud; return v[FTOU(ud()*N)]; }
	T	lpdf(const T x) const							//! returns log PDF at x.
		{ if (x <= min || x >= max) return std::numeric_limits<T>::lowest();
		  T p = hs(bin(x)); return p > 0? T(std::log(p)): std::numeric_limits<T>::lowest(); }
	T	pdf(const T x) const							//! returns PDF at x.
		{ return (x <= min || x >= max)? T(0): hs(bin(x)); }
	T	cdf(const T x) const							//! returns CDF at x.
		{ if (x <= min) return T(0); else if (x >= max) return T(1);
		  int b = bin(x); T s = 0; for (int i = 0; i <= b; i++) s += hs(i);
		  T xl = std::max(min,min+(b-1)*bw()), fr = (x-xl)/bw(); return s+fr*hs(b); }
	T	mean() const								//! returns the mean of the observations.
		{ T s = 0; for (unsigned int i = 0; i < N; i++) s += v[i]; return s/N; }
	T	var() const								//! returns the variance of the observations.
		{ if (N < 2) throw optException("empiricalDist: Invalid parameter");
		  T m = mean(), s = 0; for (unsigned int i = 0; i < N; i++) s += SQR(v[i]-m);
		  return s/(N-1); }
	T	median() const								//! returns the median of the observations.
		{ return quantile(T(0.5)); }
	T	mode() const								//! returns the mode of the observations.
		{ unsigned int m = hs.imax(); return T(m+0.5)*(max-min)/hs.N+min; }
	T	entropy() const								//! returns the entropy of the observations.
		{ T s = 0; for (unsigned int i = 0; i < hs.N; i++)
			if (hs(i)) s -= hs(i)*std::log(hs(i));
		  return s-std::log(hs.N); }
	T	quantile(const T p) const						//! returns the quantile of the observations.
		{ if (p < T(0) || p > T(1)) throw optException("empiricalDist: Invalid parameter");
		  if (isSmall<T>(p)) return min; else if (almostEqual(p,T(1))) return max; 
		  T s = 0; unsigned int i; for (i = 0; i < hs.N; i++) { s += hs(i); if (s >= p) break; };
		  s = std::max(s-hs(i),T(0)); const T x = std::max(min,min+int(i)*bw());
		  const T fr = (p-s)/hs(i); return x+fr*bw(); }
	T	dist(const empiricalDist<T>& b)
		{ assert(hs.N == b.hs.N); T d = T(0);
		  for (unsigned int i = 0; i < hs.N; i++) d += std::abs(hs(i)-b.hs(i));
		  return d; }
};

//! Implements a bounded empirical distribution

template<typename T> class boundedDist : public empiricalDist<T> {
	const T	lf;				//!< lower bound
	const T	hf;				//!< upper bound

	void	setBounds()								//! set bounds from quantiles.
		{ this->min = this->quantile(lf); this->max = this->quantile(hf); }
	void	init()									//! determines bounds and PDF.
		{ if (this->N < 1) throw optException("boundedDist: Invalid parameter"); 
		  setBounds(); this->hist(); }
public:
	boundedDist<T>(const vecD<T>& x, const T _lf, const T _hf)			//! constructs bounded distribution from a vector.
		 : empiricalDist<T>(x), lf(_lf), hf(_hf) { init(); }
	boundedDist<T>(const image<T>& f, const T _lf, const T _hf)			//! constructs bounded distribution from an image.
		 : empiricalDist<T>(f), lf(_lf), hf(_hf) { init(); }
};
#endif
