#ifndef INTEGRATOR_H
#define INTEGRATOR_H

/*
 *
 * integrator.h: integration of a one-dimensional function
 * BRIAN Software Package Version 3.0
 *
 * $Id: integrator.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (20/07/14): first implementation
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Integration of a one-dimensional function.

    \details For detailed information, refer to:
	Philip Davis, Philip Rabinowitz, Methods of Numerical Integration,
	Second Edition, Dover, 2007.
*/

#define allPoints(i)	(unsigned int i = 0; i < n; i++)

//! Computes integral of a function from a to b over n intervals using the Clenshaw-Curtis algorithm.

template<typename T> class integrator {
	const unsigned int n;
	vecD<T>	x;
	vecD<T>	w;

	void	init()
		{ if (n == 1) { x[0] = 0.0; w[0] = 2.0; return; }
		  for allPoints(i) x[i] = cos((n-1-i)*M_PI/(n-1));
		  x[0] = -1.0; x[n-1] = 1.0; if (n%2 == 1) x[(n-1)/2] = 0.0;
		  for allPoints(i) { T t = i*M_PI/(n-1); w[i] = 1.0;
			for (unsigned int j = 1; j <= (n-1)/2; j++) {
				T b = 2*j == n-1? 1.0: 2.0;
				w[i] -= b*cos(2.0*j*t)/(4*j*j-1); } }
		  w[0] /= n-1; w[n-1] /= n-1;
		  for (unsigned int i = 1; i < n-1; i++) w[i] /= 0.5*(n-1); }
	void	rescale(const T a, const T b)
		{ for allPoints(i) x[i] = 0.5*((a+b)+(b-a)*x[i]);
		  for allPoints(i) w[i] *= 0.5*(b-a); }
	virtual T func(const T t) const = 0;
public:
	integrator(const unsigned int _n, const T a, const T b)
		: n(_n), x(n), w(n)
		{ if (n == 0 || b < a) throw optException("integrator argument error");
		  init(); rescale(a,b); }
	virtual ~integrator() { }
	T work()
		{ T s = 0; for allPoints(i) s += w(i)*func(x(i)); return s; }
};

#undef allPoints
#endif


