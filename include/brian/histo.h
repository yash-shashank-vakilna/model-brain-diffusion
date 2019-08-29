#ifndef HISTO_H
#define HISTO_H

/*
 *
 * histo.h: definitions for histogram-based classification
 * BRIAN Software Package Version 3.0
 *
 * $Id: histo.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (20/01/10): adapted to BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.41 (20/01/14): comparsion added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions for histogram-based classification methods.
 */

//! Collects distribution parameters for a Gaussian mixture.

struct parameter {
	float	mu;					//!< class mean
	float	sigma2;					//!< class variance
	float	pp;					//!< prior probability (fraction)

	parameter(const float _m = 0.0f, const float _s = 0.0f, const float _p = 0.0f)	//! allocates a parameter record with mean m, variance s, and fraction p.
		 : mu(_m), sigma2(_s), pp(_p) { }
	void	set(const float _mu, const float _sigma2, const float _p)		//! sets mean m, variance s, and fraction p.
		{ mu = _mu; sigma2 = _sigma2; pp = _p; }
	void	print(const unsigned int c) const					//! prints parameters preceeded by c.
		{ printf("class %u: %7.2f %7.2f %6.4f\n", c, mu, std::sqrt(sigma2), pp); }
	float	dist(const float v) const						//! returns normalized distance of v from mean.
		{ return 0.5f*SQR(v-mu)/sigma2; }
	float	logprob(const float v) const						//! returns -std::log PDF(v).
		{ return dist(v)+std::log(std::sqrt(float(2.0*M_PI*sigma2))); }
	float	prob(const float v) const						//! returns PDF(v).
		{ return std::exp(-logprob(v)); }
	bool operator< (const parameter &b) const					//! comparison by increasing mean.
		{ return (this->mu < b.mu); }
};

using vparam = std::vector<parameter>;

inline unsigned int findClass(const vparam& theta, const float v)			//! returns best class for v in parameter vector theta.
{	if (v <= theta[0].mu) return 0;	
	unsigned int nc = ITOU(theta.size());
	if (v >= theta[nc-1].mu) return nc-1;
	unsigned int cmin = 0; float dmin = std::numeric_limits<float>::max();
	for allClasses(c) {
		float d = SQR(v-theta[c].mu)/theta[c].sigma2;
		if (d < dmin) { dmin = d; cmin = c; } };
	return cmin;
}

#endif 
