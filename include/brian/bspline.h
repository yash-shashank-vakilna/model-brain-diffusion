#ifndef BSPLINE_H
#define BSPLINE_H

/*
 *
 * bspline.h: b-splines in 3d
 * BRIAN Software Package Version 3.0
 *
 * $Id: bspline.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (16/05/10): initial version
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for b-splines in 3d.
*/

//! Represents a context for spline interpolation of images.
class bspline4  {
	const unsigned int degree = 4;
	fimage	coeff;					//!< coefficient image
	float	pole[2];				//!< poles
	float	lambda;					//!< smoothness

	float	weight4(const float x) const
		{ if (x < 0.5f) { const float x2 = x*x; return x2*(x2*0.25f-0.625f)+115.0f/192.0f; };
		  if (x < 1.5f) return x*(x*(x*(5.0f/6.0f - x*(1.0f/6.0f))-1.25f)+5.0f/24.0f)+55.0f/96.0f;
		  if (x < 2.5f) { float x2 = x-2.5f; x2 *= x2; return x2*x2*(1.0f/24.0f); };
		  return 0.0f; }
	float	weight3(const float x) const
		{ if (x < 1.0f) return (x*x*(x-2.0f)*3.0f+4.0f)/6.0f;
		  if (x < 2.0f) { const float x2 = 2.0f-x; return x2*x2*x2/6.0f; };
		  return 0.0f; }
	float	lastCoeff(const fvecD& c, const float z) const 
		{ return (z/(z*z-1))*(z*c(c.N-2)+c(c.N-1)); }
	void	setWeights3(float *wt, const float p) const				//! set weights for 3-point interpolation.
		{ for (unsigned int i = 0; i <= degree; i++)
			wt[i] = weight3(std::abs(p-i+degree/2)); }
	void	setWeights4(float *wt, const float p) const				//! set weights for 4-point interpolation.
		{ for (unsigned int i = 0; i <= degree; i++)
			wt[i] = weight4(std::abs(p-i+degree/2)); }
	float	firstCoeff(const fvecD& c, const float z) const;
	void	computeCoefficients(fvecD& c) const;
	void	setPositions(unsigned int *pi, const int p0, const unsigned int np) const;
public:
	bspline4(const fimage& src);
	float	interpolate(const fvec3 &p) const;
	fvec3	gradient(const fvec3 &p) const;
};

#endif

