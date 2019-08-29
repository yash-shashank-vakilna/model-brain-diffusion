#ifndef SPHARMDWI_H
#define SPHARMDWI_H

/*
 *
 * spharmDWI.h: represent DWI fibers using spherical harmonics
 * BRIAN Software Package Version 3.0
 *
 * $Id: spharmDWI.h 414 2016-10-08 00:02:53Z frithjof $
 *
 * 0.10 (03/08/14): first implementation
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Provides a light-weight class for a real spherical harmonics transformation.
*/

#define allDegs() (int l = 0; l <= int(deg); l += 2) for (int m = -l; m <= l; m++, i++)

//! Implements a light-weight class for a real spherical harmonics transformation.

class realSHT {	
	static const unsigned int nsp = 724;	//!< vertex# of standard sphere
	static const fvec3 sp724[];		//!< points of standard sphere
	static const unsigned int maxit = 50;	//!< maximum number of iterations in peak estimation
	const unsigned int deg;			//!< transformation degree
public:
	realSHT(const unsigned int _deg)
	//! allocates a context for a real spherical harmonics transformation of degree deg.
		: deg(_deg) { }
	realSHT(const fvecD& cf)
	//! allocates a context for a real spherical harmonics transformation for coefficient vector cf.
		: deg(degree(cf.N)) { }
	void	toPolar(const fvec3& v, float& theta, float& phi) const			//! converts v to polar coordinates.
		{ theta = float(acos(v.z)); phi = float(atan2(v.y, v.x)); }
	fvec3	toCart(const float theta, const float phi) const			//! converts (theta,phi) to cartesian coordinates.
		{ return fvec3(float(cos(phi)*sin(theta)),float(sin(phi)*sin(theta)),float(cos(theta))); }
	unsigned int degree(const unsigned int nc) const				//! returns order of SH decomposition with nc coefficients.
		{ return FTOU(std::sqrt(8.0f*nc+1.0f)-3)/2; }
	unsigned int nc() const								//! returns number of coefficients.
		{ return (deg+1)*(deg+2)/2; }
	float	radiusAt(const fvecD& cf, const float t, const float p)	const		//! computes SH radius at polar coordinates (t,p).
		{ unsigned int i = 0; float r = 0; for allDegs() r += cf(i)*YlmR(l,m,t,p);
		  return std::abs(r); }
	float	radiusAt(const fvecD& cf, const fvec3& q) const				//! computes SH radius at cartesian coordinates q.
		{ float t, p; toPolar(q,t,p); return radiusAt(cf,t,p); }
	fmatD	toSpace(const vfvec3& pt) const;
	fmatD	toCoeffs(const vfvec3& pt, const float lambda = 0) const;
	fmatD	sphereTransform() const;
	fmatD	sharpeningTransform(const float ratio) const;
	fvecD	sharpen(const fvecD& c, const fmatD& R, const fmatD& B, const float lb, const float tau) const;
	fvec3	maximizeRadius(const fvecD& cf, const fvec3& q) const;			
	vfvec3	searchPeaks(const fvecD& cf, const fmatD& B, const float alim = 0.5, const float thr = 0.25) const;
	fmatD	frt() const;
	float	generalizedFA(const fvecD& cf) const;
	float	varianceGFA(const fvecD& cf, const fmatD& B) const;
	fiberMixtureDesc* getFMD(const fvecD& cf, const fmatD& B) const;
};

// see also:
// Ivanic J., Ruedenberg K. (1996)
// Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion.
// J. Phys. Chem. 100, 6342-5347.
// Additions and Corrections (1998)
// J. Phys. Chem. 102, 9099.
// based on an implementation by D. Williamson.

//! Implements a context for rotating SHT coefficients in 3D.

class shtRotation {
	std::vector<fmatD> R;			//!< vector of rotation matrices

	float&	mat(const int i, const int j, const int k)				//! returns modifiable matrix coefficent.
		{ assert(i+j >= 0); assert(i+k >= 0);
		  return R[ITOU(i)](ITOU(j+i),ITOU(k+i)); }
	float	mat(const int i, const int j, const int k) const			//! returns constant matrix coefficent.
		{ assert(i+j >= 0); assert(i+k >= 0);
		  return R[ITOU(i)](ITOU(j+i),ITOU(k+i)); }
	float	delta(const int m, const int n) const					//! returns m == n.
		{ return m == n? 1: 0;	}
	float	P(const int i, const int l, const int a, const int b) const		//! returns P coefficient for matrix.
		{ float ri1 = mat(1,i,1), rim1 = mat(1,i,-1), ri0 = mat(1,i,0);
		  if (b == -l) return ri1*mat(l-1,a,-l+1)+rim1*mat(l-1,a,l-1);
		  else if (b == l) return ri1*mat(l-1,a,l-1)-rim1*mat(l-1,a,-l+1);
		  else return ri0*mat(l-1,a,b); }					// |b|<l
	float	U(const int l, const int m, const int n) const				//! returns U coefficient for matrix.
		{ return P(0,l,m,n); }
	float	V(const int l, const int m, const int n) const				//! returns V coefficient for matrix.
		{ if (m == 0) return P(1,l,1,n)+P(-1,l,-1,n);
		  else if (m > 0) { float p0 = P(1,l,m-1,n), p1 = P(-1,l,-m+1,n);
			float d = delta(m,1); return p0*std::sqrt(1+d)-p1*(1-d); }
		  else { float p0 = P(1,l,m+1,n), p1 = P(-1,l,-m-1,n);
			float d = delta(m,-1); return p0*(1-d)+p1*std::sqrt(1+d); } }
	float	W(const int l, const int m, const int n) const				//! returns W coefficient for matrix.
		{ if (m == 0) { assert(false); return 0; }
		  else if (m > 0) return P(1,l,m+1,n)+P(-1,l,-m-1,n);
		  else return P(1,l,m-1,n)-P(-1,l,-m+1,n); }
	float	M(const int l, const int m, const int n) const				//! returns M coefficient for matrix.
		{ float d = delta(m, 0); int ma = std::abs(m);
		  float den = std::abs(n) == l? (2*l)*(2*l-1): (l+n)*(l-n);
		  float u = std::sqrt((l+m)*(l-m)/den);
		  float v = 0.5f*std::sqrt((1+d)*(l+ma-1)*(l+ma)/den)*(1.0f-2.0f*d);
		  float w = -0.5f*std::sqrt((l-ma-1)*(l-ma)/den)*(1-d);
		  if (isSmall<float>(u) == false) u *= U(l,m,n);
		  if (isSmall<float>(v) == false) v *= V(l,m,n);
		  if (isSmall<float>(w) == false) w *= W(l,m,n);
		  return u+v+w; }
public:
	shtRotation(const unsigned int deg)
		//! allocates a context for rotating spherical harmonics coefficients of degree deg.
		: R(deg+1)
		{ for (unsigned int i = 0; i <= deg; i++) R[i] = fmatD(2*i+1,2*i+1); }
	void	computeRotation(const fmat3& r);
	const fmatD& level(const unsigned int i) const					//! returns coefficient matrix at level i.
		{ return R[i]; }
	fvecD	realRotate(const fvecD& cf) const;					//! rotates SPHARM coefficients cf and returns rotated coefficients.
};
#endif

