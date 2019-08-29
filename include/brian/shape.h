#ifndef SHAPE_H
#define SHAPE_H

/*
 *
 * shape.h: shape descriptors of binary objects
 * BRIAN Software Package Version 3.0
 *
 * $Id: shape.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (26/02/01): initial release
 * 0.20 (28/11/09): adapted for BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Shape descriptors of binary objects.
*/

//! Collects shape information about a binary object

class shape {
public:
	float	nv;				//!< number of pixels in the entire image
	float	nobj;				//!< number of pixels in the object
	float	pc_vol;				//!< ratio nobj / npix
	fvec3	cent;				//!< center of object
	fvec3	sigma;				//!< variance
	fvec3	voxel;				//!< voxel dimensions
	fvec3	sigmaM;				//!< variance on principle axis
	fvec3	theta;				//!< angle between axes
	float	m000,m001,m010,m100,m011;	//!< un-normalized moments
	float	m101,m110,m002,m020,m200;	//!< un-normalized moments
	float	mu000,mu001,mu010,mu100,mu011;	//!< normalized moments
	float	mu101,mu110,mu002,mu020,mu200;	//!< normalized moments

	void	print(FILE *outf, const unsigned int l) const				//! prints shape information to stdout.
		{ if (nobj > 0) fprintf(outf, "%6.0u %8.0f %8.1lf %8.1lf %8.1lf %8.2lf %8.2lf %8.2lf\n",
			l, nobj, cent.x, cent.y, cent.z, sigmaM.x, sigmaM.y, sigmaM.z); }
	void	add(const uvec3& s)							//! adds voxel s to shape.
		{ m000 += 1.0f; m001 += s.x; m010 += s.y; m100 += s.z;
		  m011 += s.y * s.x; m101 += s.z * s.x;
		  m110 += s.z * s.y; m002 += s.x * s.x;
		  m020 += s.y * s.y; m200 += s.z * s.z; }
	void	eval(const float np)							//! evaluates a shape record.
		{ // whole volume, object volume and percentage of object
		  nv = np; nobj = m000; pc_vol = nobj*100/nv; if (m000 == 0) return;
		  cent = fvec3(m001, m010, m100)/m000;					// object center and second-order moments
		  mu002 = m002-m000*cent.x*cent.x;
		  mu020 = m020-m000*cent.y*cent.y;
		  mu200 = m200-m000*cent.z*cent.z;
		  mu011 = m011-m000*cent.y*cent.x;
		  mu101 = m101-m000*cent.z*cent.x;
		  mu110 = m110-m000*cent.z*cent.y;
		  sigma.x = std::sqrt(mu002/m000); sigma.y = std::sqrt(mu020/m000); sigma.z = std::sqrt(mu200/m000);
		  theta.x = almostEqual(mu020,mu002)? float(0.25*M_PI): 0.5f*std::atan(2*mu011/(mu002-mu020));
		  theta.y = almostEqual(mu200,mu020)? float(0.25*M_PI): 0.5f*std::atan(2*mu110/(mu020-mu200));
		  theta.z = almostEqual(mu200,mu002)? float(0.25*M_PI): 0.5f*std::atan(2*mu101/(mu002-mu200));
		  float cos2 = std::cos(theta.x)*std::cos(theta.x);
		  float sin2 = std::sin(theta.x)*std::sin(theta.x);
		  float sin_2 = std::sin(2*theta.x);
		  sigmaM.x = std::sqrt((cos2*mu002+sin2*mu020+sin_2*mu011)/m000);	// computes the moments variance on both principal axes
		  sigmaM.y = std::sqrt((sin2*mu002+cos2*mu020-sin_2*mu011)/m000);
		  cos2 = std::cos(theta.z)*std::cos(theta.z);
		  sin2 = std::sin(theta.z)*std::sin(theta.z);
		  sin_2 = std::sin(2*theta.z);
		  sigmaM.z = std::sqrt((sin2*mu002+cos2*mu200-sin_2*mu101)/m000); }
	float	size() const								//! returns the size of a shape.
		{ return nobj; }
	shape()										//! allocates an empty shape record.
		 : nv(0),nobj(0),pc_vol(0),cent(0),sigma(0),voxel(1),sigmaM(0),theta(0),
		  m000(0),m001(0),m010(0),m100(0),m011(0), m101(0),m110(0),m002(0),m020(0),
		  m200(0),mu000(0),mu001(0),mu010(0),mu100(0),mu011(0),mu101(0),mu110(0),
		  mu002(0),mu020(0),mu200(0) { }
};

#endif
