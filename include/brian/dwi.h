#ifndef DWI_H
#define DWI_H

/*
 *
 * dwi.h: common definitions for computations on diffusion-weighted images
 * BRIAN Software Package Version 3.0
 *
 * $Id: dwi.h 466 2017-01-19 23:50:04Z kruggel $
 *
 * 0.10 (21/04/14): released version 2.6
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Common definitions for diffusion-weighted images.
*/

#define allDir(i)	(unsigned int i = 0; i < nd; i++)

//! Collects information about a diffusion-weighted image.

struct dwImage {
	vfimage& src;				//!< vector of (nd+1) source images
	const unsigned int nd;			//!< number of gradient directions
	const unsigned int nv;			//!< number of gradient directions
	const uvec3 ex;				//!< image extent
	const unsigned int nc;			//!< number of gradient coils
	const vfvec3 g;				//!< vector of gradient directions
	const fvecD b;				//!< vector of gradient strengths

	fvec3	gradientDirection(const unsigned int d) const				//! returns gradient direction from attribute.
		{ fvec3 v; const char *s = src[d].at.lookup("direction");
		  if (s && sscanf(s, "%f%f%f", &v.x, &v.y, &v.z) == 3) return v;
		  throw optException("Invalid/missing direction attribute"); }
	float	gradientStrength(const unsigned int d) const				//! returns gradient strength from attribute.
		{ float v = 0; const char *s = src[d].at.lookup("gradient_strength");
		  if (s && sscanf(s, "%f", &v) == 1) return v;
		  throw optException("Invalid/missing gradient_strength attribute"); }
	float	getDiffusionTime(const unsigned int d) const				//! returns diffusion time from attribute.
		{ float v = 0; const char *s = src[d].at.lookup("diffusion_time");
		  if (s && sscanf(s, "%f", &v) == 1) return v;
		  throw optException("Invalid/missing diffusion_time attribute"); }
	float	getPulseWidth(const unsigned int d) const				//! returns gradient pulse width from attribute.
		{ float v = 0; const char *s = src[d].at.lookup("pulse_width");
		  if (s && sscanf(s, "%f", &v) == 1) return v;
		  throw optException("Invalid/missing pulse_width attribute"); }
	vfvec3	initGradients(const unsigned int flipCode) const			//! returns a vector of all gradient vectors and optionally flips them.
		{ vfvec3 v(nd);
		  for allDir(d) { fvec3 gi = gradientDirection(d);
		  	if (flipCode == 1) gi.x = -gi.x;
		  	else if (flipCode == 2) gi.y = -gi.y;
		  	else if (flipCode == 3) gi.z = -gi.z;
		  	else if (flipCode == 4) { gi.x = -gi.x; gi.y = -gi.y; }
		  	else if (flipCode == 5) { gi.x = -gi.x; gi.z = -gi.z; }
		  	else if (flipCode == 6) { gi.y = -gi.y; gi.z = -gi.z; };
			v[d] = gi; };
		  return v; }
	fvecD	initStrengths() const							//! returns a vector of all gradient strengths.
		{ fvecD v(nd);
		  for allDir(d) v[d] = gradientStrength(d)/1000.0f;
		  return v; }
	fmatD	initBm() const								//! initializes the linear tensor estimation matrix.
		{ fmatD B(nd, 6);
		  for allDir(d) { 
			B(d,0) = b(d)*g[d].x*g[d].x; B(d,1) = 2*b(d)*g[d].x*g[d].y;
			B(d,2) = 2*b(d)*g[d].x*g[d].z; B(d,3) = b(d)*g[d].y*g[d].y;
			B(d,4) = 2*b(d)*g[d].y*g[d].z; B(d,5) = b(d)*g[d].z*g[d].z; }
		  return B; }
	fvecD	logData(const fvec3& s) const						//! returns negative normalized std::log data at site s.
		{ fvecD sm(nd); sm = 0.0f;
		  if (src[nd](s) <= 0.0f) return sm; 
		  const float ls0 = std::log(src[nd](s));
		  for allDir(d) sm[d] = -std::log(src[d](s))+ls0;
		  return sm; }
	fvecD	data(const fvec3& s) const						//! returns normalized data at site s.
		{ fvecD sm(nd); float sc = src[nd](s); sc = sc > 0.0f? 1.0f/sc: 0.0f;
		  for allDir(d) sm[d] = src[d](s)*sc;
		  return sm; }
	float	estimateNoise(const float alpha = 0.1f, const bool verbose = false) const;
	bimage	getMask() const;
	dwImage(vfimage& _src, const unsigned int _nc, const unsigned int flipCode)
		: src(_src), nd(ITOU(src.size()-1)), nv(src[nd].nel()), ex(src[nd].extent()),
		nc(_nc), g(initGradients(flipCode)), b(initStrengths())
		{  for allDir(d) { for allVoxels(i) { 
			if (src[d](i) <= 0) src[d](i) = 0.01f; } } }			// fix negative values inside mask
};

#endif

