#ifndef IMAGEPYRAMID_H
#define IMAGEPYRAMID_H

/*
 *
 * imagePyramid.h: build a multi-resolution image pyramid
 * BRIAN Software Package Version 3.0
 *
 * $Id: imagePyramid.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (14/03/09): initial version
 * 0.11 (22/10/09): filter sigma changed to 1.0
 * 0.20 (28/10/09): revised for changes in image format
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements a multi-resolution pyramid for images.
*/

//! Implements a multi-resolution pyramid for images

class imagePyramid {
	const float sigma;				//!< smoothing factor
	const unsigned int nl;				//!< number of levels
	vfimage lv;					//!< vector of multiresolution images

	fimage	downsample(fimage& src)							//! computes an image downsampled by factor 2.
		{ return src.filterGaussian(sigma).scale(src.extent()/2); }
public:
	imagePyramid(const fimage& src, const unsigned int _nl, const float _sigma)	//! constructs a multi-resolution pyramid with nl levels.
		 : sigma(_sigma), nl(_nl), lv(_nl)
		{ for (unsigned int i = 0; i < nl; i++)
			lv[i] = i? downsample(lv[i-1]): src; }
	const fimage& operator[] (const unsigned int i) const				//! access image at level i in pyramid.
		{ return lv[i]; }
};

#endif
