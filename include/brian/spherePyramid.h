#ifndef SPHEREPYRAMID_H
#define SPHEREPYRAMID_H

/*
 *
 * spherePyramid.h: linear interpolator for spherical triangle meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: spherePyramid.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (11/01/10): first implementation
 * 0.20 (21/01/10): spherePyramid added
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (20/12/13): documented
 * v302 (19/03/16): rewritten for sphereMesh
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Linear interpolator for spherical triangle meshes.
*/

//! Implements a multi-resolution pyramid for spherical meshes

class spherePyramid {
	const unsigned int nl;			//!< number of levels
	std::vector<sphereInterpolator> level;	//!< contains the meshes
public:
	spherePyramid(const sphereMesh& src, const unsigned int _nl)
	//! allocates a multi-resolution pyramid for mesh src with nl levels.
		 : nl(_nl), level()
		{ level.push_back({src});
		  for (unsigned int i = 1; i < nl; i++)
			level.push_back({level[i-1].downsample(0.5)}); }
	sphereInterpolator& operator[](const unsigned int i)				//! accesses level i.
		{ return level[i]; }
	const sphereInterpolator& operator[](const unsigned int i) const		//! accesses level i.
		{ return level[i]; }
	vertex* operator()(const unsigned int l, const unsigned int i) const
		{ return level[l].vertices(i); }
};

#endif

