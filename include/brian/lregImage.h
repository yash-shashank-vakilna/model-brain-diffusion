#ifndef __LREGIMAGE_H
#define __LREGIMAGE_H

/*
 *
 * lregImage.h: basic definitions & functions for linear image registration
 * BRIAN Software Package Version 3.0
 *
 * $Id: lregImage.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (28/06/00): first release
 * 0.20 (09/08/02): improved final optimization - FK
 * 0.21 (31/10/07): different interface to genetic optimizer & cross correlation function corrected - FK
 * 0.22 (23/04/09): rewritten for BRIAN 1.20 - FK
 * 0.23 (16/10/09): adapted to revised vector format
 * 0.24 (28/10/09): adapted to revised image format
 * 0.25 (01/12/09): tr.z capturing range changed
 * 0.26 (22/06/10): rt.x & rt.y capturing range reduced
 * 0.27 (31/08/10): optimization removed from manual plan
 * 0.28 (06/05/11): bugfix for scaling & payload implemented
 * 0.29 (21/02/12): combined with linregDWI()
 * 0.30 (06/01/13): for BRIAN 2.4
 * 0.40 (18/06/13): linregDWI revised for better handling of iso images FK
 * 0.50 (19/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions & functions for linear image registration.
*/

class lregImage;
using costfn = float (lregImage::*)(const fvecD&) const;

//! Implements a context for linear image registration in 3D.

class lregImage  {
	const unsigned int nlevel = 4;		//!< multi-resolution levels
	const imagePyramid sp;			//!< source multi-resolution pyramid
	const imagePyramid rp;			//!< reference multi-resolution pyramid
	const bool verbose;			//!< print life information
	const bool scf;				//!< scaling wanted
	fvecD	pos;				//!< transformation parameter vector
	fvecD	range;				//!< transformation parameter range
	unsigned int level;			//!< current resolution level
	unsigned int np;			//!< number of parameters
	fmat3	ms;				//!< source matrix
	fmat3	mr;				//!< reference matrix
	fmat3	mri;				//!< inverse reference matrix
	costfn	cost;				//!< similarity function

	void	param2mat(const fvecD& p, fmat3& m) const				//! converts parameter vector pos to matrix m.
		{ m.id(); fvec3 t(p(0), p(1), p(2)); m.translate(t);
		  fvec3 r(p(3), p(4), p(5)); m.rotate(r);
		  if (scf) { fvec3 s(p(6), p(7), p(8)); m.scale(s); };	  
		  m = (ms*m)*mri; float f = std::pow(2.0f,-float(level)); m.scaleTranslation(f); }
	void	print() const;
	void	geneticOptimize(const unsigned int s, const unsigned int n);
	void	simplexOptimize();
	void	optimize();
public:
	lregImage(const fimage& _src, const fimage& _ref, const bool _scf = false,
			const bool _verbose = false);
	float	similarityCC(const fvecD& p) const					//! cross-correlation similarity function.
		{ fmat3 m; param2mat(p,m); return sp[level].cc(rp[level],m); }
	float	similarityNMI(const fvecD& p) const					//! normalized mutual information similarity function.
		{ fmat3 m; param2mat(p,m); return sp[level].nmi(rp[level],m); }
	void	work(const unsigned int plan, const fvec3& tr, const fvec3& rt, const fvec3& sc);
	fvecD	map(const fvecD& val)
		const { return val*range+pos; }
	fimage	transform(const fimage& src) const;
	limage	transformLabel(const limage& pay) const;
	fmat3	getTransformation() const
		{ fmat3 m; param2mat(pos, m); return m; }
};

#endif
