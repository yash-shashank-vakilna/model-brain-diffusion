#ifndef COOC_H
#define COOC_H

/*
 *
 * cooc.h: functions for cooccurrence matrices
 * BRIAN Software Package Version 3.0
 *
 * $Id: cooc.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (20/01/10): adapted to BRIAN2
 * 0.20 (10/12/12): const added
 * 0.40 (16/12/13): documented
 * 0.41 (13/10/15): minor bug fixes, asVector() added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Provides a class for a multi-sort co-occurrence matrix.
*/

//! Implements a multi-sort co-occurrence matrix.

class cooc {
	const unsigned int ni;				//!< intensity bins
	const unsigned int ng;				//!< gradient bins
	const unsigned int na;				//!< angle bins
	unsigned int ns;				//!< number of samples
	float 	li, hi;					//!< intensity limits
	float 	lg, hg;					//!< gradient limits
	fmatD	ii;					//!< intensity-intensity matrix
	fmatD	gg;					//!< gradient-gradient matrix
	fvecD	a;					//!< angle vector

	unsigned int d2bin(const float d) const						//! convert angle d to bin.
		{ const float d0 = CLAMP(d,-1.0f,1.0f);
		  unsigned int b = FTOU((acos(d0)*na)/M_PI);
		  if (b >= na) b = na-1; return b; };
	void 	init()									//! (re)initialize structures.
		{ ns = 0; ii = 0.0, gg = 0.0; a = 0; };
	bool	intbin(const float i, unsigned int& b) const				//! convert intensity i to bin b.
		{ if (i < li || i >= hi) return false;
		  b = FTOU((ni*(i-li))/(hi-li)); return true; };
	bool	intbin(const float i1, const float i2, unsigned int& b) const		//! converts intensity (i1,i2) to bin b.
		{ unsigned int b1, b2;
		  if (intbin(i1, b1) == false && intbin(i2, b2) == false) return false;
		  b = b1*ni+b2; return true; };
	bool	gradbin(const float g, unsigned int& b) const				//! convert gradient g to bin b.
		{ if (g < lg || g >= hg) return false;
		  b = FTOU((ng*(g-lg))/(hg-lg)); return true; };
	int	angbin(const fvec3& v1, const fvec3& v2) const				//! convert angle (v1,v2) to bin.
		{ fvec3 d1 = v1.normalize(), d2 = v2.normalize();
		  return d2bin(dot(d1,d2)); };
public:
	cooc(const unsigned int _ni, const unsigned int _ng, const unsigned int _na, 
		const float _li, const float _hi, const float _lg, const float _hg)
		//! init a co-occurrence matrix. Parameters
		//! ni (ng, na) correspond to the number of intensity (gradient, angle) bins
		//! and (li,hi) to the intensity and (lg,hg) to the gradient ranges.
		: ni(_ni), ng(_ng), na(_na), ns(0), li(_li), hi(_hi), lg(_lg), hg(_hg), 
		ii(ni, ni), gg(ng, ng), a(na) { init(); };
	unsigned int sampleAt(const fimage& src, const fimage& gmag, const fvec3image& grad,
			const uvec3& s, const unsigned int w, const bimage& mask);
	void	print() const;
	void	print(const unsigned int l) const					//! print co-occurrence matrix, preceeded by l.
		{ printf("%1d ", l); print(); printf("\n"); }
	float	getint(const unsigned int i, const unsigned int j) const		//! returns number of intensity co-occurrences at (i,j).
		{ return (i >= ni || j >= ni)? 0: ii(i, j); };
	float	getgrad(const unsigned int i, const unsigned int j) const		//! returns number of gradient co-occurrences at (i,j).
		{ return (i >= ng || j >= ng)? 0: gg(i, j); };
	float	getang(const unsigned int i) const					//! returns number of angle co-occurrences at (i,j).
		{ return (i >= na)? 0: a(i); };
	void	add(const cooc& c)							//! add a co-occurrence matrix to this one.
		{ ii += c.ii; gg += c.gg; a += c.a; };
	void	div(const float d)							//! normalize co-occurrence matrix by d.
		{ if (d == 0) return; ii /= d; gg /= d; a /= d; };
	float	dist(const cooc& c) const;
	void	entropy(float& ei, float& eg, float& ea) const;
	float	map(const fmatD& r, const fvecD& m, const fmatD& c) const;
	float	jointMap(const cooc& c, const fmatD& r, const fvecD& mn, const fmatD& cv) const;
	float	map(fmatD& r, fvecD& m, fmatD& c) const;
	fvecD	asVector() const;
};

static unsigned int delta[][3] = { { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, 
	{ 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 }};

unsigned int cooc::sampleAt(const fimage& src, const fimage& gmag, const fvec3image& grad,
			const uvec3& s, const unsigned int w, const bimage& mask)
//! samples a co-occurrence matrix from intensity image src, gradient magnitude image gmag,
//! and gradient image grad in a window of half-width w around site s, taking only samples
//! within a binary mask into account. All images are std::expected to have the same dimensions.
{	init(); if (mask(s) == 0) return 0;
	unsigned int ib1, ib2, gb1, gb2;
	for (unsigned int z = s.z-w; z <= s.z+w; z++)  {				// for a window w around s
		for (unsigned int y = s.y-w; y <= s.y+w; y++)  {
			for (unsigned int x = s.x-w; x <= s.x+w; x++)  {
				if (mask(x,y,z) == false) continue;			// check site against mask
				const float i1 = src(x,y,z);
				if (intbin(i1, ib1) == false) continue;			// check site's intensity range
				const float g1 = gmag(x,y,z);
				if (gradbin(g1, gb1) == false) continue;		// check site's gradient range
				for (unsigned int d = 0; d < 7; d++)  {			// for all right-sided neighbors
					const unsigned int dx = x+delta[d][0]; if (dx >= src.nx) continue;
					const unsigned int dy = y+delta[d][1]; if (dy >= src.ny) continue;
					const unsigned int dz = z+delta[d][2]; if (dz >= src.nz) continue;
					if (mask(dx,dy,dz) == false) continue;		// check neighbor against mask
					float i2 = src(dx,dy,dz);
					if (intbin(i2, ib2) == false) continue;		// check neighbor's intensity range
					float g2 = gmag(dx,dy,dz);
					if (gradbin(g2, gb2) == false) continue;	// check neighbor's gradient range
					ii(ib1, ib2) += 1.0; gg(gb1, gb2) += 1.0;	// sample this co-occurrence
					float d1 = dot(grad(x,y,z),grad(dx,dy,dz));
					a(d2bin(d1)) += 1.0; ns++; } } } };
	div(ns); return ns;								// normalize by sample count
}

float cooc::dist(const cooc& c) const
//! computes the distance between two co-occurrence matrices.
{	if (ns < 120) return 0; float si = 0, sg = 0, sa = 0; 
	for (unsigned int i = 0; i < ni; i++)
		for (unsigned int j = 0; j < ni; j++) si += std::abs(ii(i,j)-c.ii(i,j));
	for (unsigned int i = 0; i < ng; i++)
		for (unsigned int j = 0; j < ng; j++) sg += std::abs(gg(i,j)-c.gg(i,j));
	for (unsigned int i = 0; i < na; i++) sa += std::abs(a(i)-c.a(i));
	return si+sg+sa;
}

void cooc::print() const
//! prints a co-occurrence matrix.
{	for (unsigned int i = 0; i < ni; i++)
		for (unsigned int j = 0; j < ni; j++) printf("%6.4f ", ii(i,j));
	for (unsigned int i = 0; i < ng; i++)
		for (unsigned int j = 0; j < ng; j++) printf("%6.4f ", gg(i,j));
	for (unsigned int i = 0; i < na; i++) printf("%6.4f ", a(i));
}

void cooc::entropy(float& ei, float& eg, float& ea) const
//! returns the entropy of the intensity, gradient and angle observations.
{	ei = 0; eg = 0; ea = 0; 
	for (unsigned int i = 0; i < ni; i++)
		for (unsigned int j = 0; j < ni; j++) if (ii(i,j)) ei -= std::log(ii(i,j));
	for (unsigned int i = 0; i < ni; i++)
		for (unsigned int j = 0; j < ni; j++) if (gg(i,j)) eg -= std::log(gg(i,j));
	for (unsigned int i = 0; i < na; i++) if (a(i)) ea -= std::log(a(i));
}

float cooc::map(const fmatD& r, const fvecD& m, const fmatD& c) const
//! maps this co-occurrence matrix onto a vector and computes distance to (r,m,c).
{	fvecD v(ni*ni+ng*ng+na); unsigned int k = 0;
	for (unsigned int i = 0; i < ni; i++) 
		for (unsigned int j = 0; j < ni; j++) v(k++) = ii(i,j);
	for (unsigned int i = 0; i < ng; i++) 
		for (unsigned int j = 0; j < ng; j++) v(k++) = gg(i,j);
	for (unsigned int i = 0; i < na; i++) v(k++) = a(i);
	const fvecD d = c*(r*v-m); return norm2(d);
}

float cooc::jointMap(const cooc& c, const fmatD& r, const fvecD& mn, const fmatD& cv) const
//! computes similarity between this map and c with respect to (r,m,c).
{	fvecD v(2*(ni*ni+ng*ng+na)); unsigned int k = 0;
	for (unsigned int i = 0; i < ni; i++) 
		for (unsigned int j = 0; j < ni; j++) v(k++) = ii(i,j);
	for (unsigned int i = 0; i < ng; i++) 
		for (unsigned int j = 0; j < ng; j++) v(k++) = gg(i,j);
	for (unsigned int i = 0; i < na; i++) v(k++) = a(i);
	for (unsigned int i = 0; i < ni; i++) 
		for (unsigned int j = 0; j < ni; j++) v(k++) = c.ii(i,j);
	for (unsigned int i = 0; i < ng; i++) 
		for (unsigned int j = 0; j < ng; j++) v(k++) = c.gg(i,j);
	for (unsigned int i = 0; i < na; i++) v(k++) = c.a(i);
	const fvecD d = cv*(r*v-mn); return norm2(d);
}

fvecD cooc::asVector() const
//! returns co-occurence matrix as feature vector.
{	fvecD v(ni*ni+ng*ng+na); unsigned int t = 0;
	for (unsigned int i = 0; i < ni; i++)
		for (unsigned int j = 0; j < ni; j++) v[t++] = ii(i,j);
	for (unsigned int i = 0; i < ng; i++)
		for (unsigned int j = 0; j < ng; j++) v[t++] = gg(i,j);
	for (unsigned int i = 0; i < na; i++) v[t++] = a(i);
	return v;
}

#endif
