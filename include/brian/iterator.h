#ifndef ITERATOR_H
#define ITERATOR_H

/*
 *
 * iterator.h: neighborhood iterators for images
 * BRIAN Software Package Version 3.0
 *
 * $Id: iterator.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (26/10/09): initial version
 * 0.11 (13/11/09): corrections, implementation of local operators changed
 * 0.20 (10/04/10): box iterator added
 * 0.21 (27/04/10): c6p iterator added
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v319 (28/04/16): permIterator added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements neighborhood iterators for images.
*/

static const int dx4[]  = { -1,  0,  1,  0 };
static const int dy4[]  = {  0, -1,  0,  1 };
static const int dz4[]  = {  0,  0,  0,  0 };
static const int dx8[]  = { -1,  0,  1, -1,  1, -1,  0,  1 };
static const int dy8[]  = { -1, -1, -1,  0,  0,  1,  1,  1 };
static const int dz8[]  = {  0,  0,  0,  0,  0,  0,  0,  0 };
static const int dx6[]  = {  0,  0, -1,  1,  0,  0 };
static const int dy6[]  = {  0, -1,  0,  0,  1,  0 };
static const int dz6[]  = { -1,  0,  0,  0,  0,  1 };
static const int dx6p[] = {  0,  0, -1,  1,  0,  0, -1,  1, -1,  1, -1,  1, -1,  1 };
static const int dy6p[] = {  0, -1,  0,  0,  1,  0, -1, -1,  1,  1, -1, -1,  1,  1 };
static const int dz6p[] = { -1,  0,  0,  0,  0,  1, -1, -1, -1, -1,  1,  1,  1,  1 };
static const int dx18[] = {  0, -1,  0,  1,  0, -1,  0,  1, -1,  1, -1,  0,  1,  0, -1,  0,  1,  0 };
static const int dy18[] = { -1,  0,  0,  0,  1, -1, -1, -1,  0,  0,  1,  1,  1, -1,  0,  0,  0,  1 };
static const int dz18[] = { -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1 };
static const int dx26[] = { -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1 };
static const int dy26[] = { -1, -1, -1,  0,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  1,  1,  1, -1, -1, -1,  0,  0,  0,  1,  1,  1 };
static const int dz26[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1 };
static const int cex[] = { 0, 1, 0, 1, 0, 1, 0, 1 };
static const int cey[] = { 0, 0, 1, 1, 0, 0, 1, 1 };
static const int cez[] = { 0, 0, 0, 0, 1, 1, 1, 1 };

//! Helper class to keeps a voxel's site & depth

struct dvoxel {
	uvec3	s;				//!< voxel site
	float	d;				//!< voxel depth
	dvoxel()
		 : s(), d(0) { }
	dvoxel(const uvec3 _s, const float _d)						//! allocates a (site,depth) record.
		 : s(_s), d(_d) { }
	bool operator< (const dvoxel &b) const						//! comparison by increasing depth.
		{ return (this->d < b.d); }
};

//! Helper class to keeps a voxel properties

struct vxprop {
	unsigned int i;				//!< voxel index
	unsigned int c;				//!< voxel class
	float d;				//!< voxel depth or cost

	vxprop(const unsigned int _i = 0)						//! allocates an empty record for index i.
		 : i(_i), c(0), d(0) { }
	vxprop(const unsigned int _i, const float _d)					//! allocates a record for index i and depth d.
		 : i(_i), c(0), d(_d) { }
	bool operator< (const vxprop &b) const						//! comparison by increasing depth.
		{ return (this->d < b.d); }
};

//! Base class for image iterators

class imageIterator {
protected:
	uvec3	beg;				//!< first site
	uvec3	end;				//!< last site
	uvec3	loc;				//!< current site
	bool	forward(uvec3& s)							//! step forward by one site.
		{ if (loc.x <= end.x) { s = loc; loc.x++; return true; }
		  else { loc.x = beg.x; loc.y++; };
		  if (loc.y <= end.y) { s = loc; if (loc.x == beg.x) loc.x++;
			return true; }
		  else { loc.y = beg.y; loc.z++; };
		  if (loc.z <= end.z) { s = loc; if (loc.x == beg.x) loc.x++;
			return true; }
		  else return false; }
	bool	backward(uvec3& s)							//! step backward by one site.
		{ if (loc.x != beg.x-1) { s = loc; loc.x--; return true; }
		  else { loc.x = end.x; loc.y--; };
		  if (loc.y != beg.y-1) { s = loc; if (loc.x == end.x) loc.x--;
			return true; }
		else { loc.y = end.y; loc.z--; };
		  if (loc.z != beg.z-1) { s = loc; if (loc.x == end.x) loc.x--;
			return true; }
		  else return false; }
	imageIterator(const uvec3& ex, const unsigned int b)				//! allocates a image iterator for extent ex and boundary b.
		 : beg(b), end(ex-uvec3(b+1)), loc(0)
		{ if (ex.z == 1) { beg.z = 0; end.z = 0; }; }				// correct beg and end for 2D images
	virtual ~imageIterator() { }
};

//! Forward image iterator

struct fwIterator : public imageIterator { 
	fwIterator(const uvec3& ex, const unsigned int b)				//! allocates a forward iterator for extent ex and boundary b.
		 : imageIterator(ex, b) { loc = beg; }
	bool	operator()(uvec3& s)							//! returns current site in s and steps one site forward.
		{ return forward(s); }
};

//! Backward image iterator

struct bwIterator : public imageIterator {
	bwIterator(const uvec3& ex, const unsigned int b)				//! allocates a backward iterator for extent ex and boundary b.
		 : imageIterator(ex, b) { loc = end; }
	bool	operator()(uvec3& s)							//! returns current site in s and steps one site backward.
		{ return backward(s); }
};

//! Sub-sampling forward image iterator

class ssIterator : public imageIterator {	
	const unsigned int step;		//!< step size
public:
	ssIterator(const uvec3& ex, const unsigned int b, const unsigned int s)		//! allocates a subsampling iterator for extent ex, boundary b and step size s.
		 : imageIterator(ex, b), step(s) { loc = beg; }
	bool	operator()(uvec3& s)							//! returns current site in s and steps s sites forward.
		{ if (loc.x <= end.x) { s = loc; loc.x += step; return true; }
		  else { loc.x = beg.x; loc.y += step; };
		  if (loc.y <= end.y) { s = loc; if (loc.x == beg.x) loc.x += step;
			return true; }
		  else { loc.y = beg.y; loc.z += step; };
		  if (loc.z <= end.z) { s = loc; if (loc.x == beg.x) loc.x += step;
			return true; }
		  else return false; }
};

//! Iterates within a window

class boxIterator : public imageIterator {
public:
	boxIterator(const uvec3& ex, const unsigned int b, const uvec3& c, const unsigned int w) //! allocates a box iterator for extent ex, boundary b, around site c within window w.
		 : imageIterator(ex, b)
		{ unsigned int lx = c.x < w? 0: c.x-w;
		  beg.x = std::max(beg.x, lx); end.x = std::min(end.x, c.x+w);
		  unsigned int ly = c.y < w? 0: c.y-w;
		  beg.y = std::max(beg.y, ly); end.y = std::min(end.y, c.y+w);
		  unsigned int lz = c.z < w? 0: c.z-w;
		  beg.z = std::max(beg.z, lz); end.z = std::min(end.z, c.z+w); loc = beg; }
	bool	operator()(uvec3& s)							//! returns current site in s and steps one site forward.
		{ return forward(s); }
};

//! Base class for all neighborhood iterators

class nbIterator {
protected:
	const int *dx;				//!< neighborhood offset in x
	const int *dy;				//!< neighborhood offset in y
	const int *dz;				//!< neighborhood offset in z
	uvec3	ext;				//!< image extent
	uvec3	loc;				//!< current site
	unsigned int i;				//!< current index
	unsigned int n;				//!< neighborhood size

	void	setloc(const unsigned int t)						//! sets current site from index i.	
		{ const unsigned int r = t%(ext.x*ext.y);
		  loc.x = r%ext.x; loc.y = r/ext.x; loc.z = t/(ext.x*ext.y); }
	void	get(uvec3& s) const							//! returns neighbor site s.
		{ s.x = ITOU(int(loc.x)+dx[i]); s.y = ITOU(int(loc.y)+dy[i]); s.z = ITOU(int(loc.z)+dz[i]); }
	void	init(const connectivity cn)						//! initialize for connectivity cn.
		{ switch (cn) {
		  case connectivity::c4: n = 4; dx = dx4; dy = dy4; dz = dz4; break;
		  case connectivity::c8: n = 8; dx = dx8; dy = dy8; dz = dz8; break;
		  case connectivity::c6: n = 6; dx = dx6; dy = dy6; dz = dz6; break;
		  case connectivity::c6p: n = 14; dx = dx6p; dy = dy6p; dz = dz6p; break;
		  case connectivity::c18: n = 18; dx = dx18; dy = dy18; dz = dz18; break;
		  case connectivity::c26: n = 26; dx = dx26; dy = dy26; dz = dz26; break; 
          	  case connectivity::ce: n = 8; dx = cex; dy = cey; dz = cez; break; };
		  i = 0; }
	nbIterator(const uvec3& ex)							//! allocate a neighborhood iterator for extent ex.
		 : dx(0), dy(0), dz(0), ext(ex), loc(0), i(0), n(0) { }
	nbIterator(const uvec3& ex, const connectivity cn)				//! allocate a neighborhood iterator for extent ex and connectivity cn.
		 : dx(0), dy(0), dz(0), ext(ex), loc(0), i(0), n(0) { init(cn); }
public:
	nbIterator(const uvec3& ex, const uvec3& s, const connectivity cn)		//! allocate a neighborhood iterator for extent ex at site s for connectivity cn.
		 : dx(0), dy(0), dz(0), ext(ex), loc(s), i(0), n(0) { init(cn); }
	nbIterator(const nbIterator& b)							//! allocate a neighborhood iterator from model b.
		 : dx(b.dx), dy(b.dy), dz(b.dz), ext(b.ext), loc(b.loc), i(b.i), n(b.n) { }
	virtual ~nbIterator() { }
	bool	operator()(uvec3& s)							//! returns current neighbor in s and steps one neighbor forward.
		{ if (i >= n) return false; get(s); i++; return true; }
	bool	valid(const uvec3& s) const						//! checks if site s is within image extent.
		{ return s.x < ext.x && s.y < ext.y && s.z < ext.z; }
	nbIterator& operator=(const nbIterator& b)					//! assigns a neighborhood iterator from model b.
		{ if (this != &b) { i = b.i; n = b.n; dx = b.dx; dy = b.dy; dz = b.dz;
		  loc = b.loc; ext = b.ext; }; return *this; }
};

//! Special neighborhood iterator for a c6p configuration

class c6pIterator : public nbIterator {
	bool	valid(const unsigned int x, const unsigned int y, const unsigned int z) const	//! checks if site (x,y,z) is within image extent.
		{ return x < ext.x && y < ext.y && z < ext.z; }
	bool	c6pConfiguration(const image<bool>& f) const 				//! checks if current site has a c6p configuration in binary image f.
		{ if (i < 6) return true;						// normal c6
		  unsigned int x = ITOU(int(loc.x)+dx[i]),
			y = ITOU(int(loc.y)+dy[i]), z = ITOU(int(loc.z)+dz[i]);
		  if (valid(x,y,z) == false) return false;				// check if extreme position valid
		  if (f(x,loc.y,loc.z)^f(loc.x,y,z)) return false;			// if voxels on three diagonals
		  if (f(loc.x,y,loc.z)^f(x,loc.y,z)) return false;			// do not match, return
		  if (f(loc.x,loc.y,z)^f(x,y,loc.z)) return false;
		  return f(x,loc.y,loc.z)+f(loc.x,y,loc.z)+f(loc.x,loc.y,z) == 2; }	// one c6-connected voxel must be 0
public:
	c6pIterator(const uvec3& ex, const uvec3& s)					//! allocate a c6p neighborhood iterator for extent ex at site s.
		 : nbIterator(ex, s, connectivity::c6p) { }
	bool	operator()(uvec3& s, const image<bool>& f)				//! returns current neighbor in s and steps one neighbor forward.
		{ for ( ; i < n; i++) {
			if (c6pConfiguration(f)) { get(s); i++; return true; } }
		  return false; }
	bool	valid(const uvec3& s) const						//! checks if site s is within image extent.
		{ return s.x < ext.x && s.y < ext.y && s.z < ext.z; }
};

//! Iterator that finds foreground boundary points in binary images

class boundaryIterator : public nbIterator {
	bool	evaluate(const image<bool>& f)						//! checks if current site is on the boundary.
		{ if (f(loc) == false) return false; uvec3 s = loc;
		  for (i = 0; i < n; i++) {
			get(s); if (valid(s) && f(s) == false) return true; };
		  return false; }
public:
	boundaryIterator(const uvec3& ex, const connectivity cn)			//! allocates a boundary iterator for extent ex and connectivity cn.
		 : nbIterator(ex, cn) { }
	bool	operator()(const uvec3& s, const image<bool>& f)			//! checks if site s is on a boundary in binary image f.
		{ loc = s; return evaluate(f); }
	bool	operator()(const unsigned int i, const image<bool>& f)			//! checks if index i is on a boundary in binary image f.
		{ setloc(i); return evaluate(f); }
};

//! Implements an iterator for binary 8-cell

class c8Iterator : public nbIterator {
	bool	valid(const fvec3& p) const						//! checks if point p is within extent.	
		{ return p.x >= 0 && p.x < ext.x-1 && p.y >= 0 && p.y < ext.y-1
			&& p.z >= 0 && p.z < ext.z-1; }
	bool	valid(const uvec3& s) const						//! checks if site s is within extent.
		{ return s.x < ext.x-1 && s.y < ext.y-1 && s.z < ext.z-1; }
public:
	c8Iterator(const uvec3& ex)							//! allocates a cell iterator for extent ex.
		 : nbIterator(ex) { n = 8; dx = cex; dy = cey; dz = cez; }
	bool	interpolateAt(const fvec3& p, const image<bool>& f)			//! interpolates a binary image f at point p.
		{ if (valid(p) == false) return false; float wt[8], v = 0; uvec3 s;
		  loc = uvec3(FTOU(p.x),FTOU(p.y),FTOU(p.z));
		  float dx1 = p.x-float(loc.x), dx2 = 1.0f-dx1, dy1 = p.y-float(loc.y);
		  float dy2 = 1.0f-dy1, dz1 = p.z-float(loc.z), dz2 = 1.0f-dz1;
		  wt[0] = dx2*dy2*dz2; wt[1] = dx1*dy2*dz2; wt[2] = dx2*dy1*dz2;
		  wt[3] = dx1*dy1*dz2; wt[4] = dx2*dy2*dz1; wt[5] = dx1*dy2*dz1;
		  wt[6] = dx2*dy1*dz1; wt[7] = dx1*dy1*dz1;
		  for (i = 0; i < 8; i++) { get(s); v += wt[i]*f(s); };
		  return (v >= 4.0); }
	float	average(uvec3 s, const image<float>& f)					//! returns the average of binary image f around site s.
		{ if (valid(s) == false) return false; float v = 0; loc = s;
		  for (i = 0; i < 8; i++) { get(s); v += f(s); }; return v/8.0f; }
	int	euler(uvec3 s, const connectivity cn, const image<bool>& f)		//! returns the Euler number at site s in binary image f for connectivity cn.
		{ if (valid(s) == false) return false; int v[8], n[8], p = 0, q = 0;
		  loc = s; for (i = 0; i < 8; i++) { get(s); v[i] = f(s); n[i] = 1-v[i]; p += v[i]; };
		  if (p == 8 || p == 0) return 0;
		  switch (cn)  {
		  case connectivity::c6: p = (p == 8);
			q = v[0]*(1-v[1]*(1-v[4]*v[5])-v[2]*(1-v[1]*v[3])-v[4]*(1-v[2]*v[6]))-p;
			break;
		  case connectivity::c6p: p = (p == 8);
			q = v[0]*(1-v[1]*(1-v[4]*v[5])-v[2]*(1-v[1]*v[3])-v[4]*(1-v[2]*v[6]))-p;
			q += n[0]*v[1]*v[2]*v[3]*v[4]*v[5]*v[6]*n[7];
			q += v[0]*n[1]*v[2]*v[3]*v[4]*v[5]*n[6]*v[7];
			q += v[0]*v[1]*n[2]*v[3]*v[4]*n[5]*v[6]*v[7];
			q += v[0]*v[1]*v[2]*n[3]*n[4]*v[5]*v[6]*v[7]; break;
		  case connectivity::c18: p = (p == 0);
			q = n[0]*(1-n[1]*(1-n[4]*n[5])-n[2]*(1-n[1]*n[3])-n[4]*(1-n[2]*n[6]))-p;
			q += v[0]*n[1]*n[2]*n[3]*n[4]*n[5]*n[6]*v[7];
			q += n[0]*v[1]*n[2]*n[3]*n[4]*n[5]*v[6]*n[7];
			q += n[0]*n[1]*v[2]*n[3]*n[4]*v[5]*n[6]*n[7];
			q += n[0]*n[1]*n[2]*v[3]*v[4]*n[5]*n[6]*n[7]; break;
		  case connectivity::c26: p = (p == 0);
			q = n[0]*(1-n[1]*(1-n[4]*n[5])-n[2]*(1-n[1]*n[3])-n[4]*(1-n[2]*n[6]))-p;
			break;
		  case connectivity::c4:
		  case connectivity::c8:
		  case connectivity::ce:
			throw optException("Invalid connectivity %d", static_cast<int>(cn)); };
		  return q; }
};

//! Implements an iterator for generating permutations of unsigned integers

class permIterator {
	std::list<unsigned int> l;
	bool	first;
public:
	permIterator(const unsigned int n)						//! allocates a permutation iterator for integers { 0, n-1 }.
		: l(), first(true)
		{ for (unsigned int i = 0; i < n; i++) l.push_back(i); }
	bool	operator()(std::list<unsigned int>& p)					//! returns a permutated list in p & and true if permutation is valid.
		{ const bool ok = first? true: std::next_permutation(l.begin(),l.end());
		  first = false; p = l; return ok; }
};


#endif





