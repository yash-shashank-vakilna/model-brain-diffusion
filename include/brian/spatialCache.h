#ifndef SPATIALCACHE_H
#define SPATIALCACHE_H

/*
 *
 * spatialCache.h: provide spatial queries on meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: spatialCache.h 422 2016-10-19 00:27:19Z frithjof $
 *
 * 0.10 (19/11/09): initial version
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (20/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for spatial queries on meshes.
*/

//! Implements a spatial cache for triangular meshes

template <class T> class spatialCache {
protected:
	std::list<T>**	data;			//!< pointer to cell simplices
	float		cellLength;		//!< length of a cell
	fvec3		origin;			//!< origin of mesh
	fvec3		ex;			//!< extent of mesh
	unsigned int	nx;			//!< extent in x direction
	unsigned int	ny;			//!< extent in y direction
	unsigned int	nz;			//!< extent in z direction
	fvec3		bmin;			//!< search range minimum of spatial query
	fvec3		bmax;			//!< search range maximum of spatial query
	unsigned int	ixmin;			//!< bounding box minimum index x	
	unsigned int	iymin;			//!< minimum index y	
	unsigned int	izmin;			//!< minimum index z	
	unsigned int	ixmax;			//!< bounding box maximum index x
	unsigned int	iymax;			//!< maximum index y	
	unsigned int	izmax;			//!< maximum index z	

	unsigned int nel() const							//! returns number of cells in cache.
		{ return nx*ny*nz; }
	unsigned int index(const unsigned int x, const unsigned int y, const unsigned int z)	//! returns index of cell (x,y,z).
		{ return (z*ny+y)*nx+x; }
	unsigned int index(const uvec3& s)						//! returns index of cell at site s.
		{ return (s.z*ny+s.y)*nx+s.x; }
	void	boundTo(const fvec3& a)							//! sets bounding box from position a.
		{ if (a.x < bmin.x) bmin.x = a.x;
		  if (a.y < bmin.y) bmin.y = a.y;
		  if (a.z < bmin.z) bmin.z = a.z;
		  if (a.x > bmax.x) bmax.x = a.x;
		  if (a.y > bmax.y) bmax.y = a.y;
		  if (a.z > bmax.z) bmax.z = a.z; }
	void	boundTo2(const fvec3& a, const fvec3& b)				//! sets bounding box from positions (a,b).
		{ bmin = std::numeric_limits<float>::max();
		  bmax = bmin*(-1.0f); boundTo(a); boundTo(b); }
	void	boundTo3(const fvec3& a, const fvec3& b, const fvec3& c)		//! sets bounding box from positions (a,b,c).
		{ boundTo2(a, b); boundTo(c); }
	void	boundTo4(const fvec3& a, const fvec3& b, const fvec3& c, const fvec3& d) //! sets bounding box from positions (a,b,c,d).
		{ boundTo3(a, b, c); boundTo(d); }
	void	setBounds();
	void	addSimplexToCell(const T t, const unsigned int i);
	void	removeSimplexFromCell(const T t, const unsigned int i);
	virtual void triangleBounds(const T t) = 0;
public:
	spatialCache(const fvec3 &_min, const fvec3 &_max, const float len)		//! allocates an empty cache with cell length len.
		 : data(nullptr), cellLength(len), origin(_min), ex(_max-_min),
		nx(FTOU(ex.x/cellLength+1)), ny(FTOU(ex.y/cellLength+1)), nz(FTOU(ex.z/cellLength+1)),
		bmin(), bmax(), ixmin(0), iymin(0), izmin(0), ixmax(0), iymax(0), izmax(0)
		{ data = new std::list<T>* [nx*ny*nz];
		  for (unsigned int i = 0; i < nx*ny*nz; i++) data[i] = nullptr; }
	spatialCache(const spatialCache& s) = delete;
	spatialCache(spatialCache&& s) = delete;
	virtual ~spatialCache() { clear(); }
	spatialCache& operator=(const spatialCache& s) = delete;
	spatialCache& operator=(spatialCache&& s) = delete;
	void	clear()									//! empties the cache.
		{ if (data == nullptr) return;
		  for (unsigned int i = 0; i < nel(); i++) delete data[i];
		  delete [] data; data = nullptr; }
	void	addTriangle(const T t);
	void	removeTriangle(const T t);
	std::list<T> collectTriangles();
};

template<typename T>
void spatialCache<T>::setBounds()
//! sets indices imin & imax from bounds.
{	int t = int((bmin.x-origin.x)/cellLength); ixmin = t < 0? 0: ITOU(t);
	t = int((bmin.y-origin.y)/cellLength); iymin = t < 0? 0: ITOU(t);
	t = int((bmin.z-origin.z)/cellLength); izmin = t < 0? 0: ITOU(t);
	t = int((bmax.x-origin.x)/cellLength); ixmax = t >= int(nx)? nx-1: ITOU(t);
	t = int((bmax.y-origin.y)/cellLength); iymax = t >= int(ny)? ny-1: ITOU(t);
	t = int((bmax.z-origin.z)/cellLength); izmax = t >= int(nz)? nz-1: ITOU(t);
}

template<typename T>
void spatialCache<T>::addSimplexToCell(const T t, const unsigned int i)
//! adds triangle t to cell i.
{	if (data[i] == nullptr) data[i] = new std::list<T>;				// allocate new list at i
	const auto di = data[i]; const auto it = std::find(di->begin(), di->end(), t);	// find element t in list
	if (it == di->end()) di->push_front(t);						// not found, so add to list
}

template<typename T>
void spatialCache<T>::removeSimplexFromCell(const T t, const unsigned int i)
//! removes triangle t from cell i.
{	const auto di = data[i]; if (di == nullptr) return;				// get triangle list at i
	const auto it = std::find(di->begin(), di->end(), t);				// find element t in list
	if (it != di->end()) di->erase(it);						// erase element
	if (di->empty()) { delete data[i]; data[i] = nullptr; };			// if this list is empty, remove it
}

template<typename T>
void spatialCache<T>::addTriangle(const T t)
//! adds triangle t to cache.
{	triangleBounds(t);
	for (unsigned int iz = izmin; iz <= izmax; iz++)				// for all lists in the bounding box
		for (unsigned int iy = iymin; iy <= iymax; iy++)
			for (unsigned int ix = ixmin; ix <= ixmax; ix++)
				addSimplexToCell(t, index(ix, iy, iz));			// add triangle t to list
}

template<typename T>
void spatialCache<T>::removeTriangle(const T t)
//! removes triangle t from cache.
{	triangleBounds(t);
	for (unsigned int iz = izmin; iz <= izmax; iz++)				// for all lists in the bounding box
		for (unsigned int iy = iymin; iy <= iymax; iy++)
			for (unsigned int ix = ixmin; ix <= ixmax; ix++)
				removeSimplexFromCell(t, index(ix, iy, iz));		// remove triangle t from list
}

template<typename T>
std::list<T> spatialCache<T>::collectTriangles()
//! returns a list of all triangles in the bounding box
{	setBounds(); std::list<T> tl;							// allocate empty list
	for (unsigned int iz = izmin; iz <= izmax; iz++)				// for all lists in the bounding box
		for (unsigned int iy = iymin; iy <= iymax; iy++)
			for (unsigned int ix = ixmin; ix <= ixmax; ix++) {
				const auto di = data[index(ix, iy, iz)]; if (di == nullptr) continue;
				tl.insert(tl.end(), di->begin(), di->end()); };		// join this list with tl
	tl.sort(); tl.unique();	return tl;						// return list of unique triangles
}

#endif


