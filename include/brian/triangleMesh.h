#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

/*
 *
 * triangleMesh.h: Basic definitions and classes for handling triangulated meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: triangleMesh.h 461 2017-01-15 19:43:59Z frithjof $
 *
 * v400 (20/09/16): extracted for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions and classes for handling triangulated meshes.
*/

/*! \enum smootherType
    \brief Symbolic constants for smoothing operators for triangular meshes.
*/
enum class smootherType {
	msExpLaplacian,				///< iterative smoothing using Laplacian weights
	msExpTaubin,				///< iterative smoothing using Taubins method
	msExpFujiwara,				///< iterative smoothing using edge weights
	msImpFujiwara,				///< direct smoothing using edge weights
	msExpDesbrun,				///< iterative smoothing using curvature weights
	msImpDesbrun,				///< direct smoothing using curvature weights
	msInflate,				///< inflate mesh a la "freesurfer"
};

enum class enType { Friedel, Degener, Clarenz };
enum class cvType { Gaussian, Mean, ShapeIndex };


//! Helper structure for a mesh simplification command.

struct cmd {
	char	type;				//!< command type
	unsigned int id0, id1;			//!< parameters

	cmd(const char t, const unsigned int i0, const unsigned int i1 = 0)		//! allocates a command record.
		 : type(t), id0(i0), id1(i1) { }
};

using cmdlist = std::deque<cmd>;

//! Represents triangular surface meshes in 3D.

class sphereMesh;

struct triangleMesh : public mesh {
// --- const functions assessing connectivity
	bool	cornerPoints(const unsigned int a, const unsigned int b, 
			unsigned int &p1, unsigned int &p2) const; 
	unsigned int oppositeVertex(const unsigned int i, const unsigned int j,
			const unsigned int a) const;
	std::vector<unsigned int> cornerPoints(const unsigned int a,
			const unsigned int b) const;
	bool	trianglesAtEdge(const unsigned int a, const unsigned int b,
			unsigned int& tid1, unsigned int& tid2) const;
	unsigned int countBoundaryTrianglesAt(const unsigned int a) const;
	bool	vertexHasTriangle(const unsigned int i, const unsigned int a,
			const unsigned int b) const;
 	unsigned int triangleAtEdge(const unsigned int a, const unsigned int b) const;
// --- const functions assessing vertex positions
	double	cotanAlpha(const unsigned int i, const unsigned int j,
			const unsigned int a) const;
	double	conformalWeight(const unsigned int i, const unsigned int j) const;	// returns weight for conformal mapping
	double	authalicWeight(const unsigned int i, const unsigned int j) const;	// returns weight for authalic mapping
	fvec3	center(const vertex* v) const;
	float	LaplacianAt(const vertex* v) const;
	fvec3	LaplacianFujiwara(const unsigned int i) const;
	void	curvatureTensor(const unsigned int id, const fvec3& ni,
			float& k1, float& k2, fvec3& v1, fvec3& v2) const;
	void	getCurvature(fvecD& K, fvecD& H, const vfvec3& nm,
			const float pc = 0.002f, const float d = 4.0f) const;
	float	ringArea(const unsigned int i) const;					// computes area of triangles in one-ring
	float	vertexArea(const unsigned int id) const;
	bool	checkRing(const unsigned int i) const;
	void	checkTriangleQuality(const float area_lim, const float ang_lim) const;
	fvecD	computeDepth(const fimage& dtf) const;
// --- functions changing mesh connectivity
	void	linkTriangleEdges(const bool verbose);
	bool	testEdgeCount(const bool verbose);
	bool	splitLongestEdge(const unsigned int id0);
	bool	flipEdgeAt(const unsigned int a, const unsigned int b);
	bool	weedTriangle(const unsigned int i);
public:
   	triangleMesh()									//! allocates an empty triangular mesh.
		: mesh() { }
   	triangleMesh(const vgraph& v, const pgraph& t)					//! allocates a triangular mesh from a vertex and triangle graph.
		: mesh(v, t) { }
   	triangleMesh(const triangleMesh& m)						//! copies from triangle mesh m.
		: mesh(m) { }
   	triangleMesh(triangleMesh&& m)							//! moves a triangular mesh from m.
		: mesh(std::move(m)) { }
   	triangleMesh(FILE* fp)								//! allocates an triangular mesh from FILE* fp.
		: mesh() { read(fp); }
 	triangleMesh& operator=(const triangleMesh& b)					//! assigns triangle mesh from m.
		{ if (this != &b) mesh::operator=(b); return *this; }
 	triangleMesh& operator=(triangleMesh&& b)					//! move assigns triangle mesh from m.
		{ assert(this != &b); mesh::operator=(std::move(b)); return *this; }
	triangle* tr(const unsigned int i) const					//! returns triangle i.
		{ return reinterpret_cast<triangle*>(primitives(i)); }
	void	read(FILE* fp);
	void	save(FILE* fp);
	void	read(const char* fname)							//! reads mesh from file named fname.
		{ mesh::read(fname); }
	void	save(const char* fname)							//! reads mesh from file named fname.
		{ mesh::save(fname); }
// --- const functions operating on connectivity
	int	genus() const;
// --- functions assessing vertex positions
	fvec3	center() const								//! returns the center of the mesh.
		{ fvec3 c = 0.0f; unsigned int n = 0; 
		  for allVertices(i) { if (vertices(i)) { c += vertices(i)->pos; n++; } };
		  return n != 0? c/n: c; }
	float	totalArea() const							//! returns the surface area.
		{ float a = 0.0f; 
		  for allPrimitives(i) { if (tr(i)) a += tr(i)->area(vertices); };
		  return a; }
	float	angleError(const triangleMesh& ref) const;
	float	lengthError(const triangleMesh& ref) const;
	float	areaError(const triangleMesh& ref) const;
	unsigned int checkProperties(const bool verbose = false) const;
	void	averageMeasure(fvecD& C, const unsigned int nit) const;
	bool	checkTriangleProperties(const unsigned int i, const bool verbose = false) const;
	fvecD 	shapeIndex();
	fvecD 	curvatureK();
	fvecD 	curvatureH();
	fvecD	curvatureR();
	vfvec3	computeNormals() const;
	void	computeBasis(vfvec3& b0, vfvec3& b1) const;
	bimage	voxelize(const connectivity cn, const uvec3 ex) const;
	fvecD	vertexArea();
	fvecD	lumpedArea() const;
	float	surfaceDist(const unsigned int s, const unsigned int e) const;
	fvecD	smoothedCurvature(const cvType& cv, const unsigned int nit = 5);
	triangleMesh regAffine(const triangleMesh& ref, const unsigned nit = 100) const;
// --- functions changing vertex positions
	void	unitScale(const float s = 4.0f*M_PI)					//! maps to origin and scales to unit sphere.
		{ const fvec3 c = center(); const float f = std::sqrt(s/totalArea());
		  for allVertices(i) { vertex* v = vertices(i); if (v == nullptr) continue;
		  	v->pos = (v->pos-c)*f; } };
	void	move(const fvec3& c)							//! moves mesh to center c.
		{ for allVertices(i) { vertex* v = vertices(i); if (v) v->pos -= c; } }
	void	addVertex(vertex* v, const unsigned int i)				//! adds vertex* v at i.
		{ vertices.addNodeAt(v,i); }
	void	addTriangle(triangle* t, const unsigned int i)				//! adds triangle* t at i.
		{ primitives.addNodeAt(t,i); t->linkVertices(vertices); }
	triangleMesh& smoothExpLaplacian(const float lambda, const unsigned int nit);
	triangleMesh& smoothExpDesbrun(const float lambda, const unsigned int nit);
	triangleMesh& smoothImpDesbrun(const float lambda);
	triangleMesh& smoothExpFujiwara(const float lambda, const unsigned int nit);
	triangleMesh& smoothImpFujiwara(const float lambda);
	triangleMesh& smoothExpTaubin(const float lambda, const unsigned int nit);
	triangleMesh& smooth(const smootherType t, const float lambda, const unsigned int nit = 50);
	triangleMesh& deform(const fimage& src, const float sigma, const float ilim,
			const float wt, const float tol, const bool check, const bool verbose);
	triangleMesh& transform(const fvec3& tr, const fvec3& rt, const fvec3& sc);
	triangleMesh& inflate(const unsigned int nit, const bool verbose = false);
	triangleMesh convexHull(const bool verbose = false) const;
	triangleMesh regNonRigid(const triangleMesh& ref, const unsigned nit = 100) const;
// --- functions changing mesh connectivity
  	void	compact();
	void	linkFully();
	bool	test(const bool verbose);
   	triangleMesh& orient(const bool swap = false);
	void	simplify(const unsigned int nf, const bool check);
	void	simplify(const unsigned int nf, cmdlist* h, triangleMesh& m);
	void	isolateVertex(const unsigned int a);
	void	removeTriangle(const unsigned int i);
	triangleMesh& selectBig();
	triangleMesh& correctManifoldAt(const unsigned int a, const bool verbose);
	unsigned int fixTriangleQuality(const float minArea, const float minAng);
	void	weedTriangles();
// --- functions operating on scalar values
	bool	smoothValues(const float lambda, const unsigned int nit);
	fvecD	computeDivergence(const vfvec3& grad) const;
	fvecD	euclidianDTF(const fvecD& col) const;
	fvecD	findPeaks(const fvecD& col, const bool min, const float lim) const;
	fvecD	geodesicDepth();
	fvecD	geodesicDepth(const bimage& wm) const;
	bvecD	binarize(const fvecD& src, const float vmin, const float vmax) const;
	bvecD	binarize(const uvecD& src, const unsigned int vmin, const unsigned int vmax) const;
	bvecD	erode(const bvecD& src, const unsigned int nit) const;
	bvecD	dilate(const bvecD& src, const unsigned int nit) const;
	uvecD	label(const bvecD& src) const;
	vfvecD	transformLB(const unsigned int dim);					// NB commented out in transform.C
	fmatS	heatKernel(const float gamma) const;
	vfvec3	computeFaceGrad(const fvecD& col) const;
	dmatD	sph(const sphereMesh& sp, const unsigned int lmax) const;
	void	isph(const dmatD& C);
	void	alignWithAxes(const triangleMesh& m);
	void	markDefects(sphereMesh& s, const bool verbose);
};

triangleMesh generateIsosurface(const fimage& src, const float lim);
triangleMesh generateIsosurface(const bimage& src, const connectivity cn);

using vmesh = std::vector<triangleMesh>;

#endif

