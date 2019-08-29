
#ifndef SPHEREMESH_H
#define SPHEREMESH_H

/*
 *
 * triangleMesh.h: Basic definitions and classes for handling spherical meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: sphereMesh.h 438 2016-11-10 00:58:28Z frithjof $
 *
 * v400 (20/09/16): extracted for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions and classes for handling spherical meshes.
*/

//! Represents spherical surface meshes in 3D.

struct sphereMesh : public triangleMesh {
	unsigned int addIfNew(const fvec3& p);
	void	centerVertex(vertex* v)							//! sets vertex position to center of one-ring.
		{ dvec3 d = 0.0;
		  for (const auto& e : v->edges) { dvec3 p = pos(e.id); d += p; };
		  v->pos = d.normalize(); }
	double	TuetteEnergy() const;
	void	MoebiusTransformation(vfvec3& dt);
	void	projectOntoSphere();
	void	subdivideEdges(const unsigned int n);
	fvec3	balanceAreasAt(const unsigned int id) const;
	void	promoteSolution(const triangleMesh& m, const cmdlist& h);
	bool	intersectionAt(const unsigned int id) const;
	bool	edgeMeshIntersection(const unsigned int b0, const unsigned int b1) const;
	double	interpSplineAt(const std::list<unsigned int>& nb, const fvec3& p, const float lambda = 0.0f) const;
	double	pointInTriangle(const dvec3& p, const unsigned int id) const;
	bool	pointInTriangle(const unsigned int i, const dvec3& p) const;
	bool	pointOnEdge(const unsigned int i, const dvec3& p) const;
	std::vector<unsigned int> withinRadius(const unsigned int id, const float rad) const;
	std::vector<unsigned int> withinRadius1(const unsigned int id, const float rad) const;
public:
   	sphereMesh()									//! allocates an empty sphere mesh.
		 : triangleMesh() { }
   	sphereMesh(const vgraph& v, const pgraph& p)					//! allocates a volumetric mesh from a vertex and cell graph.
		: triangleMesh(v,p) { }
	sphereMesh(const triangleMesh& m)						//! copies mesh from m.
		 : triangleMesh(m) { }
   	sphereMesh(const sphereMesh& m)							//! copies from sphere mesh m.
		 : triangleMesh(m) { }
	sphereMesh(sphereMesh&& m)							//! moves mesh from m.
		 : triangleMesh(std::move(m)) { }
	sphereMesh(const unsigned int n);
 	sphereMesh& operator=(const sphereMesh& b)					//! assigns sphere mesh from m.
		{ if (this != &b) triangleMesh::operator=(b); return *this; }
 	sphereMesh& operator=(sphereMesh&& b)						//! move assigns sphere mesh from m.
		{ assert(this != &b); triangleMesh::operator=(std::move(b)); return *this; }
	void	read(const char* fname)							//! reads mesh from file named fname.
		{ mesh::read(fname); }
	void	save(const char* fname)							//! reads mesh from file named fname.
		{ mesh::save(fname); }
	void	read(FILE* fp)								//! reads mesh from file fp.
		{ triangleMesh::read(fp);
		  const char* s = pmAttributes().lookup("graphRepn");
		  if (s && strcmp(s,"sphere")) throw rtException("Not a spherical mesh"); }
	void	save(FILE* fp)								//! reads mesh from file fp.
		{ triangleMesh::save(fp); }
	void	mapTuette(const unsigned int nit = 1000, const bool verbose = false);
	unsigned int mapAreas(const unsigned int nit = 1000, const bool verbose = false);
	void	mapConformal(const bool verbose = false);
	void	mapIsometric(const unsigned int n, const float alpha_c, 
			const float alpha_l, const float alpha_a, const float incnf,
			const enType tp, const bool verbose);
	void	mapRemesh(const bool verbose);
	bool	unfoldToSphere(const bool verbose = false);
	uvecD	detectDefects(const bool verbose = false) const;
	unsigned int vertexAt(const dvec3& p) const;
	float	filterGaussianAt(const unsigned int id, const float sigma, const float rad) const;
	void	filterGaussian(const float sigma, const float rad);
	float	interpolateNNAt(const unsigned int id) const;				//! returns scalar at vertex id.
	float	interpolateAt(const unsigned int id, const dvec3& p) const;		//! returns scalar interpolated at pi.
	float	interpolateLabelAt(const unsigned int id, const dvec3& p) const;
	float	interpolateSplineAt(const unsigned int id, const dvec3& p, const float rad, const float lambda) const;
	unsigned int triangleAt(const dvec3& p, const unsigned int i0) const;
	dvec3	baryWeights(const dvec3& p, const unsigned int id) const;
	sphereMesh transform(const sphereMesh& src, const fmat3& m, const bool label) const;
	sphereMesh deform(const sphereMesh& src, const vfvec3& def, const bool label) const;
};

//! Provides a context for interpolation on spherical surface meshes.

class sphereInterpolator : public kdTree<float>, public sphereMesh {
	float	getData(const unsigned int i, const unsigned int d) const		//! returns scalar field d of data iterm i.
		{ dvec3 p(pos(i)); return d == 0? p.x: d == 1? p.y: p.z; }
	float	distance(const fvecD& v, const unsigned int i) const			//! returns distance of vector v and data iterm i.
		{ const fvec3 p = pos(i); return SQR(p.x-v(0))+SQR(p.y-v(1))+SQR(p.z-v(2)); }
public:
	sphereInterpolator()
		//! Allocates an empty context for nearest neighbor search.
		: kdTree<float>(3), sphereMesh()
		{ }
	sphereInterpolator(const sphereMesh& _m)
		//! Allocates a context for nearest neighbor search in Euclidean space for spherical mesh _m.
		: kdTree<float>(3), sphereMesh(_m)
		{ const unsigned int n = nVertices(); uvecD ind(n);
		  for (unsigned int i = 0; i < n; i++) ind[i] = i;
		  build(ind); }
	sphereInterpolator(const sphereInterpolator& si)				//! copy-constructs interpolator from si.
		: kdTree(si), sphereMesh(si) { }
	sphereInterpolator(sphereInterpolator&& si)					//! move-constructs interpolator from si.
		: kdTree(std::move(si)), sphereMesh(std::move(si)) { }
	sphereInterpolator& operator=(const sphereInterpolator& b)			//! copy-assigns a tree from b.
		{ if (this != &b) { kdTree::operator=(static_cast<const kdTree&>(b));
		  	sphereMesh::operator=(static_cast<const sphereMesh&>(b)); };
		  return *this; }
	sphereInterpolator& operator=(sphereInterpolator&& b)				//! move-assigns a tree from b.
		{ assert(this != &b); 
		  kdTree::operator=(std::move(static_cast<kdTree&>(b)));
		  sphereMesh::operator=(std::move(static_cast<sphereMesh&>(b)));
		  return *this; }
	sphereMesh downsample(const float fac);						//! returns a spherical mesh downsampled by factor fac.
	unsigned int vertexAt(const fvec3& p) const					//! returns vertex id of nearest neighbor to point p.
		{ fvecD q(3); q[0] = p.x; q[1] = p.y; q[2] = p.z;
		  return nearest(q); }
	unsigned int triangleAt(const fvec3& p) const					//! returns triangle id of that contains point p.
		{ const unsigned int id = vertexAt(p); assert(id);
		  return sphereMesh::triangleAt(p,id); }
	std::vector<unsigned int> within(const fvec3& p, const float rad) const 	//! returns list of vertices within radius rad of point p.
		{ fvecD q(3); q[0] = p.x; q[1] = p.y; q[2] = p.z;
		  return kdTree::within(q,rad); }
	float	interpolateNNAt(const fvec3& p) const					//! returns scalar at point p using nearest neighbor interpolation.
		{ const unsigned int id = vertexAt(p);
		  return sphereMesh::interpolateNNAt(id); }
	float	interpolateAt(const fvec3& p) const					//! returns scalar at point p using linear interpolation.
		{ const unsigned int id = vertexAt(p);
		  return sphereMesh::interpolateAt(id,p); }
	float	interpolateLabelAt(const fvec3& p) const				//! returns scalar at point p using label interpolation.
		{ const unsigned int id = vertexAt(p);
		  return sphereMesh::interpolateLabelAt(id,p); }
	float	interpolateSplineAt(const fvec3& p, const float rad,
			const float lambda = 1e-8) const				//! returns scalar at point p using spline interpolation.
		{ const unsigned int id = vertexAt(p);
		  return sphereMesh::interpolateSplineAt(id,p,rad,lambda); }
	float	similaritySSD(const sphereMesh& ref, const fquat& q, const bool label) const;
	float	similarityXOR(const sphereMesh& ref, const fquat& q) const;
	float	similarityNMI(const sphereMesh& ref, const fquat& q, const bool label) const;
	float	similarityCC(const sphereMesh& ref, const fquat& q, const bool label) const;
	float	similarityCCRE(const sphereMesh& ref, const fquat& q, const bool label) const;
	fvec3	rotate(const fvec3& a, const float r, const fvec3& v) const;
	fvec3	rotate(const fvec3& s, const fvec3& d, const fvec3& t) const;
	void	smoothDeformationField(vfvec3& def, const float lambda,
			const unsigned int nit) const;
	fvec3	interpolateAt(const fvec3& p, const vfvec3& fld) const;
	void	deformField(vfvec3& fld, const vfvec3& def) const;
	fquat	linearRegistration(const sphereMesh& ref) const;
};

#endif

