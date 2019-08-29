#ifndef PRIMITIVE_H
#define PRIMITIVE_H

/*
 *
 * primitive.h: spatial structures in a primitive graph
 * BRIAN Software Package Version 3.0
 *
 * $Id: primitive.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (02/04/11): rewritten for BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v400 (20/09/16): rewritten for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions for spatial structures in a primitive graph.
*/

using mmap = std::multimap<unsigned int,unsigned int>;

#define allVtx(i) (unsigned int i = 0; i < vtx.N; i++)

//! Base class for all primitives.

struct primitive : public node {
   	uvecD	vtx;				//!< vertex indices

	primitive(const unsigned int n = 0)						//! allocates a primitive with n vertices
		: node(), vtx(n) { }
   	primitive(const primitive& p)							//! copies from primitive p.
		: node(p), vtx(p.vtx) { }
   	primitive(primitive&& p)							//! moves from primitive p.
		: node(std::move(p)), vtx(std::move(p.vtx)) { }
	primitive& operator=(const primitive& p)					//! assigns from  primitive p including links.
		{ if (this != &p) { node::operator=(p); vtx = p.vtx; }; return *this; }
	primitive& operator=(primitive&& p)						//! move assigns from  primitive p including links.
		{ assert(this != &p); node::operator=(p);
		  vtx = std::move(p.vtx); return *this; }
        bool	isDegenerate() const  							//! checks a primitive for duplicated vertices.
		{ for allVtx(i) 
			for (unsigned int j = i+1; j < vtx.N; j++)
				if (vtx(i) == vtx(j)) return true;
		  return false; }
	template<typename N> bool checkVertices(const graph<N>& g) const 		//! checks a primitive for missing vertices.
		{ for allVtx(i) { if (g(vtx(i)) == nullptr) return false; };
		  return true; }
	void	removeVertexToprimitiveEntries(mmap& m, const unsigned int id)		//! removes references from all vertices to a primitive.
		{ for allVtx(i) {
			for (auto it = m.find(i); it != m.end(); it++) {
				if (it->first != i) break;
				if (it->second == id) { m.erase(it); break; } } } }
	bool	hasVertex(const unsigned int v) const					//! returns true if triangle contains vertex v.
		{ for allVtx(i) { if (vtx(i) == v) return true; }; return false; }
  	template<typename N> N* getVertex(const graph<N>& g, const unsigned int i) const //! returns vertex i.
		{ assert(i < vtx.N); const unsigned int id = vtx(i); return g(id); }
   	template<typename N> fvec3 pos(const graph<N>& g, const unsigned int i) const	//! returns position of vertex i.
		{ N* v = getVertex(g,i); assert(v); return v->pos; }
        bool	vertexPosition(unsigned int& i0, const unsigned int id) const		//! given vertex id, returns vertex array index.
		{ for allVtx(i) { if (vtx(i) == id) { i0 = i; return true; } };
		  return false; }
	bool	replaceVertex(const unsigned int a, const unsigned int b)		//! replaces vertex a by b.
		{ unsigned int i; const bool succ = vertexPosition(i,a);
		  if (succ) { vtx[i] = b; }; return succ; }
	void	read(is& in, const bool hasWeights = false)				//! reads a primitive from stream in.
		{ node::read(in,hasWeights); unsigned int n = 0; in.read(n);
		  vtx.resize(n); for allVtx(i) { in.read(n); vtx[i] = n; } }
	void	save(os& out, const bool hasWeights = false) const			//! saves a primitive to stream out.
		{ node::save(out,hasWeights); unsigned int n = vtx.N; out.save(n);
		  for allVtx(i) { n = vtx(i); out.save(n); } }
	void	print(const bool hasWeights = false) const
		{ printf("("); for (unsigned int i = 0; i < vtx.N-1; i++)
			printf("%6d ", vtx(i));
		  printf("%6d) : ", vtx(vtx.N-1)); node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::primitive; }
	void	addVertices(mmap& m, const unsigned int pid) const			//! adds vertices of this primitive to mmap
		{ for allVtx(i) { bool succ = false; const unsigned int vid = vtx(i);
			for (auto it = m.find(vid); it != m.end() && it->first == vid; it++)
				if (it->second == pid) { succ = true; break; };
			if (succ == false) m.insert({vid,pid}); } }
	primitive* clone() const							//! clones this instance.
		{ return new primitive(*this); }
};

//! Implements a triangle node

struct triangle : public primitive {
   	triangle()
		: primitive(3) { }
   	triangle(const unsigned int a, const unsigned int b, const unsigned int c)	//! allocates a triangle node with vertex ids (a,b,c).
		: primitive(3) { vtx[0] = a; vtx[1] = b; vtx[2] = c; }
	void	testOrientation(const unsigned int v0, const unsigned int v1)		//! checks orientation of triangle w.r.t. edge (v0,v1).
		{ unsigned int i0,i1;
		  if ((vertexPosition(i0,v0) & vertexPosition(i1,v1)) == false) return; 
		  if ((i0 > i1 && i0-i1 == 1) || (i0 < i1 && i1-i0 == 2)) return;
		  std::swap(vtx[i0],vtx[i1]); }
	template<typename N> fvec3 computeNormal(const graph<N>& g) const		//! returns plane normal.
		{ const fvec3 p0 = pos(g,0), e1 = pos(g,1)-p0, e2 = pos(g,2)-p0;
		  return cross(e1,e2).normalize(); }
	template<typename N> fvec3 center(const graph<N>& g) const			//! returns barycenter of triangle.
		{ return (pos(g,0)+pos(g,1)+pos(g,2))/3.0f; }
  	bool	hasEdge(const unsigned int v0, const unsigned int v1) const		//! returns true if triangle contains edge (v0,v1).
		{ return hasVertex(v0) && hasVertex(v1); }
	template<typename N> double cornerAngle(const graph<N>& g, const unsigned int i) const
		//! returns angle at vertex i.
		{ unsigned int ip = (i == 0)? 2: i-1, in = (i == 2)? 0: i+1;
		  const dvec3 p0(pos(g,i)), p1(pos(g,ip)), p2(pos(g,in));
		  const double d = dot((p1-p0).normalize(),(p2-p0).normalize());
		  return acos(CLAMP(d,-1.0,1.0)); }
	unsigned int otherVertex(const unsigned int v0, const unsigned int v1) const	//! given two vertices v0, v1, finds the other vertex.
		{ for allVtx(i) { if (vtx(i) != v0 && vtx(i) != v1) return vtx(i); }
		  return 0; }
	template<typename N> double area(const graph<N>& g) const			//! returns triangle area.
		{ dvec3 p0(pos(g,0)), p1(pos(g,1)), p2(pos(g,2));
		  return 0.5*norm(cross(p1-p0,p2-p0)); }
        bool	isDegenerate() const  							//! checks a triangle for duplicated vertices.
		{ return vtx(0) == vtx(1) || vtx(0) == vtx(2) || vtx(1) == vtx(2); }
	bool	sharesEdge(const triangle& t) const  					//! checks if triangle shares an edge with t.
		{ for allVtx(i) { if (t.hasEdge(vtx(i),vtx((i+1)%vtx.N))) return true; };
		  return false; }
	template<typename N> bool check(const graph<N>& g, const float area_l, const float ang_l) const
		//! check if triangle is small or obtuse.
		{ for allVtx(i) { if (cornerAngle(g,i) < ang_l) return false; };
		  return area(g) > area_l; } 
	template<typename N> void bound(const graph<N>& g, fvec3& bmin, fvec3& bmax) const
		//! checks bounds of triangle.
		{ for allVtx(i) { fvec3 p = pos(g,i); bmin = bmin.min(p); bmax = bmax.max(p); } }
	template<typename N> bool outwardNormal(const graph<N>& g) const
		//! returns true if position vector and normal point into the same half-sphere.
		{ const fvec3 p = pos(g,0), n = computeNormal(g); return dot(n,p) > 0.0f; }
	template<typename N> void linkVertices(graph<N>& g) const 			//! links vertices in a triangle.
		{ g.linkNodesBi(vtx(0),vtx(1)); g.linkNodesBi(vtx(0),vtx(2));
		  g.linkNodesBi(vtx(1),vtx(2)); }
	float	intersectsRay(const fvec3& org, const fvec3& dir) const			//! returns true if ray (org,dir) intersects with this triangle.
		{ const float m = std::numeric_limits<float>::max();
		  const fvec3 e1 = vtx(1)-vtx(0), e2 = vtx(2)-vtx(0), h = cross(dir,e2);
		  const float a = dot(e1,h); if (std::abs(a) < 0.00001f) return m;
		  const fvec3 s = org-vtx(0); const float f = 1.0f/a, u = f*dot(s,h);
		  if (u < 0.0 || u > 1.0) return m;
		  const fvec3 q = cross(s,e1); const float v = f*dot(dir,q);
		  if (v < 0.0 || u+v > 1.0) return m;
		  const float t = f*dot(e2,q); return t > 0.00001f? t: m; }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::triangle; }
	triangle* clone() const								//! clones this instance.
		{ return new triangle(*this); }
};


/*! \enum tesselationType
    \brief Symbolic constants for cell types in volumetric meshes.
*/

enum class tesselationType {
	tetra5,					///< tetrahedron from 5-division of 8-cell
	tetra6a,				///< tetrahedron from 6-division of 8-cell, type a
	tetra6b,				///< tetrahedron from 6-division of 8-cell, type b
	hexa					///< hexahedron
};

/*! \enum materialType
    \brief Symbolic constants for material types in volumetric meshes.
*/

enum class materialType {
	undef = 0,				///< undefined material (background,air)
	whiteMatter = 1,			///< white matter
	muscle = 2,				///< skeletal muscle
	meninges = 3,				///< meninges
	greyMatter = 4,				///< grey matter
	skin = 5,				///< skin
	cartilage = 6,				///< cartilage
	csf = 7,				///< cerebro-spinal fluid
	fluid = 8,				///< fluid material (modeled as fluid)
	bone = 9,				///< bone
	fat = 10,				///< soft tissue fat
};

//! Collects mechanical properties of a material

struct materialData {
	const materialType id;			//!< material id
	const char *name;			//!< description
	const float rho;			//!< density
	const float nu;				//!< Poisson's ratio (viscosity for fluids)
	const float E;				//!< Young's modulus
	const float kappa;			//!< bulk modulus
};

extern const materialData& lookupMaterial(const materialType id);

//! Implements a volumetric node (aka cell)

struct cell : public primitive {
	materialType mat;			//!< material code

	cell(const unsigned int n = 0)							//! allocates an empty cell.
		: primitive(n), mat(materialType::undef) { }
   	cell(const cell& c)								//! copies from cell c.
		: primitive(c), mat(c.mat) { }
   	cell(cell&& c)									//! moves from cell c.
		: primitive(std::move(c)), mat(c.mat) { }
	cell& operator=(const cell& c)							//! assigns from cell c including links.
		{ if (this != &c) { primitive::operator=(c); mat = c.mat; };
		  return *this; }
	cell& operator=(cell&& c)							//! move assigns from node b including links.
		{ assert(this != &c); primitive::operator=(c); mat = c.mat;
		  return *this; }
	bool	hasFace(const unsigned int a, const unsigned int b, const unsigned int c) const
		//! returns true if cell shares face (a,b,c).
		{ return hasVertex(a) && hasVertex(b) && hasVertex(c); } 
	void	read(is& in, const bool hasWeights = false)				//! reads a cell from stream in.
		{ primitive::read(in,hasWeights); in.read(mat); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a cell to stream out.
		{ primitive::save(out,hasWeights); out.save(mat); }
	void	print(const bool hasWeights = false) const
		{ printf("("); for (unsigned int i = 0; i < vtx.N-1; i++)
			printf("%6d ", vtx(i));
		  printf("%6d), [%u] : ", vtx(vtx.N-1), static_cast<unsigned int>(mat));
		  node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::cell; }
	cell*	clone() const								//! clones this instance.
		{ return new cell(*this); }
	bool	isSolid() const
		{ return mat != materialType::fluid; }
	const unsigned int* edgeList(unsigned int& ec) const;
	const unsigned int* faceList(unsigned int& fc) const;
	float 	getPossionRatio() const
		{ return lookupMaterial(mat).nu; }
	float	getYoungModulus() const
		{ return lookupMaterial(mat).E; }
	float	getDensity() const
		{ return lookupMaterial(mat).rho; }
	float	getViscosity() const
		{ return lookupMaterial(mat).nu; }
	const char* getName() const
		{ return lookupMaterial(mat).name; }
};

//! Implements a 4-node tetrahedron - specialized from cell

struct tetrahedron : public cell {
	tetrahedron() : cell(4) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a tetrahedron from stream in.
		{ cell::read(in,hasWeights); in.read(mat); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a tetrahedron to stream out.
		{ cell::save(out,hasWeights); out.save(mat); }
	void	print(const bool hasWeights = false) const
		{ cell::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::tetrahedron; }
	tetrahedron* clone() const							//! clones this instance.
		{ return new tetrahedron(*this); }
};

//! Implements a 10-node tetrahedron - specialized from cell

struct tetrahedron10 : public cell {
	tetrahedron10() : cell(10) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a tetrahedron from stream in.
		{ cell::read(in,hasWeights); in.read(mat); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a tetrahedron to stream out.
		{ cell::save(out,hasWeights); out.save(mat); }
	void	print(const bool hasWeights = false) const
		{ cell::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::tetrahedron10; }
	tetrahedron10* clone() const							//! clones this instance.
		{ return new tetrahedron10(*this); }
};

//! Implements a 6-node hexahedron - specialized from cell

struct hexahedron : public cell {
	hexahedron() : cell(8) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a hexahedron from stream in.
		{ cell::read(in,hasWeights); in.read(mat); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a hexahedron to stream out.
		{ cell::save(out,hasWeights); out.save(mat); }
	void	print(const bool hasWeights = false) const
		{ cell::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::hexahedron; }
	hexahedron* clone() const							//! clones this instance.
		{ return new hexahedron(*this); }
};

//! Implements a 20-node hexahedron - specialized from cell

struct hexahedron20 : public cell {
	hexahedron20() : cell(20) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a hexahedron from stream in.
		{ cell::read(in,hasWeights); in.read(mat); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a hexahedron to stream out.
		{ cell::save(out,hasWeights); out.save(mat); }
	void	print(const bool hasWeights = false) const
		{ cell::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::hexahedron20; }
	hexahedron20* clone() const							//! clones this instance.
		{ return new hexahedron20(*this); }
};

#undef allVtx

#endif
