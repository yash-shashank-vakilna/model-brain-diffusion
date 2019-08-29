#ifndef MESH_H
#define MESH_H

/*
 *
 * mesh.h: Basic definitions and classes for handling meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: mesh.h 451 2016-12-05 01:29:53Z frithjof $
 *
 * 0.10 (17/02/99): initial release
 * 0.20 (12/09/09): adapted for BRIAN 1.20 FK
 * 0.21 (28/03/11): renamed faces to triangles & added cells
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v400 (20/09/16): rewritten for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions and classes for handling meshes.
*/

using vgraph = graph<vertex>;
using pgraph = graph<primitive>;

#undef allVertices
#undef allPrimitives
#undef allTrianglesAt

#define allVertices(i)	(unsigned int i = 0; i < nVertices(); i++)
#define allPrimitives(i) (unsigned int i = 0; i < nPrimitives(); i++)
#define allTrianglesAt(id) (auto it = vtop.find(id); it->first == id; it++)

//! Collects information about a mesh.

struct mesh {
   	vgraph	vertices;			//!< contains vertices of a mesh
   	pgraph	primitives;			//!< contains geometric primitives of a mesh
	mmap	vtop;				//!< maps vertices to primitives
protected:
	bool	testPrimitiveVertices(const bool verbose) const;
	void	linkPrimitiveNeighbours();
public:
	mesh()										//! allocates an empty mesh.
		 : vertices(), primitives(), vtop(mmap()) { }
   	mesh(const vgraph& v, const pgraph& p)						//! allocates a mesh from a vertex and primitive graph.
		 : vertices(v), primitives(p), vtop(mmap()) { }
   	mesh(const vgraph& v)								//! allocates a mesh from a vertex graph.
		 : vertices(v), primitives(), vtop() { }
	mesh(const mesh& m)								//! copies mesh from m.
		 : vertices(m.vertices), primitives(m.primitives), vtop(m.vtop) { }
	mesh(mesh&& m)									//! moves mesh from m.
		 : vertices(std::move(m.vertices)),
		   primitives(std::move(m.primitives)), 
		   vtop(std::move(m.vtop)) { }
  	virtual ~mesh();
 	mesh& operator=(const mesh& m)							//! assigns mesh from m.
		{ if (this != &m) { vertices = m.vertices; primitives = m.primitives; };
		  vtop = m.vtop; return *this; }
 	mesh& operator=(mesh&& m)							//! move assigns mesh from m.
		{ assert(this != &m);
		  vertices = std::move(m.vertices); primitives = std::move(m.primitives);
		  vtop = std::move(m.vtop); return *this; }
	void	clear()									//! remove all nodes and primitives. 
		{ vertices.clear(); primitives.clear(); vtop.clear(); }
	void	read(const char* fname)							//! reads mesh from file named fname.
		{ FILE* fp = openFile(fname,"r"); read(fp); closeFile(fp); }
	void	save(const char* fname)							//! saves mesh to file named fname.
		{ FILE* fp = openFile(fname,"w"); save(fp); closeFile(fp); }
	attrList& vtxAttributes()							//! returns vertex attribute list.
		{ return vertices.attributes(); }
	attrList& pmAttributes()							//! returns primitive attribute list.
		{ return primitives.attributes(); }
	const attrList& vtxAttributes()	const						//! returns vertex attribute list.
		{ return vertices.attributes(); }
	const attrList& pmAttributes() const 						//! returns primitive attribute list.
		{ return primitives.attributes(); }
	unsigned int nVertices() const							//! returns number of vertices in mesh.
		{ return vertices.size(); }
	unsigned int nPrimitives() const						//! returns number of primitives in mesh.
		{ return primitives.size(); }
	bool	isCompact() const 							//! checks if mesh is compact.
		{ return vertices.isCompact() && primitives.isCompact(); }
	const fvec3 pos(const unsigned int i) const					//! returns vertex position i.
		{ return vertices(i)->pos; }
	fvec3&	pos(const unsigned int i)						//! returns reference to vertex position i.
		{ return vertices(i)->pos; }
	float	distance(const unsigned int i, const unsigned int j) const		//! returns distance between vertex positions (i,j).
		{ return norm(pos(i)-pos(j)); }
	void	copyAttributes(const mesh& m)						//! copies attributes from mesh m.
		{ vertices.copyAttributes(m.vtxAttributes());
		  primitives.copyAttributes(m.pmAttributes()); }
	void	addVertex(vertex* v, const unsigned int i)				//! adds vertex at i
		{ vertices.addNodeAt(v,i); }
	unsigned int nComponents() const;
	unsigned int nEdges() const;
	float	edgeLength() const;
	void	edgeLengths();
	float	shortestEdge() const;
	float	longestEdge() const;
	void	boundingBox(fvec3& bmin, fvec3& bmax) const;
	void	range(float& vmin, float& vmax) const;
	void	setValues(const fimage& im);
	void	setValues(const fvecD& val);
	void	setLabels(const uvecD& lbl);
	void	setValues(const bvecD& ind);
	fvecD	getValues() const;
	uvecD	getLabels() const;
	void	vertexToPrimitiveMap();
	unsigned int countPrimitivesAtEdge(const unsigned int i, const unsigned int j) const;
	void	linkPrimitivesAtEdge(const unsigned int i, const unsigned int j);
	void	unlinkPrimitive(const unsigned int i, const unsigned int t);
	std::vector<unsigned int> findPrimitivesAtEdge(const unsigned int i, const unsigned int j) const;
	unsigned int vertexAt(const fvec3& p) const;
	std::list<unsigned int> withinRadius(const unsigned int id, const float rad) const;
	std::list<unsigned int> nRing(const unsigned int id, const unsigned int n) const;
	virtual void read(FILE* fp) = 0;
	virtual void save(FILE* fp) = 0;
  	virtual void compact() = 0;
  	virtual void linkFully() = 0;
  	virtual vfvec3 computeNormals() const = 0;
};
#endif
