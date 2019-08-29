#ifndef VMESH_H
#define VMESH_H

/*
 *
 * volumeMesh.h: Basic definitions and classes for handling volumetric meshes
 * BRIAN Software Package Version 3.0
 *
 * $Id: volumeMesh.h 488 2017-02-16 00:40:01Z kruggel $
 *
 * v400 (20/09/16): extracted for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions and classes for handling volumetric meshes.
*/

//! Represents volumetric meshes in 3D.

struct volumeMesh : public mesh {
	void	linkCellEdges(const bool verbose);
	bool	sharesTriangle(const cell* c, const unsigned int p, const unsigned int q,
			const unsigned int r) const;
	void	walkStrip(const cell* c, const unsigned int* pt, const unsigned int n,
			vfvec3& nm) const;
public:
   	volumeMesh()									//! allocates an empty volumetric mesh.
		: mesh() { }
   	volumeMesh(const vgraph& v, const pgraph& p)					//! allocates a volumetric mesh from a vertex and cell graph.
		: mesh(v,p) { }
	volumeMesh(const volumeMesh& m)							//! copy-constructs volumetric mesh from m.
		: mesh(m) { }
	volumeMesh(volumeMesh&& m)							//! move-constructs volumetric mesh from m.
		: mesh(m) { }
   	volumeMesh(FILE* fp)								//! allocates a volumetric mesh from FILE* fp.
		: mesh() { read(fp); }
 	volumeMesh& operator=(const volumeMesh& m)					//! assigns volume mesh from m.
		{ if (this != &m) mesh::operator=(m); return *this; };
 	volumeMesh& operator=(volumeMesh&& m)						//! move assigns volume mesh from m.
		{ assert(this != &m); mesh::operator=(std::move(m)); return *this; }
	attrList& clAttributes()							//! returns primitive attribute list.
		{ return primitives.attributes(); }
	unsigned int nVertices() const							//! returns number of vertices in mesh.
		{ return vertices.size(); }
	void	read(FILE *fp);
	void	save(FILE *fp);
	void	read(const char* fname)							//! reads mesh from file named fname.
		{ mesh::read(fname); }
	void	save(const char* fname)							//! reads mesh from file named fname.
		{ mesh::save(fname); }
	cell*	cl(const unsigned int i) const						//! returns cell at i.
		{ return reinterpret_cast<cell*>(primitives(i)); }
	void	addCell(cell* c, const unsigned int i)					//! adds cell c at i.
		{ primitives.addNodeAt(c,i); }
	bool	generateGrid(limage& src, const tesselationType t,
			const unsigned int cmin, const unsigned int cmax, const unsigned int order);
  	void	compact();
   	vfvec3	computeNormals() const;
   	void	linkFully();
	std::vector<unsigned int> verticesAt(const unsigned int id) const;
	unsigned int cellAt(const fvec3& p) const;
	bool	interpolateAt(const fvec3& p, float& v) const;
};
#endif
