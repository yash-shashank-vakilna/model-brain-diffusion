#ifndef VERTEX_H
#define VERTEX_H

/*
 *
 * vertex.h: nodes in a vertex graph
 * BRIAN Software Package Version 3.0
 *
 * $Id: vertex.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (02/04/11): rewritten for BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v400 (20/09/16): rewritten for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for vertex variants.
 */

//! Implements a vertex in 3D

struct vertex : public node {								// former type 1 vertex
	fvec3	pos;				//!< vertex position

	vertex(const fvec3& p = 0.0f)							//! allocates a vertex with position p.
		: node(), pos(p) { }
	vertex(const float x, const float y, const float z)				//! allocates a vertex with position (x,y,z).
		: node(), pos(x,y,z) { }
	vertex(const vertex& v)
		: node(v), pos(v.pos) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertex from stream in.
		{ node::read(in,hasWeights); pos.read(in); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertex to stream out.
		{ node::save(out,hasWeights); pos.save(out); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z);
		  node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertex; }
	vertex* clone() const								//! clones this instance.
		{ return new vertex(*this); }
	virtual float getValue() const
		{ return 0.0f; }
	virtual void setValue(const float)
		{ }
	virtual fvec3 getNormal() const
		{ return fvec3(0.0f); }
	virtual fvec3& normal()
		{ throw rtException("No normal in this vertex"); }
	virtual void setNormal(const fvec3&)
		{ }
	virtual fvec3 getColor() const
		{ return fvec3(0.0f); }
	virtual void setColor(const fvec3&)
		{ throw rtException("No color in this vertex"); }
};

struct vertexV : public vertex {							// former type 4 vertex
	float	val;				//!< scalar

	vertexV(const fvec3& p = 0.0f, const float v = 0.0f)				//! allocates a vertex with position p and scalar v.
		: vertex(p), val(v) { }
	vertexV(const vertex& v)
		: vertex(v), val(0.0f) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertexV from stream in.
		{ vertex::read(in,hasWeights); in.read(val); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertexV to stream out.
		{ vertex::save(out,hasWeights); out.save(val); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z); 
		  printf("[%8.4f] : ", val); node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertexV; }
	vertexV* clone() const								//! clones this instance.
		{ return new vertexV(*this); }
	float	getValue() const
		{ return val; }
	void	setValue(const float v)
		{ val = v; }
};

struct vertexN : public vertex {							// former type 2 vertex
	fvec3	norm;				//!< vertex normal

	vertexN(const fvec3& p = 0.0f, const fvec3& n = 0.0f)				//! allocates a vertex with position p and normal n.
		: vertex(p), norm(n) { }
	vertexN(const vertex& v)
		: vertex(v), norm(0.0f) { }
	vertexN(const vertexN& v)
		: vertex(v), norm(v.norm) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertexN from stream in.
		{ vertex::read(in,hasWeights); norm.read(in); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertexN to stream out.
		{ vertex::save(out,hasWeights); norm.save(out); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z); 
		  printf("(%8.4f %8.4f %8.4f) : ", norm.x,norm.y,norm.z); 
		  node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertexN; }
	vertexN* clone() const								//! clones this instance.
		{ return new vertexN(*this); }
	fvec3	getNormal() const
		{ return norm; }
	fvec3&	normal()
		{ return norm; }
	void	setNormal(const fvec3& n)
		{ norm = n; }
};

struct vertexC : public vertex {							// former type 3 vertex
	fvec3	col;				//!< RGB color

	vertexC(const fvec3& p = 0.0f, const fvec3& c = 0.0f)				//! allocates a vertex with position p and RGB color c.
		 : vertex(p), col(c) { }
	vertexC(const vertex& v)
		: vertex(v), col(0.0f) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertexC from stream in.
		{ vertex::read(in,hasWeights); col.read(in); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertexC to stream out.
		{ vertex::save(out,hasWeights); col.save(out); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z); 
		  printf("(%6.4f %6.4f %6.4f) : ", col.x,col.y,col.z); 
		  node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertexC; }
	vertexC* clone() const								//! clones this instance.
		{ return new vertexC(*this); }
	fvec3	getColor() const
		{ return col; }
	void	setColor(const fvec3& c)
		{ col = c; }
};

struct vertexNV : public vertexN {							// former type 5 vertex
	float	val;				//!< scalar

	vertexNV(const fvec3& p = 0.0f, const fvec3& n = 0.0f, const float v = 0.0f)	//! allocates a vertex with position p, normal n and scalar v.
		 : vertexN(p,n), val(v) { }
	vertexNV(const vertexN& v)
		: vertexN(v), val(0.0f) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertexN from stream in.
		{ vertexN::read(in,hasWeights); in.read(val); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertexN to stream out.
		{ vertexN::save(out,hasWeights); out.save(val); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z); 
		  printf("(%8.4f %8.4f %8.4f), ", norm.x,norm.y,norm.z); 
		  printf("[%8.4f] : ", val); node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertexNV; }
	vertexNV* clone() const								//! clones this instance.
		{ return new vertexNV(*this); }
	float	getValue() const
		{ return val; }
	void	setValue(const float v)
		{ val = v; }
};

struct vertexNC : public vertexN {							// former type 6 vertex
	fvec3	col;				//!< RGB color

	vertexNC(const fvec3& p = 0.0f, const fvec3& n = 0.0f, const fvec3& c = 0.0f)	//! allocates a vertex with position p, normal n and RGB color c.
		 : vertexN(p,n), col(c) { }
	vertexNC(const vertexN& v)
		: vertexN(v), col(0.0f) { }
	void	read(is& in, const bool hasWeights = false)				//! reads a vertexN from stream in.
		{ vertexN::read(in,hasWeights); col.read(in); }
	void	save(os& out, const bool hasWeights = false) const			//! saves a vertexN to stream out.
		{ vertexN::save(out,hasWeights); col.save(out); }
	void	print(const bool hasWeights = false) const
		{ printf("(%8.4f %8.4f %8.4f), ", pos.x,pos.y,pos.z); 
		  printf("(%8.4f %8.4f %8.4f), ", norm.x,norm.y,norm.z); 
		  printf("(%6.4f %6.4f %6.4f) : ", col.x,col.y,col.z); 
		  node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::vertexNC; }
	vertexNC* clone() const								//! clones this instance.
		{ return new vertexNC(*this); }
	fvec3	getColor() const
		{ return col; }
	void	setColor(const fvec3& c)
		{ col = c; }
};
#endif
