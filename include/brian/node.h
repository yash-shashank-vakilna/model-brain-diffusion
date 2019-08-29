#ifndef NODE_H
#define NODE_H

/*
 *
 * node.h: classes and definitions for nodes in a graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: node.h 506 2017-03-19 01:39:21Z frithjof $
 *
 * 0.10 (02/04/11): rewritten for BRIAN2
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v400 (20/09/16): rewritten for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

enum class nodeID {
	node,
	vertex,
	vertexV,
	vertexN,
	vertexC,
	vertexNV,
	vertexNC,
	primitive,
	triangle,
	cell,
	tetrahedron,
	tetrahedron10,
	hexahedron,
	hexahedron20,
	descriptor,
	stream,
	tensor,
	spharm,
	fiberMixture,
	fiberZCD };

//! Collects information about a edge between two nodes

struct edge {
	unsigned int id;			//!< node reference
	float	weight;				//!< weight of this node

	edge(const unsigned int i = 0, const float w = 0.0f)				//! allocates a edge to node i with weight w.
		 : id(i), weight(w) { }
	bool	operator==(const edge& b) const						//! equal if both ids equal.
		{ return id == b.id; }
	void	print(const bool hasWeights) const					//! prints information about this edge
		{ printT(id); if (hasWeights) printT(weight); }
};

//! Collects generic information about a node in a graph

struct node {
	float	weight;				//!< weight of this node
	std::list<edge> edges;			//!< the list of edges

	node(const float w = 0.0f)							//! allocate a node with weight w.
		 : weight(w), edges() { }
	node(const node& n)
		: weight(n.weight), edges(n.edges) { }
	node(node&& n)
		: weight(n.weight), edges(n.edges) { }
	virtual ~node()	{ }								//! base class deconstructor.
	node&	operator=(const node& n)						//! assigns from bundle b.
		{ if (this != &n) { clear(); weight = n.weight; edges = n.edges; };
		  return *this; }
	node&	operator=(node&& n)							//! move assigns from bundle b.
		{ assert(this != &n); weight = n.weight; edges = std::move(n.edges);
		  return *this; }
	void 	clear()									//! clears all data in this edge.
		{ weight = 0.0f; edges.clear(); }
	void	addEdge(const edge& e)							//! adds edge to this node.
		{ edges.push_back(e); }	
	bool	hasEdge(const unsigned int id) const					//! returns true if this node is linked to node id.
		{ for (const auto& e : edges) if (e.id == id) return true;
		  return false; }
	void	unlink(const unsigned int id)						//! unlinks this node from node id.
		{ for (auto it = edges.begin(); it != edges.end(); it++)
			if (it->id == id) { edges.erase(it); return; } }
	unsigned int nEdges() const							//! returns number of edges.
		{ return edges.size(); }
	unsigned int firstNeighbor() const						//! returns id of first neighbor.
		{ return edges.size()? edges.front().id: ~0u; }
	void	clearWeights()								//! clears all weights.
		{ weight = 0.0f; for (auto& e : edges) e.weight = 0.0f; }
	float	totalDegree() const							//! returns weight of this node including all edge weights.
		{ float w = weight; for (auto& e : edges) { w += e.weight; }; return w; }
	void	addWeight(const unsigned int id, const float w)				//! adds weight w to edge id.
		{ for (auto& e : edges) if (e.id == id) { e.weight += w; return; } }
	virtual void read(is& in, const bool hasWeights = false)			//! reads node and edges from byte stream in.
		{ unsigned int ne = 0; in.read(ne);					// read number of edges
		  while (ne--) { unsigned int id = 0; float wt = 0.0f; 			// for all edges
			in.read(id); if (hasWeights) in.read(wt);			// read vertex id and weight
			edges.push_back({id,wt}); };					// save edge
		  if (hasWeights) in.read(weight); else weight = 0.0f; }		// read node weight
	virtual void save(os& out, const bool hasWeights = false) const			//! saves node and edges into byte stream out.
		{ const unsigned int ne = edges.size(); out.save(ne);			// save number of edges
		  for (const auto& e : edges) { 					// for all edges
			out.save(e.id); if (hasWeights) out.save(e.weight); };		// save vertex id and weight
		  if (hasWeights) out.save(weight); }					// save node weight
	virtual void print(const bool hasWeights = false) const				//! prints information about this node.
 		{ if (hasWeights) printT(weight);
		  for (const auto& e : edges) e.print(hasWeights);
		  printf("\n"); }
	virtual nodeID getID() const							//! identifies this class for serialization.
		{ return nodeID::node; }
	virtual node* clone() const							//! clones this instance.
		{ return new node(*this); }
};

template<typename T> node* createNode() { return new T; }				// templated function that creates nodes
using nodeCreator = node*(*)();								// defines a pointer to this function
using nodeCreatorMap = std::map<nodeID,nodeCreator>;

#endif
