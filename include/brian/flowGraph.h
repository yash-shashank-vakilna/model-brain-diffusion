#ifndef FLOWGRAPH_H
#define FLOWGRAPH_H

/*
 *
 * flowGraph.h: optimization using max flow graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: flowGraph.h 459 2017-01-15 02:16:35Z frithjof $
 *
 * 0.10 (21/06/15): initial version
 * v406 (28/09/16): bumped to version 3.0
 *
 * Implements the algorithm in:
 * Boykov Y, Kolmogorov V (2004)
 * An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision.
 * IEEE T PAMI 26, 1124-1137.
 *
 * The energy class is based on:
 * Kolmogorov V., Zabih R. (2004)
 * What Energy Functions can be Minimized via Graph Cuts?
 * IEEE T PAMI 26, 147-159. 
 *
 * This code is based on an implementation by Vladimir Kolmogorov but was completely re-written.
 * The original code is copyrighted by Vladimir Kolmogorov and Yuri Boykov under a GPL license.
 *
 */

/*! \file
    \brief Templated max flow graph and energy optimization classes.
*/


using ulist = std::list<unsigned int>;

//! Provides an implementation for a max flow / min cut graph.

template <class T> class flowGraph {
	static const unsigned int NOID	= std::numeric_limits<unsigned int>::max();		//!< points nowhere
	static const unsigned int TERMINAL = std::numeric_limits<unsigned int>::max()-1;	//!< to terminal
	static const unsigned int ORPHAN = std::numeric_limits<unsigned int>::max()-2;		//!< orphan

	/*! \enum termtype
    	    \brief Symbolic constants for the connectivity of nodes in a flowGraph.
	*/
	enum class termtype {
		none = 0,			///< node is not connected
		source = 1,			///< node is connected to flow source
		sink = 2 };			///< node is connected to flow sink
	//! Represents an node in the flow graph
	struct Node {
		ulist	nb;			//!< outgoing edges
		unsigned int edge;		//!< node's parent edge
		unsigned int it;		//!< iteration showing when distance was computed
		unsigned int dist;		//!< distance to the terminal
		bool	sink;			//!< true if the node is in the sink tree
		T	cap;			//!< cap > 0 denotes residual capacity (source to node) else (node to sink)
		Node() : nb(), edge(NOID), it(0), dist(0), sink(0), cap(0) { }
	};
	//! Represents an edge in the flow graph
	struct Edge {
		unsigned int node;		//!< node the edge points to
		unsigned int rev;		//!< reverse edge
		T	cap;			//!< residual capacity
		Edge(const unsigned int h, const unsigned int s, const T c)
			: node(h), rev(s), cap(c) { }
	};
	const T eps = 1e-6;			//!< small residual capacity
	std::vector<Node> nd;			//!< nodes of this graph
	std::vector<Edge> ed;			//!< edges of this graph
	unsigned int it;			//!< number of iterations
	T	flow;				//!< total flow

	bool	nodeSat(const unsigned int n) const					//! returns true if node capacity is saturated.
		{ return std::abs(nd[n].cap) < eps; }
	bool	edgeSat(const unsigned int e) const					//! returns true if edge capacity is saturated.
		{ return std::abs(ed[e].cap) < eps; }
	unsigned int getDist(const unsigned int n0)					//! returns distance to terminal node (or max).
		{ unsigned int d = 0, e; 
		  for (unsigned int n = n0; ; n = ed[e].node) {
			Node& t = nd[n]; if (t.it == it) { d += t.dist; break; }
			e = t.edge; d++; if (e == ORPHAN) return std::numeric_limits<unsigned int>::max();
			if (e == TERMINAL) { t.it = it; t.dist = 1; break; } };
		  return d; }
	void	setMarks(const unsigned int e, unsigned int d);
	T	bottleneckCap(const unsigned int e);
	void	augmentSource(ulist& O, const unsigned int e, const T bn);
	void	augmentSink(ulist& O, const unsigned int e, const T bn);
	unsigned int growSource(ulist& A, const Node& t);
	unsigned int growSink(ulist& A, const Node& t);
	void	processSourceOrphan(ulist& A, ulist& O, const unsigned int n);
	void	processSinkOrphan(ulist& A, ulist& O, const unsigned int n);
	void	init(ulist& A);
public:
	flowGraph<T>() :
		nd(), ed(), it(0), flow(0)
		{ }
	~flowGraph() = default;
	int	addNode(const unsigned int n = 1)					//! adds n nodes to the graph.
		{ for (unsigned int i = 0; i < n; i++) nd.push_back({});
		  return nd.size()-1; }
	void	addEdge(const unsigned int i, const unsigned int j, const T s,
			const T t)							//! adds edge (i,j) to the graph with capacities (s,t).
		{ assert(i < nd.size() && j < nd.size() && i != j);
		  assert(s >= 0 && t >= 0); unsigned int n = ed.size();
		  ed.push_back({j,n+1,s}); ed.push_back({i,n,t});
		  nd[i].nb.push_front(n); nd[j].nb.push_front(n+1); }
	void	addTweights(const unsigned int i, T s, T t)				//! links node i to terminal with capacities (s,t).
		{ assert(i < nd.size()); const T d = nd[i].cap;
		  if (d > 0) s += d; else t -= d;
		  flow += s < t? s: t; nd[i].cap = s-t; }
	bool	segmentAt(const unsigned int i)	const					//! returns true if node i belongs to the sink tree.
		{ return nd[i].sink; }
	T	work();
};

template <typename T>
T flowGraph<T>::bottleneckCap(const unsigned int m)
//! returns the bottleneck capacity for a path starting at edge m.
{	T bn = ed[m].cap; unsigned int n = ed[ed[m].rev].node;
	while (true) { unsigned int e = nd[n].edge; if (e == TERMINAL) break;		// now in the source tree
		bn = std::min(bn, ed[ed[e].rev].cap); n = ed[e].node; };
	bn = std::min(bn, nd[n].cap); n = ed[m].node;
	while (true) { unsigned int e = nd[n].edge; if (e == TERMINAL) break;		// now in the sink tree
		bn = std::min(bn, ed[e].cap); n = ed[e].node; };
	bn = std::min(bn,-nd[n].cap); return bn;
}

template <typename T>
void flowGraph<T>::augmentSource(ulist& O, const unsigned int m, const T bn)
//! fills the source path starting at edge m with bottleneck capacity bn.
{	ed[ed[m].rev].cap += bn; unsigned int e;					// augment flow in the source tree
	for (unsigned int n = ed[ed[m].rev].node; ; n = ed[e].node) {
		Node& t = nd[n]; e = t.edge;
		if (e == TERMINAL) { t.cap -= bn;
			if (nodeSat(n)) { t.edge = ORPHAN; O.push_front(n); };		// add node to adoption list
			return; };
		ed[e].cap += bn; ed[ed[e].rev].cap -= bn;
		if (edgeSat(ed[e].rev)) { t.edge = ORPHAN; O.push_front(n); } };	// add node to adoption list
}

template <typename T>
void flowGraph<T>::augmentSink(ulist& O, const unsigned int m, const T bn)
//! fills the sink path starting at edge m with bottleneck capacity bn.
{	ed[m].cap -= bn; unsigned int e;
	for (unsigned int n = ed[m].node; ; n = ed[e].node) {
		Node& t = nd[n]; e = t.edge;
		if (e == TERMINAL) { t.cap += bn;
			if (nodeSat(n)) { t.edge = ORPHAN; O.push_front(n); };		// add node to adoption list
			return; };
		ed[e].cap -= bn; ed[ed[e].rev].cap += bn;
		if (edgeSat(e)) { t.edge = ORPHAN; O.push_front(n); } };		// add node to adoption list
}

template <typename T>
void flowGraph<T>::setMarks(const unsigned int m, unsigned int d)
//! sets marks along the path starting at edge m.
{	unsigned int e = NOID;
	for (unsigned int n = ed[m].node; nd[n].it != it; n = ed[e].node) {		// set marks along the path
		Node& t = nd[n]; t.it = it; t.dist = d--; e = t.edge; };
}

template <typename T>
void flowGraph<T>::processSourceOrphan(ulist& A, ulist& O, const unsigned int o)
//! finds a new parent for orphan node o in the source tree.
{	Node& t = nd[o]; unsigned int emin = NOID, dmin = std::numeric_limits<unsigned int>::max();
	for (const auto e : t.nb) { if (edgeSat(ed[e].rev)) continue;
		const unsigned int n = ed[e].node;
		const Node& u = nd[n]; if (u.sink || u.edge == NOID) continue;		// check the origin of n
		const unsigned int d = getDist(n); if (d == std::numeric_limits<unsigned int>::max()) continue;
		if (d < dmin) { emin = e; dmin = d; }; setMarks(e,d); };		// n originates from the source - done
	t.edge = emin; if (emin != NOID) { t.it = it; t.dist = dmin+1; return; }
	for (const auto e : t.nb) { const unsigned int n = ed[e].node; 			// try to find a new edge
		Node& u = nd[n];  if (u.sink || u.edge == NOID) continue;
		if (edgeSat(ed[e].rev) == false) A.push_back(n);
		const unsigned int s = u.edge;
		if (s != TERMINAL && s != ORPHAN && ed[s].node == o) {
			u.edge = ORPHAN; O.push_back(n); } };				// add node to the adoption list
}

template <typename T>
void flowGraph<T>::processSinkOrphan(ulist& A, ulist& O, const unsigned int o)
//! finds a new parent for orphan node o in the sink tree.
{	Node& t = nd[o]; unsigned int emin = NOID, dmin = std::numeric_limits<unsigned int>::max();
	for (const auto e : t.nb) { if (edgeSat(e)) continue;				// try to find a new edge
		const unsigned int n = ed[e].node;
		const Node& u = nd[n]; if (!u.sink || u.edge == NOID) continue; 	// check the origin of n
		const unsigned int d = getDist(n); if (d == std::numeric_limits<unsigned int>::max()) continue;
		if (d < dmin) { emin = e; dmin = d; }; setMarks(e,d); };
	t.edge = emin; if (emin != NOID) { t.it = it; t.dist = dmin+1; return; }
	for (const auto e : t.nb) { const unsigned int n = ed[e].node;
		Node& u = nd[n]; if (!u.sink || u.edge == NOID) continue;
		if (edgeSat(e) == false) A.push_back(n);
		const unsigned int s = u.edge;
		if (s != TERMINAL && s != ORPHAN && ed[s].node == o) {
			u.edge = ORPHAN; O.push_back(n); } };				// add node to the adoption list
}

template <typename T>
unsigned int flowGraph<T>::growSource(ulist& A, const Node& t)
//! grows source tree starting at node& t; returns path ending at edge e (or NOID).
{	for (const auto e : t.nb) {							// try to find a new edge
		const unsigned int r = ed[e].rev; if (edgeSat(e)) continue;
		const unsigned int n = ed[e].node; Node& u = nd[n];
		if (u.edge == NOID) { u.sink = 0; u.edge = r; 
			u.it = t.it; u.dist = t.dist+1; A.push_back(n); }
		else if (u.sink) return e;
		else if (u.it <= t.it && u.dist > t.dist) {				// try to shorten the distance from j to the source
			u.edge = r; u.it = t.it; u.dist = t.dist+1; } };
	return NOID;
}

template <typename T>
unsigned int flowGraph<T>::growSink(ulist& A, const Node& t)
//! grows sink tree starting at node& t; returns path ending at edge e (or NOID).
{	for (const auto e : t.nb) {							// try to find a new edge
		const unsigned int r = ed[e].rev; if (edgeSat(r)) continue;
		const unsigned int n = ed[e].node; Node& u = nd[n];
		if (u.edge == NOID) { u.sink = 1; u.edge = r;
			u.it = t.it; u.dist = t.dist+1; A.push_back(n); }
		else if (!u.sink) return r;
		else if (u.it <= t.it && u.dist > t.dist) {				// try to shorten the distance from j to the source
			u.edge = r; u.it = t.it; u.dist = t.dist+1; } };
	return NOID;
}

template <typename T>
void flowGraph<T>::init(ulist& A)
//! initializes terminal links and list of active nodes A.
{	for (unsigned int n = 0; n < nd.size(); n++) { Node& t = nd[n];
		if (t.cap > 0) { t.sink = 0; A.push_back(n);
			t.edge = TERMINAL; t.dist = 1;}					// n is connected to the source
		else if (t.cap < 0) { t.sink = 1; A.push_back(n);
			t.edge = TERMINAL; t.dist = 1; }				// n is connected to the sink
		else t.edge = NOID; };
}

template <typename T>
T flowGraph<T>::work()
//! flushes the graph until the maximum capacity is reached.
{	ulist A, O;  
	for (init(A); A.empty() == false; A.pop_front()) {
		const unsigned int n = A.front(); if (nd[n].edge == NOID) continue;
		const unsigned int e = nd[n].sink? growSink(A,nd[n]): growSource(A,nd[n]);
		if (e == NOID) continue;
		T f = bottleneckCap(e); augmentSource(O,e,f); augmentSink(O,e,f);
		it++; flow += f;
		while (O.empty() == false) { const unsigned int n = O.front(); O.pop_front();
			if (nd[n].sink) processSinkOrphan(A,O,n);
			else processSourceOrphan(A,O,n); } };
	return flow;
}

//! Provides an interface for energy optimization using a max flow / min cut graph.

template <class T> class mincutEnergy: public flowGraph<T> {
	T	energy;				//!< constant energy of this problem.
public:
	mincutEnergy<T>()
		 : energy(0) { }
	unsigned int addVariable(const unsigned int n = 1)				//! adds a variable to the optimization problem.
		{ return this->addNode(n); }
	void	addConstant(const T e)							//! adds a constant energy to the problem.
		{ energy += e; }
	void	addTerm1(const unsigned int i, const T e0, const T e1)			//! adds an unary energy term to the problem.
		{ this->addTweights(i,e1,e0); }
	void	addTerm2(const unsigned int i, const unsigned int j, 
	                const T e00, const T e01, const T e10, const T e11)		//! adds an binary energy term to the problem.
		{ this->addTweights(i,e11,e00); T d1 = e01-e00, d2 = e10-e11; assert(d1+d2 >= 0);
		  if (d1 < 0) { this->addTweights(i,0,d1); this->addTweights(j,0,-d1);
			this->addEdge(i,j,0,d1+d2); }
		  else if (d2 < 0) { this->addTweights(i,0,-d2); this->addTweights(j,0,d2);
			this->addEdge(i,j,d1+d2,0); }
		  else this->addEdge(i,j,d1,d2); };
	void	addTerm3(const unsigned int i, const unsigned int j, const unsigned int k,
	               const T e000, const T e001, const T e010, const T e011,
	               const T e100, const T e101, const T e110, const T e111)		//! adds a ternary energy term to the problem.
		{ T d0 = (e000+e011+e101+e110)-(e100+e010+e001+e111);
		  if (d0 >= 0) { energy += e111-(e011+e101+e110);
			this->addTweights(i,e101,e001); this->addTweights(j,e110,e100); this->addTweights(k,e011,e010);
			const T d1 = (e010+e001)-(e000+e011); assert(d1 >= 0); this->addEdge(j,k,d1,0);
			const T d2 = (e100+e001)-(e000+e101); assert(d2 >= 0); this->addEdge(k,i,d2,0);
			const T d3 = (e100+e010)-(e000+e110); assert(d3 >= 0); this->addEdge(i,j,d3,0);
		  	if (d0 > 0) { unsigned int u = addVariable();
				this->addEdge(i,u,d0,0); this->addEdge(j,u,d0,0);
				this->addEdge(k,u,d0,0); this->addTweights(u,0,d0); } }
		  else { energy += e000-(e100+e010+e001);
			this->addTweights(i,e110,e010); this->addTweights(j,e011,e001); this->addTweights(k,e101,e100);
			const T d1 = (e110+e101)-(e100+e111); assert(d1 >= 0); this->addEdge(k,j,d1,0);
			const T d2 = (e110+e011)-(e010+e111); assert(d2 >= 0); this->addEdge(i,k,d2,0);
			const T d3 = (e101+e011)-(e001+e111); assert(d3 >= 0); this->addEdge(j,i,d3,0);
			unsigned int u = addVariable();
			this->addEdge(u,i,-d0,0); this->addEdge(u,j,-d0,0);
			this->addEdge(u,k,-d0,0); this->addTweights(u,-d0,0); } }
	T	minimize()								//! optimize the problem using maximum flow.
		{ return energy+this->work(); }
	unsigned int getVar(const unsigned int i) const					//! returns segment of variable i.
		{ return this->segmentAt(i); }
};

#endif
