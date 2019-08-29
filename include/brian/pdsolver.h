#ifndef PDSOLVER_H
#define PDSOLVER_H

/*
 *
 * pdsolver.h: Primal-dual method for optimizing Markov-random fields
 * BRIAN Software Package Version 3.0
 *
 * $Id: pdsolver.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (07/06/16): initial version FK
 * v406 (28/09/16): bumped to version 3.0
 *
 * refer to:
 *
 * Komodakis N, Tziritas G (2007)
 * Approximate Labeling via Graph-Cuts Based on Linear Programming.
 * IEEE PAMI.
 *
 * Komodakis N, Tziritas G, Paragios N (2008)
 * Performance vs Computational Efficiency for Optimizing Single and Dynamic MRFs: 
 * Setting the State of the Art with Primal Dual Strategies.
 * Computer Vision and Image Understanding.
 *
 */

/*! \file
    \brief Optimize Markov-random fields using a primal-dual method.
*/

#define allNodes(i)	(unsigned int i = 0; i < nv; i++)
#define allPairs(i)	(unsigned int i = 0; i < np; i++)
#define allLabels(i)	(unsigned int i = 0; i < nl; i++)

#define NOID	 std::numeric_limits<unsigned int>::max()	// points nowhere
#define TERMINAL std::numeric_limits<unsigned int>::max()-1	// to terminal
#define ORPHAN   std::numeric_limits<unsigned int>::max()-2	// orphan

//! Implements a context for optimizing Markov-random fields using graph cuts.

template <class T> class flowGraph {
public:
	using ulist = std::list<unsigned int>;

	//! Represents an edge in the graph cut algorithm.

	struct Edge {
		unsigned int node;		//!< node this edge points to
		unsigned int rev;		//!< reverse edge
		T	rcap;			//!< residual capacity
		T	cap;			//!< capacity

		Edge(const unsigned int n = 0, const unsigned int r = 0)		//! Allocates an edge for node n and reverse edge r.
			: node(n), rev(r), rcap(0), cap(0) { }  
	};

	//! Represents a node in the graph cut algorithm.

	struct Node {
		ulist	nb;			//!< outgoing edges
		unsigned int edge;		//!< node's parent
		unsigned int dist;		//!< distance to the terminal
		bool	sink;			//!< true iff node is connected to sink
		T	cap;			//!< residual capacity
		int	ts;			//!< timestamp showing when dist was computed
		int	conflict;		//!< timestamp when conflict occured

		Node()									//! Allocates an empty node.
			 : nb(), edge(NOID), dist(0), sink(false), cap(0), ts(0), conflict(-1)
			{ }
		void	init()								//! initializes this node as terminal.
			{ ts = 0; dist = 1; edge = TERMINAL;
			  if (cap == 0) edge = NOID; else sink = cap < 0; }
	};
	Edge*	ed;				//!< array of edges
	Node*	nd;				//!< array of nodes
	unsigned int nv;			//!< number of nodes
	ulist	A;				//!< list of active nodes
	ulist	O;				//!< list of orphan nodes
	int	time;				//!< monotonically increasing global counter
	T	flow;				//!< total flow

	void	init()									//! initializes flow problem.
		{ A.clear(); O.clear(); time = 0;
		  for (unsigned int i = 0; i < nv; i++) {
			nd[i].init(); if (nd[i].edge != NOID) A.push_back(i); } }
	unsigned int growSource(const unsigned int n);
	unsigned int growSink(const unsigned int n);
	void	augment(const unsigned int i, const bool add);
	void	process(const unsigned int n);
	void	processOrphans()							//! decides about orphan nodes.
		{ while (O.empty() == false) {
			const unsigned int n = O.front(); O.pop_front();
			process(n); } }
        flowGraph<T>(const flowGraph<T>& b) = delete;
        flowGraph<T> operator=(const flowGraph<T>& b) = delete;
	flowGraph(Node* _nd, Edge* _ed, unsigned int _nv)				//! Allocates a context for solving the flow problem on nodes nd and edges ed.
		: ed(_ed), nd(_nd), nv(_nv), A(), O(), time(0), flow(0)
		{ }
	void	addEdges(unsigned int* pairs, unsigned int np)				//! allocates np edges referenced in vector pairs.
		{ for (unsigned int i = 0; i < 2*np; i += 2) {
			const unsigned int f = pairs[i], r = pairs[i+1];
			nd[f].nb.push_back(i); nd[r].nb.push_back(i+1);
			ed[i] = Edge(r,i+1); ed[i+1] = Edge(f,i); } }
	void	addWeights(Node* n, T source, T sink)					//! adds weights to node n.
		{ T d = n->cap; if (d > 0) source += d; else sink -= d;
		  flow += source < sink? source: sink; n->cap = source-sink; }
	T	work(const bool in, const bool run);
};

template <typename T>
void flowGraph<T>::process(const unsigned int k)
//! process state of node k.
{	unsigned int emin = NOID, dmin = std::numeric_limits<unsigned int>::max();
	Node& n = nd[k]; for (const auto i : n.nb) { Edge& e = ed[i];
		T cap = n.sink? e.rcap: ed[e.rev].rcap; if (cap == 0) continue;
		unsigned int p = e.node; unsigned int a = nd[p].edge;
		if (n.sink != nd[p].sink || a == NOID) continue;
		unsigned int d = 0; while (1) {
			if (nd[p].ts == time) { d += nd[p].dist; break; }
			a = nd[p].edge; d++;
			if (a == TERMINAL) { nd[p].ts = time; nd[p].dist = 1; break; }
			if (a == ORPHAN) { d = std::numeric_limits<unsigned int>::max(); break; }
			p = ed[a].node; }
		if (d == std::numeric_limits<unsigned int>::max()) continue;
		if (d < dmin) { emin = i; dmin = d; };
		for (p = e.node; nd[p].ts != time; p = ed[nd[p].edge].node) {
			nd[p].ts = time; nd[p].dist = d--; } };
	n.edge = emin; if (emin != NOID) { n.ts = time; n.dist = dmin+1; }
	else {	n.ts = 0; for (const auto i : n.nb) { Edge& e = ed[i];
		Node& p = nd[e.node]; unsigned int a = p.edge;
		if (n.sink != p.sink || a == NOID) continue;
		T cap = n.sink? e.rcap: ed[e.rev].rcap;	if (cap) A.push_back(e.node);
		if (a != TERMINAL && a != ORPHAN && ed[a].node == k) {
			p.edge = ORPHAN; O.push_front(e.node); } } }
}

template <typename T>
void flowGraph<T>::augment(const unsigned int k, const bool add)
//! augment flow in node k.
{	unsigned int n,e; Edge& me = ed[k]; T bn = me.rcap;
	for (n = ed[me.rev].node; ; n = ed[e].node) {
		e = nd[n].edge; if (e == TERMINAL) break;
		bn = std::min(bn,ed[ed[e].rev].rcap); };
	if (bn > nd[n].cap) bn = nd[n].cap;
	for (n = me.node; ; n = ed[e].node) {
		e = nd[n].edge; if (e == TERMINAL) break;
		bn = std::min(bn,ed[e].rcap); };
	if (bn > -nd[n].cap) bn = -nd[n].cap;
	ed[me.rev].rcap += bn; me.rcap -= bn;
	for (n = ed[me.rev].node; ; n = ed[e].node) {
		e = nd[n].edge; if (e == TERMINAL) break;
		ed[e].rcap += bn; ed[ed[e].rev].rcap -= bn;
		if (ed[ed[e].rev].rcap == 0) { nd[n].edge = ORPHAN; O.push_front(n); } };
	nd[n].cap -= bn; if (nd[n].cap == 0) { nd[n].edge = ORPHAN; O.push_front(n); };
	for (n = me.node; ; n = ed[e].node) {
		e = nd[n].edge; if (e == TERMINAL) break;
		ed[ed[e].rev].rcap += bn; ed[e].rcap -= bn;
		if (ed[e].rcap == 0) { nd[n].edge = ORPHAN; O.push_front(n); } };
	nd[n].cap += bn; flow += bn;
	if (nd[n].cap == 0) { if (add) { nd[n].edge = ORPHAN; O.push_front(n); }
		else nd[n].edge = NOID; };	
}

template <typename T>
unsigned int flowGraph<T>::growSource(const unsigned int k)
//! grow source tree from node k.
{	const Node& n = nd[k]; for (const auto i : n.nb) {
		Edge& e = ed[i]; Node& p = nd[e.node]; if (e.rcap == 0) continue;
		if (p.edge == NOID) { p.sink = false; p.edge = e.rev;
			p.ts = n.ts; p.dist = n.dist+1; A.push_back(e.node); }
		else if (p.sink) return i;
		else if (p.ts <= n.ts && p.dist > n.dist) {
			p.edge = e.rev; p.ts = n.ts; p.dist = n.dist+1; } };
	return NOID;
}

template <typename T>
unsigned int flowGraph<T>::growSink(const unsigned int k)
//! grow sink tree from node k.
{	const Node& n = nd[k]; for (const auto i : n.nb) { 
		Edge& e = ed[i]; Node& p = nd[e.node]; if (ed[e.rev].rcap == 0) continue;
		if (p.edge == NOID) { p.sink = true; p.edge = e.rev;
			p.ts = n.ts; p.dist = n.dist+1; A.push_back(e.node); }
		else if (p.sink == false) return e.rev;
		else if (p.ts <= n.ts && p.dist > n.dist) {
			p.edge = e.rev; p.ts = n.ts; p.dist = n.dist+1; } };
	return NOID;
}

template <typename T>
T flowGraph<T>::work(const bool in, const bool run)
//! optimize flow problem using graph cuts.
{	if (in) init();	unsigned int k = 0;
	while (true) { if (k == 0 || nd[k].edge == NOID) { 
			while (true) { if (A.empty()) return flow;
				k = A.front(); A.pop_front();
				if (nd[k].edge != NOID) break; } }
		unsigned int e = NOID;
		if (nd[k].sink == false) e = growSource(k);
		else if (run) e = growSink(k);
		time++;
		if (e != NOID) { augment(e,run); processOrphans(); }
		else k = 0; }
	return flow;
}

template<typename T> using fg = flowGraph<T>;
template<typename T> using fgNode = typename flowGraph<T>::Node;
template<typename T> using fgEdge = typename flowGraph<T>::Edge;

//! Implements a context for optimizing Markov-random fields using a primal-dual method.

template<typename T> class PDSolver {
private:

	//! Helper structure for node information.

	struct NodeInfo {
		unsigned int label;		//!< current label.
		T	height;			//!< active height of node.
		int	time;			//!< timestamp of change.
 		int	next;    		//!< next node.
		int	prev;    		//!< previous node.
		int*	pairs;			//!< neighboring edges.
		unsigned int np;		//!< number of edges.
	};

	//! Helper structure for pair information.

	struct PairInfo {
		unsigned int i0;		//!< first node.
		unsigned int i1;		//!< second node.
		int	time;			//!< timestamp of change.
		PairInfo(const unsigned int _i0 = 0, const unsigned int _i1 = 0,
			const int _time = -1)						//! Allocates a pair for nodes (_i0,_i1) at time point time.
			: i0(_i0), i1(_i1), time(_time) { }
	};

	//! Helper structure for edge information.

	struct EdgeInfo {
		int	head;			//!< head node.
		int	tail;			//!< tail node.
		T	balance;		//!< active height.
		EdgeInfo(int _head = 0, int _tail = 0, T bal = 0)			//! Allocates an edge for nodes (_head,_tail) with height bal.
			: head(_head), tail(_tail), balance(bal) { }
	};

	const unsigned int nv;			//!< number of vertices.
	const unsigned int nl;			//!< number of levels (states).
	const unsigned int np;			//!< number of edges.
        matD<T> h;				//!< height variables.
	const matD<T>& dist;			//!< distance function for pairwise potential.
        const uvecD& pairs;			//!< list of vertex pairs.
	const vecD<T>& wcosts;			//!< weights of MRF pairwise potentials.

	matD<T> y;				//!< balance variables.
	std::vector<fgNode<T>> nodes;		//!< nodes and edges.
	std::vector<fgEdge<T>> edges;		//!< of max-flow graphs.
	std::vector<fg<T>*> graphs;		//!< vector of flow graphs.
	std::vector<NodeInfo> pinfo;		//!< vector of node information records.
	std::vector<PairInfo> pairInfo;		//!< vector of pair information records.
	std::vector<EdgeInfo> einfo;		//!< vector of edge information records.
	std::vector<unsigned int> children;	//!< children ids.
	ivecD	sourceNodes1;			//!< auxiliary lists for keeping.
	ivecD	sourceNodes2;			//!< track of source-linked nodes.
	ivecD	pairsArr;			//!< array of active pairs

	double	APF;				//!< MRF energy.
	int	time;				//!< current time point.
	int	activeList;			//!< head of active ids.
	int	APF_changeTime;			//!< time point of last change.

	bool	newLabel(const fgNode<T>& n) const					//! allocates a new label for node n.
		{ return n.edge != NOID && n.sink == false; }
	void	init();
	void	firstIteration(const unsigned int label);
	void	secondIteration(const unsigned int label);
	void	trackNodes(const unsigned int label);
	double	getDelta(const unsigned int label, NodeInfo& pi, fgNode<T>* nd);

        PDSolver<T>(const PDSolver<T>& other) = delete;
        PDSolver<T> operator=(const PDSolver<T>& other) = delete;
public:

	PDSolver(const unsigned int _nv, const unsigned int _nl, const unsigned int _np,
		const matD<T>& lcosts, const matD<T>& dist, const uvecD& pairs, 
		const vecD<T>& wcosts);
	~PDSolver()
		{ for (unsigned int i = 0; i < nl; i++) delete graphs[i]; }
	uvecD	work(const unsigned int maxit);
}; 

template <typename T>
PDSolver<T>::PDSolver(const unsigned int _nv, const unsigned int _nl, const unsigned int _np,
	const matD<T>& lcosts, const matD<T>& _dist, const uvecD& _pairs, const vecD<T>& _wcosts)
	: nv(_nv), nl(_nl), np(_np), h(lcosts),	dist(_dist), pairs(_pairs), wcosts(_wcosts),
	y(matD<T>(np,nl)), nodes(nv*nl), edges(np*nl*2), graphs(nl), pinfo(nv), pairInfo(np),
	einfo(np), children(nv), sourceNodes1(nv), sourceNodes2(nv), pairsArr(np*2), APF(0.0),
	time(-1), activeList(-1), APF_changeTime(-1)
//! Allocates a context for optimizing a Markov random field described by nv vertices, nl levels, np edges, labeling costs lcosts, edge costs _dists, edges _pairs and edge weights _wcosts.
{	sourceNodes1 = -2; sourceNodes2 = -2;
	for allLabels(i) { if (dist(i,i)) throw optException("Potential (a,a) must be 0"); 
		fg<T>* g = new fg<T>(&nodes[nv*i], &edges[np*i*2], nv);
		g->addEdges(pairs.x, np); graphs[i] = g; };
	for allNodes(i) pinfo[i].np = 0;
	for allPairs(i) { pinfo[pairs(i*2)].np++; pinfo[pairs(i*2+1)].np++; }
	unsigned int offset = 0;
	for allNodes(i) { pinfo[i].pairs = &pairsArr[offset];
		offset += pinfo[i].np; pinfo[i].np = 0; }
	for allPairs(i) {
		const unsigned int i0 = pairs(i*2), i1 = pairs(i*2+1);
		pinfo[i0].pairs[pinfo[i0].np++] =  i;
		pinfo[i1].pairs[pinfo[i1].np++] = -i;
		einfo[i] = EdgeInfo(i1,i0); pairInfo[i] = PairInfo(i0, i1); }
}

template <typename T>
void PDSolver<T>::init()
//! initializes all graph structures.
{	for allNodes(i) {
		pinfo[i].label = 0; pinfo[i].time = -1;
		pinfo[i].next = -1; pinfo[i].prev = -2; }
	y = 0; for allPairs(i) {
		const int id0 = einfo[i].tail, id1 = einfo[i].head;
		const int l0 = pinfo[id0].label, l1 = pinfo[id1].label;
		const T d = wcosts(i)*dist(l0,l1) - (y(i,l0)-y(i,l1));
		if (l0 == l1) assert(d == 0);
		else { y(i,l0) += d; h(id0,l0) += d; h(id1,l0) -= d; }
		einfo[i].balance = -y(i,l1); }
	APF = 0.0; for allNodes(i) {
		pinfo[i].height = h(i,pinfo[i].label); APF += pinfo[i].height; }
}

template <typename T>
double PDSolver<T>::getDelta(const unsigned int label, NodeInfo& pi, fgNode<T>* nd)
//! returns primal-dual difference.
{	T td = 0.0; for (unsigned int k = 0; k < pi.np; k++) {
		if (pi.pairs[k] > 0) continue;
		const unsigned int pid = std::abs(pi.pairs[k]); 
		const int i0 = pairInfo[pid].i0; if (newLabel(nd[i0])) continue;
		const unsigned int l0 = pinfo[i0].label, l1 = pi.label;
		const T d = wcosts(pid)*(dist(l0,label)+dist(label,l1)-dist(l0,l1));
		if (d < 0) { y(pid,label) -= d; td += d; nd[i0].cap += d; } };
	return td;
}

template <typename T>
void PDSolver<T>::firstIteration(const unsigned int label)
//! runs the first iteration of the optimization problem for this label.
{	fgNode<T>* nd = &nodes[nv*label]; fgEdge<T>* ed = &edges[np*label*2];
	time++; if (APF_changeTime < time-int(nl)) return;
	for allPairs(i) { EdgeInfo& ei = einfo[i];
		const unsigned int l0 = pinfo[ei.tail].label, l1 = pinfo[ei.head].label; 
		if (l1 != label && l0 != label)  {
			T d = wcosts(i)*(dist(label,l1)+dist(l0,label)-dist(l0,l1));
			const T d1 = wcosts(i)*dist(label,l1)-(y(i,label)+ei.balance), d2 = d-d1;
			if (d1 < 0 || d2 < 0)  {
				y(i,label) += d1; h(ei.tail,label) += d1; h(ei.head,label) -= d1; 
				ed[i*2].cap = ed[i*2].rcap = 0;
				if (d < 0) { d = 0; nd[ei.head].conflict = time; }
				ed[i*2+1].rcap = d; }
			else { ed[i*2].cap = ed[i*2].rcap = d1; ed[i*2+1].rcap = d2; } }
		else { ed[i*2].cap = ed[i*2].rcap = 0; ed[i*2+1].rcap = 0; } }
	T tc = 0.0;
	for allNodes(i) { const T d = pinfo[i].height-h(i,label);
		nd[i].cap = d; if (d > 0) tc += d; }
	fg<T>* g = graphs[label]; g->flow = 0; T f = g->work(true,true);
	APF -= tc-f; if (tc > f) APF_changeTime = time;
	for allPairs(i) { EdgeInfo& ei = einfo[i]; 
		if (pinfo[ei.head].label != label && pinfo[ei.tail].label != label) {
			if (newLabel(nd[ei.head])) ei.balance = -(y(i,label)+ ed[i*2].cap - ed[i*2].rcap); }
		else if (pinfo[ei.head].label != label) {
			if (newLabel(nd[ei.head])) ei.balance = -y(i,label); } }
	for allNodes(i) { NodeInfo& pi = pinfo[i]; 
		if (pi.label != label && newLabel(nd[i])) {
			if (nd[i].conflict > pi.time) {
				for (unsigned int k = 0; k < pi.np; k++) {
					if (pi.pairs[k] > 0) continue;
					const unsigned int pid = std::abs(pi.pairs[k]); 
					const int i0 = pairInfo[pid].i0;  if (newLabel(nd[i0])) continue;
					const unsigned int l0 = pinfo[i0].label, l1 = pi.label;
					const T d = wcosts(pid)*(dist(l0,label)+dist(label,l1)-dist(l0,l1));
					if (d < 0) { y(pid,label) -= d; einfo[pid].balance = -y(pid,label);
						pi.height += d; APF += d; nd[i0].cap += d; } } };
			pi.label = label; pi.height -= nd[i].cap;
			nd[i].cap = 0; pi.time = time; }
		h(i,label) = pi.height; }
}

template <typename T>
void PDSolver<T>::secondIteration(const unsigned int label)
//! runs the subsequent iterations of the optimization problem for this label.
{	fgNode<T>* nd = &nodes[nv*label]; fgEdge<T>* ed = &edges[np*label*2];
	time++; const int dt = time-int(nl); if (APF_changeTime < dt) return;
	for allPairs(i) {
		const unsigned int i0 = pairs(i*2), i1 = pairs(i*2+1);
		if (pinfo[i0].time >= dt || pinfo[i1].time >= dt) {
			if (h(i0,label) != pinfo[i0].height) {
				const T hi = h(i0,label)-nd[i0].cap;
				nd[i0].cap = pinfo[i0].height-hi;
				h(i0,label) = pinfo[i0].height; };
			if (h(i1,label) != pinfo[i1].height) {
				const T hi = h(i1,label)-nd[i1].cap;
				nd[i1].cap = pinfo[i1].height-hi;
				h(i1,label) = pinfo[i1].height; }
			const unsigned int l0 = pinfo[i0].label, l1 = pinfo[i1].label;
			if (l0 != label && l1 != label) {
				fgEdge<T>& ed1 = edges[(np*l1+i)*2];
				const T y_pq = y(i,label)+ed[i*2].cap-ed[i*2].rcap;
				const T y_qp = -(y(i,l1)+ed1.cap-ed1.rcap);
				T delta  = wcosts(i)*(dist(label,l1)+dist(l0,label)-dist(l0,l1));
				const T delta1 = wcosts(i)*dist(label,l1)-(y_pq+y_qp);
				const T delta2 = delta-delta1;
				if (delta1 < 0 || delta2 < 0)  {
					y(i,label) = y_pq+delta1; ed[i*2].rcap = ed[i*2].cap = 0;
					if (delta < 0) { delta = 0; nd[i1].conflict = time; }
					ed[i*2+1].rcap = delta;
					nd[i0].cap -= delta1; nd[i1].cap += delta1; }
				else {	y(i,label) = y_pq; ed[i*2].rcap = ed[i*2].cap = delta1;
					ed[i*2+1].rcap = delta2; } }
			else { y(i,label) += ed[i*2].cap - ed[i*2].rcap; } } }
	fg<T>* g = graphs[label]; g->flow = 0; g->work(true,false);
	double prevAPF = APF;
	for allNodes(i) { NodeInfo& pi = pinfo[i];
		if (newLabel(nd[i]) == false) continue;
		if (nd[i].conflict > pi.time)  nd[i].cap -= getDelta(label,pi,nd);
		pi.height -= nd[i].cap; 
		APF -= nd[i].cap; pi.time = time; pi.label = label;
		if (pi.prev == -2) { pi.next = activeList; pi.prev = -1;
			if (activeList >= 0) pinfo[activeList].prev = i;
			activeList = i; } }
	if (APF < prevAPF) APF_changeTime = time;
}

template <typename T>
void PDSolver<T>::trackNodes(const unsigned int label)
//! tracks changes of nodes with this label.
{	fg<T>* g = graphs[label]; g->flow = 0;
	fgNode<T>* nd = &nodes[nv*label]; fgEdge<T>* ed = &edges[np*label*2];
	time++; const int dt = time-int(nl); if (APF_changeTime < dt) return;
	int ss1 = -1, ss2 = -1;	
	for (int i = activeList; i >= 0; ) { NodeInfo& pi = pinfo[i]; int i_next = pi.next;
		if (pi.time >= dt) {
			if (h(i,label) != pi.height) {
				const T hi = h(i,label)-nd[i].cap;
				nd[i].cap = pi.height-hi; h(i,label) = pi.height; }
			if (nd[i].cap) { nd[i].edge = TERMINAL; nd[i].sink = 1; nd[i].dist = 1; }
			else nd[i].edge = NOID; }
		else {	const int prev = pi.prev;
			if (prev >= 0) { pinfo[prev].next = pi.next;
				if (pi.next >= 0) pinfo[pi.next].prev = prev; }
			else {	activeList = pi.next;
				if (activeList >= 0) pinfo[activeList].prev = -1; }
			pi.prev = -2; }
		i = i_next; }
	for (int i = activeList; i >= 0; ) { NodeInfo& pi = pinfo[i]; int i_next = pi.next;
		for (unsigned int k = 0; k < pi.np; k++) {
			unsigned int l0, l1, i0, i1, ii; int pid = pi.pairs[k];
			if (pid >= 0) {	
				if (pairInfo[pid].time == time) continue;
				i0 = i; i1 = pairInfo[pid].i1; l0 = pi.label; l1 = pinfo[i1].label; ii = i1; }
			else {	pid = -pid;
				if (pairInfo[pid].time == time) continue;
				i1 = i; i0 = pairInfo[pid].i0; l1 = pi.label; l0 = pinfo[i0].label; ii = i0; }
			pairInfo[pid].time = time; 
			if (l0 != label && l1 != label) {
				fgEdge<T>& ed1 = edges[(np*l1+pid)*2];
				const T y_pq = y(pid,label)+ed[pid*2].cap-ed[pid*2].rcap ;
				const T y_qp = -(y(pid,l1)+ed1.cap-ed1.rcap);
				T delta  = wcosts(pid)*(dist(label,l1)+dist(l0,label)-dist(l0,l1));
				const T delta1 = wcosts(pid)*dist(label,l1)-(y_pq+y_qp);
				const T delta2 = delta-delta1;
				if (delta1 < 0 || delta2 < 0) {
					y(pid,label) = y_pq+delta1; ed[pid*2].rcap = ed[pid*2].cap = 0;
					if (delta < 0) { delta = 0; nd[i1].conflict = time; }
					ed[pid*2+1].rcap = delta; nd[i0].cap -= delta1; nd[i1].cap += delta1;
					if (pinfo[ii].prev == -2 && sourceNodes2[ii] == -2) {
						sourceNodes2[ii] = ss2;
						ss2 = ii; } }
				else { y(pid,label) = y_pq; ed[pid*2].rcap = ed[pid*2].cap = delta1;
					ed[pid*2+1].rcap = delta2; } }
			else {	y(pid,label) += ed[pid*2].cap-ed[pid*2].rcap;
				ed[pid*2+1].rcap = 0; ed[pid*2].rcap = 0; ed[pid*2].cap = 0; } }
		if (nd[i].cap > 0) {
			nd[i].sink = 0; nd[i].edge = TERMINAL; nd[i].dist = 1;
			g->A.push_back(i); sourceNodes1[i] = ss1; ss1 = i; }
		else if (nd[i].cap < 0) {
			nd[i].sink = 1; nd[i].edge = TERMINAL; nd[i].dist = 1; }
		else nd[i].edge = NOID;
		i = i_next; }
	for (int i = ss2; i >= 0; ) {
		if (nd[i].cap > 0) {
			nd[i].sink = 0; nd[i].edge = TERMINAL; nd[i].dist = 1;
			g->A.push_back(i); sourceNodes1[i] = ss1; ss1 = i; }
		else if (nd[i].cap < 0) {
			nd[i].sink = 1; nd[i].edge = TERMINAL; nd[i].dist = 1; }
		else nd[i].edge = NOID;
		int tmp = i; i = sourceNodes2[i]; sourceNodes2[tmp] = -2; }
	g->work(false,false); double prevAPF = APF; unsigned int nc = 0;
	for (int i = ss1; i >= 0; ) {
		if (nd[i].edge == TERMINAL) { NodeInfo& pi = pinfo[i];
			if (nd[i].conflict > pi.time) 
				nd[i].cap -= getDelta(label,pi,nd);
			pi.height -= nd[i].cap;  APF -= nd[i].cap;
			pi.label = label; pi.time =time;
			if (pi.prev == -2)  {
				pi.next = activeList; pi.prev = -1;
				if (activeList >= 0) pinfo[activeList].prev = i;
				activeList = i; }
			for (const auto j : nd[i].nb) { const fgEdge<T>& e = ed[j];
				if (nd[e.node].edge == e.rev) children[nc++] = e.node; } }
		int tmp = i; i = sourceNodes1[i]; sourceNodes1[tmp] = -2; }
	for (unsigned int i = 0; i < nc; i++) {  
		const unsigned int c = children[i]; NodeInfo& pi = pinfo[c];
		if (nd[c].conflict > pi.time) nd[c].cap -= getDelta(label,pi,nd);
		pi.height -= nd[c].cap; APF -= nd[c].cap; pi.label = label; pi.time = time;
		if (pi.prev == -2) {
			pi.next = activeList; pi.prev = -1;
			if (activeList >= 0) pinfo[activeList].prev = c;
			activeList = c; }
		for (const auto j : nd[c].nb) { const fgEdge<T>& e = ed[j];
			if (nd[e.node].edge == e.rev) children[nc++] = e.node; } }
	if (APF < prevAPF) APF_changeTime = time;
}

template <typename T>
uvecD PDSolver<T>::work(const unsigned int maxit)
//! optimizes the MRF using the primal-dual algorithm.
{	init();
	for (unsigned int it = 0; it < maxit; it++) {
		const double prevAPF = APF; for allLabels(l) {
			if (it == 0) firstIteration(l);
			else if (it == 1) secondIteration(l);
			else trackNodes(l); };
			if (prevAPF <= APF) break; }
	uvecD lbl(nv); for allNodes(i) lbl[i] = pinfo[i].label;
	return lbl;
}

#endif

