#ifndef GRAPHPARTTIONER_H
#define GRAPHPARTTIONER_H

/*
 *
 * graphPartitioner.h: partition a graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: graphPartitioner.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (17/03/17): initial version FK
 *
 */

/*
 * based on code from gmetis
 *
 */

using uvec = std::vector<unsigned int>;

//! Implements a helper structure for partitioning information.

struct ptinfo {
	unsigned int pt;			//!< partition id
	unsigned int mask;			//!< split mask
	float	wt;				//!< partition weight

	ptinfo()									//! initializes an empty partition.
		: pt(~0u), mask(~0u), wt(0.0f) {}
	ptinfo(const float _wt)								//! initializes partition of weight wt.
		: pt(0), mask(1), wt(_wt) {}
	ptinfo(const unsigned int pn, const unsigned int pm, const float pw)		//! initializes partition with id, mask pm, and weight wt.
		: pt(pn), mask(pm), wt(pw) {}
	unsigned int splitID() const							//! returns id for splitting partition.
		{ return pt | mask; }
	float splitRatio(const unsigned int np) const
		{ unsigned int L = 0, R = 0, LM = mask-1;
		  for (unsigned int i = 0; i < np; i++) {
			if ((i & LM) != pt) continue;
			if (i & mask) R++; else L++; }
		  return wt*R/float(L+R); }
	ptinfo split()									//! splits this partition.
		{ ptinfo np(splitID(), mask << 1, 0.0f); mask <<= 1; return np; }
};

//! Implements a context for splitting a graph into partitions using recursive bisection.

template<typename N> struct bisectSplitter {
	//! Implements a helper structure for matching nodes to partitions.
	struct matchInfo {
		unsigned int id1;		//!< node id
		unsigned int id2;		//!< neighbor id
		float	gain;			//!< weight change for adding id2
	};

	const graph<N>& g;			//!< graph to be split
	const unsigned int np;			//!< number of partitions
	std::vector<ptinfo> parts;		//!< partitioning information
	uvecD	pt;				//!< partition vector

	unsigned int findSeed(const ptinfo& op) const					//! returns randomly chosen seed in partition op.
		{ float w = drand48()*op.wt; unsigned int id = 0;			// split at a random weight (=number of edges)
		  for (unsigned int i = 0; i < g.size(); i++) {				// for all nodes in partition op
			if (pt(i) == op.pt) { id = i; w -= g(i)->weight;		// if weight exceeds w
				if (w < 0.0f) return i; } };				// return cut node i
		  return id; }								// else return last node in partition
	ptinfo	bisect(ptinfo& op, uvec& nl)						//! splits partition op, generating a list of nodes nd.
		{ ptinfo sp = op.split(); std::deque<unsigned int> wl;			// split partition op
		  float tw = op.wt*0.5f; sp.wt = 0.0f;
		  do {	wl.push_back(findSeed(op));					// get a random seed in partition op.pt
			while (sp.wt < tw && !wl.empty()) { unsigned int i = wl.front(); 
				wl.pop_front(); if (pt(i) == sp.pt) continue;		// for all nodes id that are not yet in the new partition
				sp.wt += g(i)->weight; pt[i] = sp.pt; nl.push_back(i);	// else add to new partition and list
				for (const auto& e : g(i)->edges)			// and put all neighbors onto worklist
					if (pt(e.id) == op.pt) wl.push_back(e.id); } }
		  while (sp.wt < tw); op.wt -= sp.wt; return sp; }
	float	gain(const unsigned int id, const unsigned int op,
			const unsigned int np) const					//! returns weight gain for moving node id from op to np
		{ float w = 0.0f; for (const auto& e : g(id)->edges) {			// for all neighbors of node id
		  	if (!(pt(e.id) == op || pt(e.id) == np)) continue;
			w += pt(e.id) == pt(id)? -e.weight: e.weight; };		// increase gain if neighbor is in other partition
		  return w; }
	std::vector<matchInfo> matchBnd(const uvec& nl, const unsigned int op,
			const unsigned int np) const					//! finds neighbors of node id that yield the highest weight gain.
		{ std::vector<matchInfo> m;
		  for (const auto i : nl) { if (!(pt(i) == op || pt(i) == np)) continue; // for all nodes in the list
			float mg = 0.0f; unsigned int ei = ~0u;
		  	const float sg = gain(i,op,np);					// get gain for moving node id from partition op to np
			for (const auto& e : g(i)->edges) {				// for all neighbors of node id
				if (pt(i) == pt(e.id)) continue;			// skip if not on boundary
				if (!(pt(e.id) == op || pt(e.id) == np)) continue;	// between op and np
				const float ng = gain(e.id,op,np);			// get gain for neighbor
				const float tg = sg+ng-2.0f*e.weight;			// compute total gain
				if (tg > mg) { mg = tg; ei = e.id; } };			// keep gain and neighbor id for best match
			if (ei != ~0u) m.push_back({i,ei,mg}); };
		  return m; }								// return list of best matches
	void	refine(const uvec& nl, const unsigned int op, const unsigned int np)
		{ for (unsigned int it = 0; it < 100; it++) {
			const auto m = matchBnd(nl,op,np); if (m.size() == 0) return;	// find best matching boundary nodes
		  	float s = 0.0f, gm = 0.0f; unsigned int lim = ~0u;		// determine the last match that still yields a gain
		  	for (unsigned int i = 0; i < m.size(); i++) { 
				s += m[i].gain; if (s > gm) { gm = s; lim = i; } };
		  	if (gm <= 0.0f || lim == ~0u) return;				// no gain, we're done
		  	for (unsigned int i = 0; i <= lim; i++) { 			// for all matches on the list
				const unsigned int id1 = m[i].id1, id2 = m[i].id2;	// get ids of matching nodes
				const unsigned int pt1 = pt(id1), pt2 = pt(id2);	// get partitions of matching nodes
				const float dw = g(id2)->weight-g(id1)->weight;
				parts[pt1].wt += dw; pt[id1] = pt2;			// update weights and partitioning
				parts[pt2].wt -= dw; pt[id2] = pt1; } } } 
public:
	bisectSplitter(const graph<N>& _g, const unsigned int _np)
		: g(_g), np(_np), parts(np), pt(g.size()) { }
	uvecD	work(const unsigned int tw, const bool verbose)
		{ parts[0] = ptinfo(tw); pt = 0;					// initialize partition vector
		  std::stack<unsigned int> wl; wl.push(0);
		  while (!wl.empty()) {	const auto i = wl.top(); wl.pop(); 		// while the worklist is not empty
			if (parts[i].splitID() >= np) continue;				// get a partition from worklist
			uvec nl; ptinfo sp = bisect(parts[i],nl); parts[sp.pt] = sp;	// split partition, generating a set of boundary nodes
			refine(nl,parts[i].pt,sp.pt);					// refine boundary
			wl.push(sp.pt); wl.push(i); };					// push old and new partition on worklist
		  if (verbose) { for (const auto& p : parts)
				printf("pt %d %e\n", p.pt, p.wt); };
		  return pt; }								// return partition vector
};

//! Implements a context for splitting a graph into partitions using linear aggregation.

template<typename N> struct linearSplitter {
	const graph<N>& g;			//!< graph to be split
	const unsigned int np;			//!< number of partitions
	std::vector<ptinfo> parts;		//!< partitioning information
	uvecD	pt;				//!< partition vector
public:
	linearSplitter(const graph<N>& _g, const unsigned int _np)
		: g(_g), np(_np), parts(np), pt(g.size()) { }
	uvecD	work(const unsigned int tw, const bool verbose)
		{ float w = 0.0f, s = 0.0f, wl = float(tw)/np; unsigned int p = 0;
		  for (unsigned int i = 0; i < g.size(); i++) {
			const float wi = g(i)->weight; w += wi; s += wi; pt[i] = p;
			if (s > (p+1)*wl) { parts[p].pt = p; parts[p].wt = w;
				w = 0.0f; p++; } };
		  parts[np-1].pt = p; parts[np-1].wt = w;
		  if (verbose) { for (const auto& p : parts)
				printf("pt %d %e\n", p.pt, p.wt); };
		  return pt; }								// return partition vector
};

//! Implements a context for building a multi-resolution graph.

template<typename N> struct graphLevel {						//! represents a resolution level of the graph
	const graph<N> fg;			//!< fine-level graph
	graph<N> cg;				//!< coarse-level graph
	uvecD	map;				//!< maps node ids fine-to-coarse
	unsigned int size;			//!< size for next level

	using ivp = std::pair<unsigned int,float>;

	ivp	oneHopMatch(const unsigned int id) const				//! returns best matching neighbor for node id.
		{ ivp rv(id,0.0f); 
	  	  for (const auto& e : fg(id)->edges) {
			if (map(e.id) == ~0u && rv.second < e.weight) rv = {e.id,e.weight}; }
		  return rv; }
	ivp	twoHopMatch(const unsigned int id) const				//! returns best matching over-next neighbor for node id.
		{ ivp rv(id,0.0f);
	  	  for (const auto& e : fg(id)->edges) {
			const ivp tv = oneHopMatch(e.id);
			if (tv.first != id && tv.second > rv.second) rv = tv; }
		  return rv; }
	uvec	ldIndexer() const							//! returns node list sorted by increasing number of edges.
		{ std::vector<std::pair<unsigned int,unsigned int>> ar(fg.size());
		  for (unsigned int i = 0; i < fg.size(); i++) ar[i] = {fg(i)->nEdges(),i};
		  std::sort(ar.begin(),ar.end()); uvec ix(fg.size());
		  for (unsigned int i = 0; i < fg.size(); i++) ix[i] = ar[i].second;
		  return ix; }
	unsigned int matchNodes(uvec& loners, const unsigned int mode)			//! matches nodes in the fine graph and adds to coarse graph.
		{ unsigned int rem = 0; const uvec ix = ldIndexer();
		  for (const auto id : ix) { if (map(id) != ~0u) continue;		// for all unmatched nodes in the fine graph
			if (fg(id)->nEdges() == 0) loners.push_back(id);		// if number of neighbors is 0, save as loner
			const unsigned int nb = mode == 2?
				twoHopMatch(id).first: oneHopMatch(id).first;		// get best matching neighbor
			const unsigned int nn = cg.addNode(new N); map[id] = nn;	// add new node to coarse graph
			cg(nn)->weight = fg(id)->weight;
			if (nb == id) rem++;
			else { cg(nn)->weight += fg(nb)->weight; map[nb] = nn; } };	// match found
		  return rem; }
	void	fixupLoners(const uvec& loners)						//! adds unmatched nodes in loners to coarse graph cg.
		{ const unsigned int nl = loners.size();
		  for (unsigned int i = 0; i < nl; i++) {				// for all loners
			const unsigned int l = loners[i], nn = cg.addNode(new N);	// add to coarse graph
			cg(nn)->weight = fg(l)->weight; map[l] = nn;
			if (i != nl-1) { const unsigned int nx = loners[++i];		// merge with next member
                        	cg(nn)->weight += fg(nx)->weight; map[nx] = nn; } } }	// unconditionally add to coarse graph
	void	createCoarseEdges()							//! adds edges to coarse graph based on mapregation vector map.
		{ std::vector<std::map<unsigned int,float>> ed(cg.size());
		  for (unsigned int a = 0; a < map.N; a++) {				// for all parents that point to this node
			unsigned int id = map(a); auto& ei = ed[id];
			for (const auto& e : fg(a)->edges) {				// for all edges of this parent
				const auto d = map(e.id); if (d == id) continue;	// get edge destination node in coarse graph
				auto it = ei.find(d);					// find in current list of edges
				if (it == ei.end()) ei.insert({d,e.weight});		// not found, so add to list
				else it->second += e.weight; } }			// else update edge weight
		  for (unsigned int i = 0; i < cg.size(); i++) {			// for all parents that point to this node
		  	for (const auto& e : ed[i])
				cg(i)->addEdge({e.first,e.second}); } }
	unsigned int pickPartition(const uvecD& pt, const unsigned int id,
			const std::vector<ptinfo>& parts, const unsigned int maxs)	//! finds the partition node id is most connected to.
		{ fvecD ed(parts.size()); ed = 0.0f;
		  for (const auto& e : cg(id)->edges) { const unsigned int pe = pt(e.id);
			if (parts[pe].wt < maxs || pe == pt(id)) ed[pe] += e.weight; }
		  return ed.imax(); }
public:
	graphLevel<N>(const graph<N>& _fg, const bool twoHop)
			: fg(_fg), cg(), map(fg.size()), size(0)
			{ map = ~0u; uvec loners;
			  unsigned int rem = matchNodes(loners,1);
			  fixupLoners(loners);
			  if (twoHop) rem = matchNodes(loners,2);
			  createCoarseEdges(); size = fg.size()/2+rem/2; }
	uvecD	refine(uvecD& pt, std::vector<ptinfo>& parts, const unsigned int mins, const unsigned int maxs)
		{ uvecD fp(fg.size());
		  for (unsigned int id = 0; id < fg.size(); id++) fp[id] = pt(map(id)); // init partition vector from current solution
		  for (unsigned int id = 0; id < cg.size(); id++) {			// for all nodes of the coarse graph
			const unsigned int cp = pt(id);					// current partition
			const unsigned int np = pickPartition(pt,id,parts,maxs);	// select alternate partition
			if (parts[cp].wt < mins || cp == np) continue;			// keep current partition
			pt[id] = np; const float w = cg(id)->weight;			// else move coarse node to alternate partition
			parts[cp].wt -= w; parts[np].wt += w;
			for (unsigned int i = 0; i < fg.size(); i++) {			// move all nodes on the fine graph to the alternate partition
				if (map(i) == id) fp[i] = np; } };
		  return fp; }
	uvecD	edgeCount(const uvecD& pt, const unsigned int np) const
		{ uvecD ed(np); ed = 0;
		  for (unsigned int i = 0; i < pt.N; i++) { if (pt(i) == ~0u) continue;
	  	 	for (const auto& e : fg(i)->edges)
				if (pt(e.id) != ~0u && pt(e.id) != pt(i)) ed[pt(i)]++; }
		  return ed; }
};		

//! Implements a context for partitioning a graph using edge connectivity.

template<typename N> class graphPartitioner {						//! represents a resolution level of the graph
	const graph<N>&	g;			//!< graph to be partitioned
	const unsigned int np;			//!< number of partitions
	unsigned int mins;			//!< minimum partition size
	unsigned int maxs;			//!< maximum partition size
	std::vector<graphLevel<N>> levels;	//!< graph resolution levels
public:
	graphPartitioner(const graph<N>& _g, const unsigned int _np,
		const float imb = 0.01f)						//! allocates a context for graph partitioning.
		: g(_g), np(_np), mins(0), maxs(0), levels()
		{ const float mw = g.size()/float(np); 
		  mins = FTOU(mw*(1.0f-imb)); maxs = FTOU(mw*(1.0f+imb)); }
	uvecD	work(const bool verbose = false)					//! performs partitioning; returns partition vector.
		{ srand48(time(nullptr)); graph<N> fg = g; bool twoHop = false;
		  for (unsigned int id = 0; id < fg.size(); id++) { fg(id)->weight = 1.0f;
			for (auto& e : fg(id)->edges) e.weight = 1.0f; };
		  for (unsigned int size = fg.size(); size > 20*np; ) {			// create resolution levels of g
			levels.push_back({fg,twoHop}); 
			const unsigned int ns = levels.back().size;
			if (verbose) { printf("[%2zd] nv %6d %s hop\r", levels.size(),
					 ns, twoHop? "two": "single"); fflush(stdout); };
			twoHop = size*3 < ns*4; size = ns; fg = levels.back().cg;
			if (ns*4 < 20*np) size = fg.size(); };
		  bisectSplitter<N> sp(levels.back().cg,np);
		  uvecD pt = sp.work(g.size(),verbose);					// build partition vector at lowest level
		  for (unsigned int l = levels.size()-1; l != ~0u; l--)			// refine partitioning at increasing resolution
			pt = levels[l].refine(pt,sp.parts,mins,maxs);
		  if (verbose) { for (const auto& p : sp.parts)
			printf("pt %d %e\n", p.pt, p.wt); };
		  return pt; }								// return partition vector at original resolution
};
#endif
