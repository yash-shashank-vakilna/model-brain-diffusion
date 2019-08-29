#ifndef GRAPH_H
#define GRAPH_H

/*
 *
 * graph.h: graph structure
 * BRIAN Software Package Version 3.0
 *
 * $Id: graph.h 506 2017-03-19 01:39:21Z frithjof $
 *
 * 0.10 (12/12/00): initial version
 * 0.20 (08/11/09): rewritten for BRIAN2
 * 0.30 (10/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v400 (20/09/16): reimplemented for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for a generic graph structure.
 */

//! Implements a generic graph

template<typename N> class graph {
protected:
	std::vector<N*> nd;			//!< vector of nodes in this graph (must be pointers...)
	attrList at;				//!< list of other attributes
	bool	hasWeights;			//!< true if weights are used

	void 	assign(const graph& g)							//! copies a graph.
		{ clear(); nd.resize(g.size()); at = g.at; hasWeights = g.hasWeights;
		  const unsigned int nv = size();
		  for (unsigned int i = 0; i < nv; i++)
			nd[i] = g.nd[i]? g.nd[i]->clone(): nullptr; }
public:
	graph(const unsigned int n = 0, const bool w = false)				//! allocates a graph for n nodes.
		 : nd(n,nullptr), at(), hasWeights(w) { }
	graph(const graph& g)								//! copies from graph g.
		 : nd(), at(), hasWeights() { assign(g); }
	graph(graph&& g)
		 : nd(std::move(g.nd)), at(std::move(g.at)), hasWeights(g.hasWeights)
		{ std::vector<N*>().swap(g.nd); }					// empties node list in g.
	graph(FILE* fp)
		: nd(), at(), hasWeights() { read(fp); }
	virtual ~graph() { clear(); }
	graph& operator=(const graph& g)						//! assigns from graph g.
		{ if (this != &g) assign(g); return *this; }
	graph& operator=(graph&& g)							//! move assigns from graph g.
		{ assert(this != &g); clear();
		  nd = std::move(g.nd); at = std::move(g.at); hasWeights = g.hasWeights;
		  return *this; }
	N*	operator()(const unsigned int i) const					//! returns node i.
		{ return nd[i]; }
	unsigned int size() const							//! returns number of last used node position.
		{ return nd.size(); }
	const attrList& attributes() const						//! returns a graph's attribute list (for reading).
		{ return at; }
	attrList& attributes()								//! returns a graph's attribute list (for writing).
		{ return at; }
	void	copyAttributes(const attrList& list)					//! copies attribute list from list.
		{ at = list; }
	void 	clearNodes()								//! clears all nodes in graph.
		{ const unsigned int nv = size();
		  for (unsigned int i = 0; i < nv; i++) {				// delete all node pointers
			delete nd[i]; nd[i] = nullptr; };
		  nd.clear(); }								// reset size information
	unsigned int countNodes() const							//! returns count of nodes in this graph.
		{ const unsigned int nv = size(); unsigned int n = 0;
		  for (unsigned int i = 0; i < nv; i++) { if (nd[i]) n++; };
		  return n; }								// reset size information
	void 	clear()									//! clears all data in graph.
		{ at.clear(); clearNodes(); }
	bool	isCompact() const							//! returns true if all nodes are filled.
		{ for (const auto n : nd) if (n == nullptr) return false;
		  return true; }
	bool	nextNode(unsigned int& s) const						//! returns id of next allocated node.
		{ const unsigned int nv = size();
		  for (unsigned int i = s; i < nv; i++)
			if (nd[i]) { s = i; return true; }
		  return false; }
	unsigned int lastNode() const							//! returns id of the last node in this graph.
		{ for (unsigned int i = size()-1; i != ~0u ; i--) {
			if (nd[i]) return i+1; };
		  return 0; }								// reset size information
	bool	hasEdgeUni(unsigned int a, unsigned int b) const			//! checks for unidirectional edge from a to b.
		{ return nd[a]? nd[a]->hasEdge(b): false; }
	bool	hasEdgeBi(unsigned int a, unsigned int b) const				//! checks for bidirectional edge between a and b.
		{ return hasEdgeUni(a,b) && hasEdgeUni(b,a); }
	void	linkNodesUni(unsigned int a, unsigned int b, float w = 0.0f) 		//! makes a unidirectional edge from a to b.
		{ if (nd[a] && nd[a]->hasEdge(b) == false) nd[a]->addEdge({b,w}); }
	void	linkNodesBi(unsigned int a, unsigned int b, float w = 0.0f) 		//! makes a bidirectional edge between a and b.
		{ linkNodesUni(a,b,w); linkNodesUni(b,a,w); }
	void	unlinkNodesUni(unsigned int a, unsigned int b)				//! removes edge from a to b.
		{ if (nd[a]) nd[a]->unlink(b); }
	void	unlinkNodesBi(unsigned int a, unsigned int b)				//! removes edges between a and b.
		{ unlinkNodesUni(a,b); unlinkNodesUni(b,a); }
	void	unlinkAll(unsigned int a)						//! removes all links of node a.
		{ if (nd[a] == nullptr) return;
		  for (const auto& e : nd[a]->edges) unlinkNodesUni(e.id,a);
		  nd[a]->clear(); }
	unsigned int addNode(N* n)							//! adds node to graph. Always resizes container.
		{ nd.push_back(n); return size()-1; }
	void	addNodeAt(N* n, const unsigned int i)					//! inserts node n at index i. May resize container.
		{ if (i >= size()) { nd.resize(i+1,nullptr); }; 
		  delete nd[i]; nd[i] = n; }
	void	clearWeights()								//! clears weights in all nodes and edges.
		{ for (auto n : nd) if (n) n->clearWeights(); }
	float	totalDegree() const							//! returns total node weight.
		{ float s = 0.0f; for (const auto n : nd) if (n) s += n->totalDegree();
		  return s; }
	unsigned int nEdges() const 							//! returns total number of edges.
		{ unsigned int s = 0;
		  for (const auto n : nd) { if (n) s += n->nEdges(); };
		  return s; }
	bool	removeNode(unsigned int i)						//! removes node i from graph including edges.
		{ if (nd[i] == nullptr) return false;
		  unlinkAll(i); delete nd[i]; nd[i] = nullptr; return true; }
	void	resize(const unsigned int n)						//! resizes graph for n nodes.
		{ nd.resize(n,nullptr); }
	void	print()	const								//! prints information about this graph.
		{ const unsigned int nv = size();
		  for (unsigned int i = 0; i < nv; i++)
			if (nd[i]) { printf("%6d: ", i); nd[i]->print(hasWeights); }; }
	nodeID	nodeType() const							//! finds node type from first allocated node.
		{ unsigned int s = 0;
		  if (nextNode(s) == false) throw rtException("Graph is empty");
		  return nd[s]->getID(); }
	void	read(const char* fname)							//! reads vector graph from file named fname.
		{ FILE *fp = openFile(fname,"r"); read(fp); closeFile(fp); }
	void	save(const char* fname)	const						//! saves vector graph in file named fname.
		{ FILE *fp = openFile(fname,"w"); save(fp); closeFile(fp); }
	void	read(FILE* fp)								//! reads the first graph from file.
		{ bundleList bl; readFile(fp,bl); fromBundle(*bl.begin()); }
	void	save(FILE* fp, const char* tag = "vertices") const			//! saves graph into file.
		{ bundleList bl; bl.push_back(toBundle(tag)); writeFile(fp,bl); }
	void	fromBundle(const bundle& b);						// non-generic, see graph.C
	bundle	toBundle(const char* tag) const;					// generic, see below
	uvecD	compact();
	bool	testLinks(bool verbose = false);
	uvecD	label(unsigned int &lmax) const;
	uvecD	label(const bvecD& ind) const;
	bool	resizeForValues();							// specialization for vertex graphs
	bool	resizeForNormals();
	bool	resizeForColors();
};

template<typename N>
bundle graph<N>::toBundle(const char* tag) const
//! saves a graph in bundle b.
{	assert(size()); os out;	char buf[LINELEN]; 					// allocate stream buffer
	for (unsigned int i = 0; i < lastNode(); i++) {					// for all nodes
		const node* n = nd[i]; if (n == nullptr) continue;			// skip empty ids
		out.save(i); out.save(n->getID()); n->save(out,hasWeights); };		// save id, node type and node fields
	bundle b; b.copyStream(out.data()); b.copyAttributes(at);			// copy buffer and attributes to bundle
	b.repn = repnType::graph; b.updateValue(tag);					// add tag for graph type
	sprintf(buf, "%ld", 0L); b.at.update("data", buf);
	sprintf(buf, "%zd", b.length); b.at.update("length", buf);
	sprintf(buf, "%u", lastNode()); b.at.update("nodes", buf);
	if (hasWeights) b.at.update("hasWeights", "true");
	return b;
}

template<typename N>
bool graph<N>::testLinks(const bool verbose)
//! performs connectivity tests on a graph, returns false if isolated nodes were found.
{	bool succ = true;
	for (unsigned int i = 0; i < size(); i++) {					// for all nodes
		node* n = nd[i]; if (n == nullptr) continue;
repeat:		for (const auto& e : n->edges) {
			if (e.id == i) { n->unlink(i); 					// self-reference
				if (verbose) printf("Removed link to self at node %u.\n", i);
				goto repeat; }
			else { node* p = nullptr;
				if (e.id < size()) p = nd[e.id];
				if (p == nullptr) { n->unlink(e.id); 			// non-existing neighbor
					if (verbose) printf("Removed link from node %u to non-existing neighbor %u.\n", i,e.id);
					goto repeat; }
				else if (p->hasEdge(i) == false) { p->addEdge({i});	// check bidirectional link
					if (verbose) printf("Linked node %u to node %u.\n", i,e.id);
					goto repeat; } } } };
	for (unsigned int i = 0; i < size(); i++) {					// for all nodes
		if (nd[i] == nullptr || nd[i]->nEdges()) continue;
		succ = false; if (verbose) printf("Isolated node %u.\n", i); }		// report isolated nodes
	return succ;
}

template<typename N>
uvecD graph<N>::compact()
//! moves nodes in table so that all entries are consecutive.
{	uvecD map(size()); map = ~0u; unsigned int j = 0;
	for (unsigned int i = 0; i < size(); i++) { if (nd[i] == nullptr) continue;	// compact pointer table
		if (i != j) { nd[j] = nd[i]; nd[i] = nullptr; }; map[i] = j++; };
	nd.resize(j);
	for (auto& n : nd) { for (auto& e : n->edges) e.id = map(e.id); };		// reassign links
	return map;
}

template<typename N>
uvecD graph<N>::label(unsigned int& lmax) const
//! labels connected components based on node connectivity.
{	uvecD lbl(size()); lbl = 0; std::stack<unsigned int> st; unsigned int li = 0; 
	for (unsigned int i = 0; i < size(); i++) {
		if (nd[i] == nullptr || lbl(i)) continue;	
		st.push(i); lbl(i) = ++li;						// new component, push onto stack
		while (!st.empty()) { const unsigned int j = st.top(); st.pop();	// for all unvisited vertices
			if (nd[j] == nullptr) continue;
			for (const auto& e : nd[j]->edges) { const unsigned int k = e.id; // for all unvisited neighbors
				if (nd[k] == nullptr || lbl(k)) continue;		// not marked or already labeled
				st.push(k); lbl(k) = li; } } };				// push onto stack and label
	lmax = li; return lbl;
}

template<typename N>
uvecD graph<N>::label(const bvecD& ind) const
//! labels connected components based on boolean indicator ind.
{	uvecD lbl(size()); lbl = 0; std::stack<unsigned int> st; unsigned int li = 0; 
	for (unsigned int i = 0; i < size(); i++) {
		if (ind(i) == false || lbl(i)) continue;	
		st.push(i); lbl(i) = ++li;						// new component, push onto stack
		while (!st.empty()) { const unsigned int j = st.top(); st.pop();	// for all unvisited vertices
			if (nd[j] == nullptr) continue;
			for (const auto& e : nd[j]->edges) { const unsigned int k = e.id; // for all unvisited neighbors
				if (ind(k) == false || lbl(k)) continue;		// not marked or already labeled
				st.push(k); lbl(k) = li; } } };				// push onto stack and label
	return lbl;
}

//! Implements a graph variant based on a voxel grid.

template<typename N> class imageGraph : public graph<N> {
protected:
	unsigned int nx;			//!< grid extent x
	unsigned int ny;			//!< grid extent y
	unsigned int nz;			//!< grid extent z

public:
   	imageGraph()									//! allocates an empty vector graph.
		: graph<N>(), nx(0), ny(0), nz(0) { }
   	imageGraph(const fimage& src, const char* ip, const bool w = false)
		//! allocates a vector graph for image src, type *ip, and weight field w.
		: graph<N>(src.nel(),w), nx(src.nx), ny(src.ny), nz(src.nz)
		{ char buf[40]; const char *s;
		  sprintf(buf, "%u %u %u", nx, ny, nz);
		  graph<N>::attributes().add(attribute("extent", buf));
		  fvec3 v = 1.0f; src.getVoxelSize(v); sprintf(buf, "%.4f %.4f %.4f", v.x,v.y,v.z);
		  graph<N>::attributes().add(attribute("voxel", buf));
		  s = src.at.lookup("patient");
		  if (s) graph<N>::attributes().add(attribute("patient", s));
		  s = src.at.lookup("date");
		  if (s) graph<N>::attributes().add(attribute("date", s));
		  s = src.at.lookup("convention");
		  if (s) graph<N>::attributes().add(attribute("convention", s));
		  graph<N>::attributes().add(attribute("nodeType", ip)); }
   	imageGraph(const graph<N>& g) : graph<N>(g), nx(0), ny(0), nz(0)		//! copies from graph g.
		{ setExtent(); }
   	imageGraph(const uvec3& ex) : graph<N>(ex.x*ex.y*ex.z), nx(ex.x), ny(ex.y), nz(ex.z) //! allocates an empty imageGraph for extent ex.
		{ }
   	imageGraph(const imageGraph& g) : graph<N>(g), nx(g.nx), ny(g.ny), nz(g.nz) { }	//! copies from imageGraph g.
   	imageGraph(imageGraph&& g)
		 : graph<N>(std::move(g)), nx(g.nx), ny(g.ny), nz(g.nz) { }		//! moves from imageGraph g.
	imageGraph& operator=(const imageGraph& g)					//! copy assigns from imageGraph g.
		{ if (this != &g) { graph<N>::operator=(static_cast<graph<N> const&>(g));
		  	nx = g.nx; ny = g.ny; nz = g.nz; }; return *this; }
	imageGraph& operator=(imageGraph&& g)						//! move assigns from imageGraph g.
		{ assert(this != &g); graph<N>::operator=(static_cast<graph<N> const&>(g));
		  nx = g.nx; ny = g.ny; nz = g.nz; return *this; }
	uvec3	extent() const								//! returns extent of underlying image.
		{ return uvec3(nx,ny,nz); }
	void	setExtent()								//! sets extent attribute.
		{ const char* s = graph<N>::attributes().lookup("extent");
		  if (s) sscanf(s, "%u %u %u", &nx, &ny, &nz);
		  else throw rtException("Extent in imageGraph missing"); }		// check that dimensions match graph size
	unsigned int nel() const							//! returns number of sites in vector graph.
		{ return nx*ny*nz; }
	unsigned int index(const uvec3& s) const					//! returns index of site s.
		{ return (s.z*ny+s.y)*nx+s.x; }
	uvec3	site(const unsigned int i) const					//! returns site of index i.
		{ const unsigned int z = i/(nx*ny), r = i%(nx*ny);
		  return uvec3(r%nx,r/nx,z); }
	bool	getVoxelSize(fvec3& v) const						//! returns voxel size in read world dimensions.
		{ const char* s = graph<N>::attributes().lookup("voxel");
		  if (s) sscanf(s, "%f%f%f", &v.x, &v.y, &v.z);
		  return s != nullptr; }
	void	addNodeAt(N* n, const uvec3& s)						//! adds node n at site s.
		{ const unsigned int i = index(s); assert(i < graph<N>::size());
		  graph<N>::addNodeAt(n,i); }
	N*	operator()(const unsigned int i) const					//! returns node at site s.
		{ return i < graph<N>::size()? graph<N>::nd[i]: nullptr; }
	N*	operator()(const uvec3& s) const					//! returns node at site s.
		{ const unsigned int i = index(s); return imageGraph<N>::operator()(i); }
	void	fromBundle(const bundle& b)						// see graph.C
		{ graph<N>::fromBundle(b); setExtent(); }
	bundle	toBundle(const char* tag) const						// see graph.C
		{ bundle b = graph<N>::toBundle(tag);
		  char buf[LINELEN]; sprintf(buf, "%u %u %u", nx, ny, nz);
		  b.at.update("extent", buf); return b; }
	void	read(FILE* fp)
		{ bundleList bl; readFile(fp,bl);
		  if (bl.size() != 1) throw rtException("Invalid number of bundles");
		  const bundle& b = *bl.begin();
		  if (strcmp(b.value,"descriptors") == 0) { read(b); setExtent(); }	// calls specialized read() above
		  else throw rtException("Invalid node type"); }
	void	save(FILE* fp, const char* tag = "descriptors") const
		{ bundleList bl; bl.push_back(toBundle(tag)); writeFile(fp,bl); }
};

using descGraph = imageGraph<descriptor>;
using streamGraph = imageGraph<streamDesc>;
using tensorGraph = imageGraph<tensorDesc>;
using spharmGraph = imageGraph<spharmDesc>;
using fiberGraph = imageGraph<fiberMixtureDesc>;
using cylinderGraph = imageGraph<fiberZCDDesc>;


#endif
