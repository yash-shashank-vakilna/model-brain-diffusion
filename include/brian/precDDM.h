#ifndef PRECDDM
#define PRECDDM

/*
 *
 * precDDM.h: implements domain-decomposition preconditioners for sparse linear systems
 * BRIAN Software Package Version 3.0
 *
 * $Id: precDDM.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (06/12/14): first implementation
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements domain-decomposition preconditioners for sparse linear systems.
*/

#define allCols(A, t)	(unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)
#define allEq(i)	(unsigned int i = 0; i < eq.N; i++)
#define allBlocks(b)	(unsigned int b = 0; b < bl.size(); b++)
#define allDOF(d)	(unsigned int d = 0; d < dof; d++)

#include <graphPartitioner.h>

//! Implements a "strength-of-aggregation" additive Schwarz preconditioner.

template<typename T, typename U = T> class precSAS : public precS<T,U> {
	const matS<U>& A;			//!< ref to LHS
	std::vector<matD<U>> Ai;		//!< block inverses
	uvecD	off;				//!< block offsets (N+1)
	uvecD	row;				//!< maps row to block
	uvecD	id;				//!< maps row to block
	matD<U>	getBlock(const unsigned int s, const unsigned int l)			//! returns inverse sub-matrix A(s:s+l,s:s+l).
		{ matD<U> S(l,l); S = U(0);						// allocate sub-matrix
		  for (unsigned int i = 0; i < l; i++) { const unsigned int r = row(s+i); // fill sub-matrix
			for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
				const unsigned int c = A.ci[t];
		  		for (unsigned int j = 0; j < l; j++) {
					if (c != row(s+j)) continue;
					S(i,j) = A.x[t]; break; } } };
		  return inv(S); }							// return inverse
	void	init(const U th)							//! computes block inverses.
		{ const ivecD agg = A.strength(th).aggregate();				// build aggregates
		  const unsigned int na = agg(agg.amax())+1; Ai.resize(na);		// find number of blocks
		  off.resize(na+1); off[0] = 0; uvecD bs(na); bs = 0;			// allocate offset array
		  for (unsigned int i = 0; i < A.M; i++) if (agg(i) >= 0) bs[agg(i)]++; // count rows per blocks
		  for (unsigned int i = 0; i < na; i++) off[i+1] = off[i]+bs[i];	// set block offsets
		  for (unsigned int i = 0; i < na; i++) bs[i] = off[i];			// copy offsets
		  for (unsigned int i = 0; i < A.M; i++) { if (agg(i) < 0) continue;	// for all aggregates
			const unsigned int r = bs[agg(i)]++; row[r] = i; }		// set aggregate to row map
		  for (unsigned int i = 0; i < Ai.size(); i++) {
			const unsigned int s = off(i), l = off(i+1)-s; assert(l > 0);	// get inverses
		  	Ai[i] = getBlock(s,l); } }
public:
	precSAS(const matS<U>& _A, const T th)
		//! constructs a Schwarz preconditioner with strong aggregation.
		/*! Parameters are: A refers to the stiffness matrix, th to the
		    aggregation threshold (see also precAMG).
		 */
		 : precS<T,U>(), A(_A), Ai(), off(), row(A.M), id(A.M)
		{ init(th); checkFpException(); }
        vecD<U>	operator()(const vecD<U>& b) const					//! applies Schwarz preconditioner to vector b.
		{ vecD<U> x(A.M); x = U(0);
		  for (unsigned int i = 0; i < Ai.size(); i++) {			// loop over the subdomains
			const unsigned int s = off(i), l = Ai[i].M; vecD<T> v(l);
			for (unsigned int r = 0; r < l; r++) v[r] = b(row(r+s));	// set RHS
			const vecD<T> xi = Ai[i]*v;					// multiply with block inverse
			for (unsigned int r = 0; r < l; r++) x[row(r+s)] += xi(r); };	// add to x
		 return x; }
};

//! Implements a restrictive additive Schwarz preconditioner.

template<typename T, typename U = T> class precRAS : public precS<T,U> {

//! Implements a block in a restrictive additive Schwarz preconditioner.

	class block {
		T eps = T(1e-12);		//!< limit for dropping elements.
		uvecD	eq;			//!< maps global to local equations in this block (restriction operator).
		LDL<T,U> Ai;			//!< holds LDL decomposition of this block.

		matS<U>	localA(const matS<U>& A) const					//! assembles sub-domain matrix.
			{ if (eq.N == 0) return matS<U>();
		  	  vecD<U> row(A.N); row = U(0); unsigned int nz = 0, k = 0;
		  	  for allEq(i) { const unsigned int r = eq(i);			// 1st pass: determine non-zeros
				for allCols(A,t) row[A.ci[t]] = A.x[t];			// map out to full row (is faster than search)
		  		for allEq(j) if (std::abs(row[eq(j)]) >= eps) nz++;
				for allCols(A,t) row[A.ci[t]] = U(0); };
		  	  matS<U> Al(eq.N,eq.N,nz); Al.ri[0] = 0;			// allocate matrix
		 	  for allEq(i) { const unsigned int r = eq(i);			// 2nd pass: fill-in
				for allCols(A,t) row[A.ci[t]] = A.x[t];
		  		for allEq(j) { const U x = row[eq(j)];
					if (std::abs(x) >= eps) {
						Al.ci[k] = j; Al.x[k] = x; k++; } };
				for allCols(A,t) row[A.ci[t]] = U(0);
				Al.ri[i+1] = k; };
		  	assert(k == nz); return Al; }					// check fill complete
	public:
		block() : eq(), Ai() { }
		block(const matS<U>& _A, const uvecD& _eq)				//! allocates preconditioner block.
			: eq(_eq), Ai(localA(_A)) { }
		void	solve(vecD<U>& x, const vecD<U>& b, const vecD<T>& wt) const	//! applies this block to b, returns x.
			{ if (eq.N == 0) return; vecD<U> bi(eq.N);
			  for allEq(i) bi[i] = b(eq(i));				// fill local RHS
			  const vecD<U> xi = Ai.solve(bi);				// solve block
			  for allEq(i) x[eq(i)] += xi(i)*wt(eq(i)); }			// add to global solution
	};

	std::vector<block> bl;			//!< vector of preconditioner blocks.
	vecD<T>	wt;				//!< equation weights.
	uvecD	getEquations(const volumeMesh& m, const uvecD& pt, const unsigned int b, const unsigned int dof)
		//! uses partitioning pt to collect interior and boundary nodes of block b.
		{ std::vector<unsigned int> nd;
		  for (unsigned int i = 0; i < pt.N; i++) { if (pt(i) != b) continue;	// for all nodes in this partition
			nd.push_back(i); vertex* v = m.vertices(i); assert(v);		// add node to this block
			for (const auto& e : v->edges) { 				// for all neighbors of this node
				if (pt(e.id) != b) nd.push_back(e.id); } };		// add neighbor if outside this partition
		  std::sort(nd.begin(),nd.end()); 					// sort 
		  nd.erase(std::unique(nd.begin(),nd.end()),nd.end()); 			// find and remove duplicates
		  uvecD eq(nd.size()*dof); unsigned int j = 0;				// allocate equation vector
		  for (const auto n : nd) { for allDOF(d) eq[j++] = n*dof+d; };		// convert inner nodes into eq#
		  return eq; }
public:
	precRAS(const matS<U>& A, const volumeMesh& m, const unsigned int dof, const unsigned int nb, const bool verbose = false)
		//! constructs a restrictive additive overlapping Schwarz preconditioner.
		/*! Parameters are: A refers to the stiffness matrix (assumed to be symmetric),
		    m to the volumetric mesh, dof to the degrees of freedom per node, and
		    nb to the number of blocks.
		 */
		: precS<T,U>(), bl(nb), wt(A.M)
		{ const uvecD pt = graphPartitioner<vertex>(m.vertices,nb).work();	// compute partitioning
		  wt = T(0);
#pragma omp parallel for
		  for allBlocks(b) {							// for all blocks
			if (verbose) { printf("prec %3.0f%% %zd\r", 100.0*b/nb, bl.size()); fflush(stdout); };
			const uvecD eq = getEquations(m,pt,b,dof);			// get eq# for this block
			for allEq(i) wt[eq(i)] += T(1);					// count occurrence of this eq
			bl[b] = block(A,eq); };						// allocate block
		  for (unsigned int i = 0; i < wt.N; i++) wt[i] = T(1)/wt(i); }		// set weight as inverse of occurrence
	vecD<U>	operator()(const vecD<U>& b) const					//! applies preconditioner to vector b.
		{ vecD<U> x(b.N); x = U(0);
#pragma omp parallel for
		  for allBlocks(i) bl[i].solve(x,b,wt);
		  return x; }
};

//! Implements a "Generic Elliptic Operator" (GenEO) preconditioner.

template<typename T, typename U = T> class precGenEO : public precS<T,U> {

//! Implements a block in a "Generic Elliptic Operator" (GenEO) preconditioner.

	using vNodes = std::vector<unsigned int>;

	class block {
	public:
		T	eps = T(1e-12);		//!< limit for dropping elements.
		unsigned int ne;		//!< number of equations in this block.
		unsigned int bd;		//!< boundary equations start at this index.
		uvecD	eq;			//!< vector of global equation numbers.
		LDL<T,U> Ai;			//!< holds LDL decomposition of this block.
		matD<U>	Wi;			//!< local corrector matrix.

		block() : ne(0), bd(0), eq(), Ai(), Wi() { }
		block(const unsigned int dof, const vNodes& in, const vNodes& bd)	//! allocate a block for the GenEO preconditioner.
			: ne((in.size()+bd.size())*dof), bd(in.size()*dof), eq(ne), Ai(), Wi()
			{ unsigned int j = 0;
			  for (const auto i : in) { for allDOF(d) eq[j++] = i*dof+d; };	// convert nodes into eq#
			  for (const auto i : bd) { for allDOF(d) eq[j++] = i*dof+d; } }
		matS<U>	localA(const matS<U>& A) const					//! assembles sub-domain matrix.
			{ vecD<U> row(A.N); row = U(0); unsigned int nz = 0, k = 0;
			  for allEq(i) { const unsigned int r = eq(i);			// 1st pass: determine non-zeros
				for allCols(A,t) row[A.ci[t]] = A.x[t];			// map out to full row (faster than search)
			  	for allEq(j) if (std::abs(row[eq(j)]) >= eps) nz++;
				for allCols(A,t) row[A.ci[t]] = U(0); };
			  matS<U> Al(eq.N,eq.N,nz); Al.ri[0] = 0;			// allocate matrix
			  for allEq(i) { const unsigned int r = eq(i);			// 2nd pass: fill-in
				for allCols(A,t) row[A.ci[t]] = A.x[t];
			  	for allEq(j) { const U x = row[eq(j)];
					if (std::abs(x) >= eps) { Al.ci[k] = j; Al.x[k] = x; k++; } };
				for allCols(A,t) row[A.ci[t]] = U(0);
				Al.ri[i+1] = k; };
			  assert(k == nz); return Al; }					// check fill complete
		matS<U>	localB(const matS<U>& Al, const vecD<T>& wt) const		//! assembles sub-domain matrix B.
			{ matS<U> Bl = Al;
			  for allEq(r) for allCols(Bl,t) {
				Bl.x[t] *= wt(eq(r))*wt(eq(Bl.ci[t])); };
			  return Bl; }
		matD<U>	correction(const matS<U>& Al, const vecD<T>& wt) const		//! computes the coarse space correction.
			{ const matS<U> Bl = localB(Al,wt);
			  const unsigned int m = 4; const underscore _;
			  const auto D = sgsev<T>(Al,Bl,m,m*4,SMALLEST_MAGN).solve(T(100.0*eps));	// solve generalized eigenproblem
			  for (unsigned int j = 0; j < m; j++) {			// multiply partition-of-unity
				for allEq(i) D.U(i,j) *= wt(eq(i)); };
			  return D.U(_,_(0,m)); }
		void	initMatrices(const matS<U>& A, const vecD<T>& wt)		//! initializes block matrices
			{ if (ne == 0) return; const matS<U> Al = localA(A); 		// skip empty block
			  Ai = LDL<T,U>(Al); Wi = correction(Al,wt); }
		void	apply(vecD<U>& x, const vecD<U>& b, const vecD<T>& wt) const	//! applies this block to b, returns x.
			{ if (ne == 0) return; vecD<U> bi(eq.N); bi = U(0);
			  for allEq(i) bi[i] = b(eq(i));				// fill local RHS
			  const vecD<U> xi = Ai.solve(bi);				// solve block
			  for allEq(i) x[eq(i)] += xi(i)*wt(eq(i)); }			// add to global solution
	};

	std::vector<block> bl;			//!< vector of blocks
	vecD<T>	wt;				//!< partition-of-unity vector
	matS<U>	Z, Zt;				//!< coarse-to-full mapping
	matD<U> Ei;				//!< coarse-space correction matrix

	vNodes	findInteriorNodes(const uvecD& pt, const unsigned int p) const
		//! uses partitioning pt to collect interior nodes of partition p.
		{ vNodes in; for (unsigned int i = 0; i < pt.N; i++) {
			if (pt(i) == p) in.push_back(i); }; return in; }
	vNodes	findBoundaryNodes(const volumeMesh& m, const uvecD& pt, const unsigned int p, const vNodes& in) const					
		//! takes mesh and a list of interior nodes - returns a list of boundary vertices.
		{ vNodes bd; for (const auto id : in) {					// for all interior vertices
			const vertex* v = m.vertices(id); assert(v);
			for (const auto& e : v->edges) { 				// for all neighbors
				if (pt(e.id) != p) bd.push_back(e.id); } };		// add as boundary if not in this partition
		  std::sort(bd.begin(),bd.end()); 					// sort 
		  bd.erase(std::unique(bd.begin(),bd.end()),bd.end()); return bd; }	// find and remove duplicates
	void	initBlocks(const matS<U>& A, const volumeMesh& m, const uvecD& pt,
			const unsigned int dof, const bool verbose)			//! initializes preconditioner blocks.
		{ wt = T(0); for allBlocks(b) {						// for all blocks
			const vNodes in = findInteriorNodes(pt,b);			// find interior nodes
			if (in.size() == 0) continue;					// skip if none
			const vNodes bd = findBoundaryNodes(m,pt,b,in); 		// add boundary nodes
			bl[b] = block(dof,in,bd);					// allocate block & init equations
			const uvecD& eq = bl[b].eq; for allEq(i) wt[eq(i)] += T(1); };	// count occurrence of eqs in this block
		  for (unsigned int i = 0; i < wt.N; i++) wt[i] = T(1)/wt(i);		// set partition-of-unity vector
		  for allBlocks(b) { bl[b].initMatrices(A,wt);				// set up all block matrices Ai,Wi
			if (verbose) { printf("prec %3.0f%%\r", 100.0*(b+1)/bl.size());	fflush(stdout); } } };
	void	initCorrection(const matS<U>& A)					//! computes correction matrices Z,Zi,Ei.
		{ unsigned int m = 0, nz = 0, n = wt.N;					// determine space for global projector Zt
		  for (const auto& b : bl) { m += b.Wi.N; nz += b.Wi.nel(); }		// for all blocks: sample dimensions of local projector
		  Zt.resize(m,n,nz); Zt.ri[0] = 0; unsigned int t = 0, r = 0;		// allocate matrix
		  for (const auto& b : bl) {
			const matD<U>& W = b.Wi; const uvecD& eq = b.eq;		// get local to global mapping
		  	for (unsigned int j = 0; j < W.N; j++) {			// assemble Zt from local matrix
		  		for (unsigned int i = 0; i < W.M; i++) { 
					Zt.x[t] = W(i,j); Zt.ci[t] = eq(i); t++; };
				r++; Zt.ri[r] = t; } };
		  assert(nz == t); assert(r == m); Z = trp(Zt); 			// check fill complete
		  const matS<U> E = Zt*(A*Z); Ei = inv(asDense(E)); }			// and its inverse Ei
public:
	precGenEO(const matS<U>& A, const volumeMesh& m, const unsigned int dof, const unsigned int nb, const bool verbose = false)
		//! constructs a restrictive additive Schwarz preconditioner.
		: precS<T,U>(), bl(nb), wt(A.M), Z(), Zt(), Ei()
		{ const uvecD pt = graphPartitioner<vertex>(m.vertices,nb).work();	// compute partitioning
		  initBlocks(A,m,pt,dof,verbose);					// allocate local blocks and partition-of-unity vector
		  initCorrection(A); }							// compute correction matrices Ei,Z,Zt
        vecD<U>	operator()(const vecD<U>& b) const					//! applies Schwarz preconditioner to vector b.
		{ vecD<U> x = Z*(Ei*(Zt*b));
		  for allBlocks(i) bl[i].apply(x,b,wt);
		  return x; }
};
#undef allCols
#undef allEq
#undef allBlocks

#endif
