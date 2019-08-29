#ifndef COMMUNITY_H
#define COMMUNITY_H

/*
 *
 * community.h: community building in graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: community.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (10/12/16): initial version FK
 * 0.20 (12/12/16): released for BRIAN3
 *
 */

#include <set>
#include <glpk.h>

#define allNodes(n)	(unsigned int n = 0; n < N; n++)
using uvec = std::vector<unsigned int>;

/*! \file
    \brief classes and methods community building in graphs represented as matrices.
*/

template<typename T> 
T modularity(const matD<T>& A, const uvecD& map)
//! returns modularity for mapping map in graph represented by matrix A.
{	const T wt = A.sum(); if (wt == T(0)) return T(0);
	vecD<T> wr = A.rowsum(), wc = A.colsum(); T q{0};
	for (unsigned int i = 0; i < A.M; i++)
		for (unsigned int j = 0; j < A.N; j++) {
			if (map(i) == map(j)) q += A(i,j)-wr(i)*wc(j)/wt; };
	return q/wt;
}

//! Implements an algorithm for finding communities in a bipartite graph.

template<typename T> class bipCommunities {
	uniformDist<T> ud;			//!< generates uniformly distributed numbers
	const matD<T> M;			//!< the affinity matrix
	const vecD<T> srows;			//!< row sums of the affinity matrix
	const vecD<T> scols;			//!< col sums of the affinity matrix
	const T smat;				//!< total sum of the affinity matrix
	const matD<T> B;			//!< Barber's matrix
	uvecD	red;				//!< labels for first graph
	uvecD	blue;				//!< labels for second graph
	T	qbest;				//!< best modularity so far

	unsigned int maxLabel(const uvecD& v) const					//! returns max. label+1 in vector v.
		{ return v(v.imax())+1; }
	bool	findLabel(const uvecD& v, const unsigned int l) const			//! returns true if vector v contains label l.
		{ for (unsigned int i = 0; i < v.N; i++) if (v(i) == l) return true;
		  return false; }
	matD<T>	indMatrix(const uvecD& v) const						//! returns indicator matrix for label assignment v.
		{ uvecD u = unique(v); matD<T> L(u.N,v.N); L = T(0);
		  if (maxLabel(v) == 0) return L;
		  for (unsigned int i = 0; i < u.N; i++)
		  	for (unsigned int j = 0; j < v.N; j++)
				L(i,j) = u(i) == v(j);
		  return L; }
	T	modularity(const uvecD& vr, const uvecD& vb) const			//! returns modularity for label assignments (vr, vb).
		{ const matD<T> M = indMatrix(vr)*B*trp(indMatrix(vb));
		  return M.d().sum()/smat; }
	uvec	division() const							//! returns division for current label assignment.
		{ uvec div; const unsigned int n = std::min(red.N,blue.N);
		  for (unsigned int i = 0; i < n; i++) {
			if (findLabel(red,i) && findLabel(blue,i))
				div.push_back(i); }
		  return div; }
	T	checkDiv(const unsigned int d1, const unsigned int d2) const		//! returns modularity for division (d1,d2).
		{ uvecD r = red, b = blue;
		  for (unsigned int i = 0; i < red.N; i++) if (red(i) == d1) r[i] = d2;
		  for (unsigned int j = 0; j < blue.N; j++) if (blue(j) == d1) b[j] = d2;
		  return modularity(r,b); }
	bool	bestDiv(const uvec& div, const unsigned int i, const unsigned int j) const //! returns true if division (i,j) is best.
		{ const T q = checkDiv(div[i],div[j]); if (q <= qbest) return false;
		  for (const auto d : div) {
			if (checkDiv(d,div[i]) > q) return false;
			if (checkDiv(d,div[j]) > q) return false; };
		  return true; }
	void	init(unsigned int m);
	vecD<T>	initRD() const;
	vecD<T>	initBD() const;
	void	updateRD(vecD<T>& rd, const vecD<T>& bd);
	void	updateBD(const vecD<T>& rd, vecD<T>& bd);
	void	propagateLabels();
	void	agglomerateLabels(const bool verbose);
public:
	bipCommunities(const matD<T>& _M)						//! allocates a context for bipartite graph matching.
		: ud(), M(_M), srows(M.rowsum()), scols(M.colsum()), smat(M.sum()),
		B(M-outer(srows,scols)/smat), red(), blue(), qbest(std::numeric_limits<T>::lowest()) { }
	T	work(uvecD& _red, uvecD& _blue, const unsigned int m = 0,
			const bool verbose = false)
		{ init(m); propagateLabels(); agglomerateLabels(verbose);
		   _red = red; _blue = blue; return qbest; }
};

//! \relates bipCommunities
template <typename T> 
void bipCommunities<T>::init(unsigned int m)
//! initializes label assignments.
{	red.resize(M.M); blue.resize(M.N); m = std::min(m,M.M);
	qbest = std::numeric_limits<T>::lowest();
	for (unsigned int i = 0; i < red.N; i++)
		red[i] = m? static_cast<unsigned int>(std::floor(ud()*m)): i;
	for (unsigned int i = 0; i < blue.N; i++) blue[i] = -1;
}

//! \relates bipCommunities
template <typename T> 
vecD<T> bipCommunities<T>::initBD() const
//! initializes blue label degrees.
{	vecD<T> bd(blue.N); bd = T(0); if (maxLabel(blue) == 0) return bd;
	for (unsigned int j = 0; j < bd.N; j++) {
		if (blue(j) > bd.N) bd[j] = scols(j); else bd[j] += scols(j); };
	return bd;
}

//! \relates bipCommunities
template <typename T> 
void bipCommunities<T>::updateBD(const vecD<T>& rd, vecD<T>& bd)
//! updates blue label degrees.
{	const unsigned int n = maxLabel(red);
	for (unsigned int j = 0; j < blue.N; j++) {
		unsigned int b = blue(j); if (b < bd.N) bd[b] -= scols(j); 
		T dmax = std::numeric_limits<T>::lowest(); uvec ix;
		for (unsigned int s = 0; s < n; s++) { T d = -scols(j)*rd(s)/smat;
		 	for (unsigned int i = 0; i < red.N; i++)
				if (red(i) == s) d += M(i,j);
			if (d < dmax) continue;
			if (d > dmax) { dmax = d; ix.clear(); }; ix.push_back(s); };
		const unsigned int k = static_cast<unsigned int>(std::floor(ud()*ix.size()));
		b = ix[k]; blue[j] = b; bd[b] += scols(j); }
}

//! \relates bipCommunities
template <typename T> 
vecD<T> bipCommunities<T>::initRD() const
//! initializes red label degrees.
{	vecD<T> rd = srows; return rd;
}

//! \relates bipCommunities
template <typename T> 
void bipCommunities<T>::updateRD(vecD<T>& rd, const vecD<T>& bd)
//! updates red label degrees.
{	const unsigned int n = maxLabel(blue);
	for (unsigned int i = 0; i < red.N; i++) {
		unsigned int r = red(i); rd[r] -= srows(i); 
		T dmax = std::numeric_limits<T>::lowest(); uvec ix;
		for (unsigned int s = 0; s < n; s++) { T d = -srows(i)*bd(s)/smat;
		 	for (unsigned int j = 0; j < blue.N; j++)
				if (blue(j) == s) d += M(i,j);
			if (d < dmax) continue;
			if (d > dmax) { dmax = d; ix.clear(); }; ix.push_back(s); };
		const unsigned int k = static_cast<unsigned int>(std::floor(ud()*ix.size()));
		r = ix[k]; red[i] = r; rd[r] += srows(i); }
}

//! \relates bipCommunities
template <typename T> 
void bipCommunities<T>::propagateLabels()
//! step one: performs label propagation.
{	vecD<T> rd = initRD(), bd = initBD(); qbest = modularity(red,blue);
	while (true) { uvecD r = red, b = blue; vecD<T> ord = rd, obd = bd;
		updateBD(rd,bd); updateRD(rd,bd); T qn = modularity(red,blue);		// at the first step, qn may be < 0
		if (qn > T(0) && qn <= qbest) {						// so do not jump out right here
			red = r; blue = b; rd = ord; bd = obd; break; };
		qbest = qn; };
}

//! \relates bipCommunities
template <typename T> 
void bipCommunities<T>::agglomerateLabels(const bool verbose)
//! step one: performs label agglomeration.
{	uvec div = division();
	for (bool flag = true; flag; ) { 
		if (div.size() <= 1) flag = false;
		else {	unsigned int n = 0;
			for (unsigned int i = 0; i < div.size()-1; i++) {
				for (unsigned int j = i+1; j < div.size(); j++) {
					if (verbose) { printf("[%3d %3d] %6.4f\r", i,j,qbest); fflush(stdout); }
					if (bestDiv(div,i,j) == false) continue;
					for (unsigned int k = 0; k < red.N; k++)
						if (red(k) == div[i]) red[k] = div[j];
					for (unsigned int k = 0; k < blue.N; k++)
						if (blue(k) == div[i]) blue[k] = div[j];
					n++; } }
			if (verbose) printf("\n");
			if (n == 0) flag = false; };
		propagateLabels(); div = division(); };
}

/*! \brief finds communities in a bipartite graphs.
 */

/*!
\param[out] red returns label assignments of the first graph
\param[out] blue returns label assignments of the second graph
\param[in] A references the affinity matrix of both graphs
\param[in] mini corresponds to the minimal number of modules to consider
\param[in] reps corresponds to the number of iterations per additional module
\return the best modularity value
 */
//! \relates bipCommunities
template <typename T> 
T bipComm(uvecD& red, uvecD& blue, const matD<T>& A, const unsigned int mini = 0,
	const bool verbose = false)
{	const bool flip = A.M > A.N; const unsigned int nit = 10;
	bipCommunities<T> lp(flip? trp(A): A); T qbest = lp.work(red,blue,0,verbose);
	if (mini) { const uvecD u = unique(red); const unsigned int nm = u.N;
		for (unsigned int m = 1; m <= mini; m++) {
			for (unsigned int it = 0; it < nit; it++) {
				uvecD r = red, b = blue; T q = lp.work(r,b,m+nm,false);
				if (q > qbest) { red = r; blue = b; qbest = q; } };
			if (verbose) { printf("[+%3d] %6.4f\r", m, qbest); fflush(stdout); } } };
	if (flip) std::swap(red,blue);
	if (verbose) { printf("final q %6.4f for %u communities.\n", qbest, unique(red).N); };
	return qbest;
}

//! Implements an algorithm for finding communities in a graph using label propagation.

/*
 * refer to:
 * Campigotto R, Cespedes PC, Guillaume JL (2014)
 * A Generalized and Adaptive Method for Community Detection
 * arVix:1406.2518v1
 *
 */

template<typename T> class lpropCommunities {
	const unsigned int maxit = 1000;	//!< maximum number of iterations in aggregation
	matD<T>& A;				//!< input graph
	unsigned int N;				//!< number of nodes in the graph
	uvecD	n2c;				//!< maps nodes to communities
	vecD<T>	in, tot; 			//!< used to compute modularity
	T	gw;				//!< total weigth of graph

	void	remove(const unsigned int i, const unsigned int c, const T dwc)		//! removes node i from community c.
		{ in[c] -= dwc+dwc+A(i,i); tot[c] -= deg(i); n2c[i] = ~0u; }
	void	insert(const unsigned int i, const unsigned int c, const T dwc)		//! inserts node i in community c.
		{ in[c] += dwc+dwc+A(i,i); tot[c] += deg(i); n2c[i] = c; }
	T	mod(const unsigned int c) const						//! returns modularity for community c.
		{ return tot(c) > T(0)? in(c)/gw-SQR(tot(c)/gw): T(0); }
	T	deg(const unsigned int i) const						//! returns weighted degree for node i.
		{ T w{0}; for allNodes(j) w += A(i,j); return w; }
	unsigned int moveNodes(const uvecD& rn);
public:
	lpropCommunities(matD<T>& _A)
		: A(_A), N(A.M), n2c(N), in(N), tot(N), gw(T(0))
		{  for allNodes(i) { n2c[i] = i; in[i] = A(i,i);
			tot[i] = deg(i); gw += tot[i]; } }
	T	modularity() const							//! returns the modularity of the current partitioning.
		{ T q{0}; for allNodes(i) q += mod(i); return q; };
	matD<T>	partition(uvecD& map) const;
	bool	aggregate(const T eps = T(1e-6));
	void	detMove(const uvecD& map);
};

template <typename T>
unsigned int lpropCommunities<T>::moveNodes(const uvecD& rn)
//! move nodes between communities; returns number of moves.
{	vecD<T> dw(N); dw = T(-1); uvecD pos(N); pos = 0; unsigned int nc = 0, nb = N;
	for (unsigned int r = 0; r < rn.N; r++) { const unsigned int i = rn(r);		// for all nodes i
		for (unsigned int j = 0; j < nb; j++) dw[pos[j]] = T(-1);		// clear weights
		pos[0] = n2c[i]; dw[pos[0]] = T(0); nb = 1;				// init first node of this community
		for allNodes(j) { if (A(i,j) == T(0) || i == j) continue;		// for all links j of node i
			const unsigned int c = n2c[j]; const T w = A(i,j);
			if (dw[c] == T(-1)) { dw[c] = w; pos[nb++] = c; }		// add new community to neighborhood
			else dw[c] += w; };						// else update weight
		const unsigned int co = n2c[i]; remove(i,co,dw[co]);			// remove node n from community co
		const T w = deg(i); unsigned int cb = co; T gb{0};	  		// find the nearest community bc
		for (unsigned int j = 0; j < nb; j++) { const unsigned int c = pos[j];	// for all communities in the neighborhood
			const T g = dw[c]-tot[c]*w/gw; if (g > gb) { cb = c; gb = g; } }; // determine gain for adding a node to community c
		insert(i,cb,dw[cb]); if (cb != co) nc++; }				// insert node in community cb
	return nc;
}

template <typename T>
matD<T> lpropCommunities<T>::partition(uvecD& map) const
//! generates a new matrix of communities as computed by aggregate().
{	uvecD rn(N); rn = ~0u; for allNodes(i) rn[n2c(i)]++;				// collect communities
	unsigned int C = 0; for allNodes(i) if (rn[i] < N) rn[i] = C++;			// determine number of communities C
	std::vector<uvec> cm(C); for allNodes(i) cm[rn[n2c(i)]].push_back(i);		// put nodes into new communities
	matD<T> At(C,C); At = T(0); uvecD l2c(N); l2c = ~0u;				// allocate new matrix At
	for (unsigned int c = 0; c < C; c++) { 						// for all communities
		for (const auto i : cm[c]) { At(c,c) += A(i,i); l2c[i] = c;		// for all nodes i of community c
			for allNodes(j) { if (A(i,j) == T(0) || i == j) continue;	// for all links j of node i
				At(c,rn[n2c(j)]) += A(i,j); } } };			// update weight in community c
	for (unsigned int i = 0; i < map.N; i++)
		map[i] = map(i) < map.N? l2c[map(i)]: ~0u; 				// remap labels
	return At;
}

template <typename T>
bool lpropCommunities<T>::aggregate(const T eps)
//! computes communities of the graph for one level.
{	uvecD rn(N); for allNodes(i) rn[i] = i;						// permute node labels
	uniformDist<float> ud; for (unsigned int i = 0; i < N-1; i++)
		std::swap(rn[i],rn[i+std::floor(ud()*(N-i-1))]);
	bool imp = false; for (unsigned int it = 0; it < maxit; it++) { 		// while improvements made...
		const T q0 = modularity(); if (moveNodes(rn)) imp = true; 		// move nodes between communities	
		const T q1 = modularity(); if (q1-q0 < eps) break; };			// check convergence
	return imp;									// returns true if improvements were made
}

template <typename T>
void lpropCommunities<T>::detMove(const uvecD& map)
//! moves nodes between communities according to map.
{	vecD<T> dw(N); dw = T(-1); uvecD pos(N); pos = 0; unsigned int nb = N;
	for (unsigned int i = 0; i < map.N; i++) {					// for all nodes i
		if (n2c[i] == map(i)) continue;
		for (unsigned int j = 0; j < nb; j++) dw[pos[j]] = T(-1);		// clear weights
		pos[0] = n2c[i]; dw[pos[0]] = 0; nb = 1;				// init first node of this community
		for allNodes(j) { if (A(i,j) == T(0) || i == j) continue;		// for all links j of node i
			const unsigned int c = n2c[j]; const T w = A(i,j);
			if (dw[c] == T(-1)) { dw[c] = w; pos[nb++] = c; }		// add new community to neighborhood
			else dw[c] += w; };						// else update weight
		const unsigned int co = n2c[i]; remove(i,co,dw[co]);			// remove node n from community co
		const unsigned int cn = map(i);	insert(i,cn,dw[cn]); };			// add to new community cn
}


/*! \brief finds communities in a single graph using label propagation.
 */

/*!
\param[in] A is a matrix representing a graph
\param[in] eps a small number used for convergence
\return the mapping of nodes in g to communities
 */
//! \relates lpropCommunities
template <typename T>
uvecD lpropComm(const matD<T>& _A, const unsigned int nit = 10, const bool verbose = false)
//! finds communities in matrix A; returns mapping of nodes onto communities.
{	uvecD map; T qbest{0}; unsigned int nc = _A.M;
	for (unsigned int it = 0; it < nit; it++) {					// for nit tries
		matD<T> A = _A; uvecD m(A.M);						// copy matrix
		for (unsigned int i = 0; i < m.N; i++) m[i] = i;			// initialize mapping
		for (bool imp = true; imp; ) {						// while improvements made...
			lpropCommunities<T> lv(A); imp = lv.aggregate();		// aggregate nodes
			if (imp) A = lv.partition(m);					// partition into new graph
			else {	const T q = lv.modularity();				// check modularity
				if (q > qbest) { qbest = q; map = m; nc = A.M;		// keep best mapping so far
					if (verbose) { printf("[%4d] q %6.4f for %u communities.\r",
						it, qbest, nc); fflush(stdout); } } } } };
	return map;									// return best mapping
}

//! Implements an algorithm for finding communities in a graph using a statistical approach.

/*
 * refer to:
 * Wilson J.D., Wang S., Mucha P.J., Bhamidi S., Nobel A.B. (2014)
 * A testing based extraction algorithm for identifying significant communities in networks
 * Ann. Appl. Stat. 8, 1853-1891 .
 *
 */
template<typename T> class esscCommunities {
	const matD<T>&	A;
	const unsigned int N;
	const vecD<T>	wt;
	const T		tw;
	
	uvec	search(uvec& B0, const T alpha, const char* dist) const			//! searches communities row B0.
		{ vecD<T> duBs(N), pvals(N), pbh(N); uvec B1(B0.size(),0);
		  for (unsigned int it = 0; it < 30; it++) {				// for all iterations
			if (std::equal(B0.begin(),B0.end(),B1.begin())) break;		// converged
			if (it > 0) B0 = B1;
			for allNodes(n) { T s{0}; 					// get probabilities for B0
				for (const auto b : B0) { s += A(n,b); }; duBs[n] = s; };
			T pB{0}; for (const auto b : B0) pB += wt(b); pB /= tw;
			if (strcmp(dist, "Binomial") == 0) {				// either weigh against binomial distribution
				for allNodes(n)	pvals[n] = wt(n) <= T(0)? T(1):
					T(1)-binomialDist<T>(wt(n),pB).cdf(duBs(n)); }
			else if (strcmp(dist, "Poisson") == 0) {			// else weigh against Poisson distribution
				for allNodes(n) pvals[n] = wt(n) <= T(0)? T(1):
					T(1)-PoissonDist<T>(wt(n)*pB).cdf(duBs(n)); }
			const uvecD rk = rank(pvals);					// divide probabilities by their rank
			for allNodes(n) pbh[rk(n)] = (pvals(rk(n))*N)/T(n+1);
			unsigned int c = 0; T thr{0}; B1.clear();			// find threshold
			for allNodes(n) { if (pbh(n) > alpha) continue;
				thr = std::max(thr,pvals(n)); c++; };
			if (c == 0) break;
			for allNodes(n) { if (pvals[n] <= thr) B1.push_back(n); };	// save significant members
			if (B1.size() == 0) break; };
		  return B1; }								// return set of significant members
	uvecD	order() const
		{ uvecD rn(N); for allNodes(n) rn[n] = n;
		  for (unsigned int n = 0; n < N-1; n++) std::swap(rn[n],rn[n+rand()%(N-n-1)]);
		  return rn; }

public:
	esscCommunities(const matD<T>& _A)
		: A(_A), N(A.M), wt(A.rowsum()), tw(wt.sum()) { }
	uvecD	work(const T alpha = T(0.01), const char* dist = "Binomial")		//! identifies statistically significant communities in undirected networks.
		{ uvecD rn = order(); std::vector<uvec> comm(N);
		  for allNodes(n) { const unsigned int r = rn(n);			// for each row r
			uvec B0; B0.push_back(r);					// get non-zero entries of A
			for allNodes(c) { if (A(r,c) > T(0)) B0.push_back(c); };
			const auto B1 = search(B0, alpha, dist);			// search for communities
			bool unq = true; if (B1.size() == 0) continue;
                	for (unsigned int i = 0; i < n; i++) {				// check if we found this one already
				if (std::equal(comm[i].begin(),comm[i].end(),B1.begin())) {
					unq = false; break; } };
			if (unq) comm[n] = B1; };					// nope, so save it
		  uvecD map(N); map = ~0u; unsigned int k = 0;				// collect compact map of communities
		  for allNodes(n) { if (comm[n].size() == 0) continue;
			for (const auto& c : comm[n]) { map[c] = k; }; k++; };
		  return map; }								// NB: -1 denote background nodes
};

/*! \brief finds communities in a graph using a statistical approach.
 */

/*!
\param[in] A is a matrix representing a graph
\param[in] alpha a community significance threshold
\param[in] dist a null distribution ("Binomial" or "Poisson")
\return the mapping of nodes in g to communities
 */
//! \relates esscCommunities
template <typename T>
uvecD esscComm(const matD<T>& A, const T alpha, const char* dist)
//! finds communities in matrix A; returns mapping of nodes onto communities.
{	return esscCommunities<T>(A).work(alpha,dist);					// returns best mapping
}


/*! \brief finds optimal number of communities in a graph using a linear programming.
 */

/*
 * refer to:
 * Brandes U., Delling D., Gaertler M., Gorke R., Hoefer M., Nikoloski Z., Wagner D. (2008)
 * On Modularity Clustering, IEEE Transactions on Knowledge and Data Engineering 20, 172-188.
 *
 */
 
/*!
\param[in] A is a matrix representing a graph
\param[in] verbose switches logging on
\return the mapping of nodes in g to communities
 */
template <typename T>
uvecD optComm(const matD<T>& A, const bool verbose = false)
//! finds optimal communities in A. NB: uses double for interfacing to the glpk library.
{	const double coef[] = { 0.0,1.0,1.0,-2.0 }; int id[] = { 0,0,0,0 };		// transitivity constraint
	const unsigned int N = A.M, nv = N*(N+1)/2; 
	const vecD<T> wc = A.colsum(), wr = A.rowsum(); const T wt = A.sum();
#define OFF(a,b) ((b)*((b)+1)/2+(st+a))
	glp_prob* ip = glp_create_prob(); glp_set_obj_dir(ip,GLP_MAX);			// allocate GLPK object
	const unsigned int st = glp_add_cols(ip,nv);					// set number of variables
	for (unsigned int i = 0; i < nv; i++) glp_set_col_kind(ip,st+i,GLP_BV);		// set variable type to binary
	for allNodes(i) { glp_set_col_bnds(ip,OFF(i,i),GLP_FX,1.0,1.0);			// reflexivity
		for (unsigned int j = i+1; j < N; j++) {				// transitivity
			for (unsigned int k = j+1; k < N; k++) {
				const unsigned int r = glp_add_rows(ip,3);
				glp_set_row_bnds(ip,r+0,GLP_UP,0.0,1.0);		// set constraint i,j,k
				id[1] = OFF(i,j); id[2] = OFF(j,k); id[3] = OFF(i,k);
				glp_set_mat_row(ip,r+0,3,id,coef);
				glp_set_row_bnds(ip,r+1,GLP_UP,0.0,1.0);		// set constraint j,i,k
				id[1] = OFF(i,j); id[2] = OFF(i,k); id[3] = OFF(j,k);
				glp_set_mat_row(ip,r+1,3,id,coef);
				glp_set_row_bnds(ip,r+2,GLP_UP,0.0,1.0);		// set constraint i,k,j
				id[1] = OFF(i,k); id[2] = OFF(j,k); id[3] = OFF(i,j);
				glp_set_mat_row(ip,r+2,3,id,coef); } } };
	for allNodes(i) { double w = double(A(i,i)-wc(i)*wr(i)/wt);
		glp_set_obj_coef(ip,OFF(i,i),w);
		for (unsigned int j = i+1; j < N; j++) {				// objective function
			w = double(A(i,j)-wc(i)*wr(j)/wt+A(j,i)-wc(j)*wr(i)/wt);
			glp_set_obj_coef(ip,OFF(i,j),w); } };
	glp_term_out(verbose? GLP_ON: GLP_OFF); glp_iocp pm; glp_init_iocp(&pm);	// solve it
	pm.br_tech = GLP_BR_DTH; pm.bt_tech = GLP_BT_BLB; pm.presolve = pm.binarize = GLP_ON;
	glp_intopt(ip, &pm);
	uvecD map(N); map = 0; unsigned int c = 0;					// id of the last community
	for allNodes(i) { unsigned int j = 0;						// retrieve communities
		for (j = 0; j < i; j++) {
			if (glp_mip_col_val(ip,OFF(j,i)) == 1) { map[i] = map[j]; break; } };
		if (j == i) map[i] = c++; };
#undef OFF
	if (verbose) printf("q %5.3f for %u communities.\n", glp_mip_obj_val(ip)/wt, c);
	glp_delete_prob(ip); return map;
}

//! Implements an algorithm for finding communities in a graph using the infomap algorithm.

/*
 * refer to:
 * Rosvall M., Bergstrom C.T. (2010)
 * Mapping change in large networks, PLoS ONE 5, e8694
 *
 * based on code by:
 * Copyright (c) 2011 Martin Rosvall
 *
 */

//! implements a node for graphs in the infomap community detection algorithm.

template <typename T> struct infoNode {
	using infoLink = std::pair<unsigned int,T>;
	uvec	 mb;				//!< list of members (for communities)
	std::vector<infoLink> in;		//!< list of incoming links and weights
	std::vector<infoLink> out;		//!< list of outgoing links and weights
	T	self;				//!< weight of this node
	T	wt;				//!< total weight of this node (degree)
	T	niso;				//!< indicates node without links (or number of isolated nodes). 
	T	exit;				//!< exit flow
	T	size;				//!< node capacity

	infoNode() : mb(), in(), out(), self(T(0)), wt(T(0)), niso(T(0)), exit(T(0)),
		size(T(0)) { }								//! allocates an empty node for infomap.
	infoNode(const unsigned int m, const T w)
		: mb(), in(), out(), self(T(0)), wt(w), niso(T(0)), exit(T(0)),
		size(T(0)) { mb.push_back(m); }						//! allocates a node with label m and weight w for infomap.
	infoNode(const infoNode<T>& b)							//! copy-constructs a node for infomap.
		: mb(b.mb), in(b.in), out(b.out), self(b.self), wt(b.wt),
		niso(b.niso), exit(b.exit), size(b.size) { }
	const infoNode<T>& operator=(const infoNode<T>& b)				//! assigns a node for infomap.
		{ mb = b.mb; in = b.in; out = b.out; self = b.self; wt = b.wt;
		  niso = b.niso; exit = b.exit; size = b.size; return *this; }
};

//! implements a graph in the infomap community detection algorithm.

template <typename T> struct infoGraph { 
	T	alpha = T(0.15);		//!< teleportation probability
	T	beta = T(1)-alpha;
	unsigned int N;				//!< number of nodes (communities) in graph
	std::vector<infoNode<T>> node;		//!< the list of nodes (communities)
	unsigned int niso;			//!< number of isolated nodes. 
	uvec	iso;				//!< the list of isolated nodes.
	T	exit; 				//!< plogp of exit flow.
	T	exitFlow;			//!< total exit flow.
	T	logExit;			//!< sum of plogp of exit flows.
	T	logSize;			//!< sum of plogp of capactity. 
	T	logNode;			//!< plogp of total capacity.

	void	init();
	void	eigenvector();

	infoGraph(const matD<T>& A)
		: N(A.M), node(N), niso(0), iso(), exit(T(0)), exitFlow(T(0)), logExit(T(0)),
		logSize(T(0)), logNode(T(0))						//!< allocates a graph from dense matrix A.
		{ for (unsigned int n = 0; n < N; n++) node[n] = infoNode<T>(n,T(1));
		  for (unsigned int i = 0; i < A.M; i++) { 
			for (unsigned int j = 0; j < A.N; j++) {
				const T w = A(i,j); if (w == T(0) || j < i) continue;
				if (i == j) node[i].wt = w;
				else {	node[i].out.push_back({j,w});
					node[j].in.push_back({i,w}); } } };
		  init(); }
	infoGraph(const infoGraph& g)							//! copy-constructs a graph for infomap.
		: N(g.N), node(g.node), niso(g.niso), iso(g.iso), exit(g.exit),
		exitFlow(g.exitFlow), logExit(g.logExit), logSize(g.logSize),
		logNode(g.logNode) { }
	infoGraph(const infoGraph& g, const uvec& mb);		//!< allocates a graph from list of nodes mb.
	void	calibrate();
	T	partition();
};

template <class T>
inline T plogp(const T x)
{ return x > T(0)? (x)*std::log(x): T(0); }

template <class T>
infoGraph<T>::infoGraph(const infoGraph& g, const uvec& mb)
	: N(mb.size()), node(N), niso(0), iso(), exit(0), exitFlow(0), logExit(0), logSize(0), logNode(0)
//! construct a graph by extracting a subgraph from the given g.
{	for (unsigned int n = 0; n < N; n++) node[n] = infoNode<T>(n,T(1));
	std::set<int> smb; for (auto m : mb) smb.insert(m);
	uvecD ren(g.N); ren = ~0u; unsigned int n = 0;	
	for (const auto m : smb) { const infoNode<T>& fn = g.node[m]; infoNode<T>& nd = node[n];
		nd.wt = fn.wt; nd.self = fn.self; ren[m] = n;
		for (const auto& l : fn.out) {
			int t = l.first, r = ren[t]; if (t >= m) continue; 
			if (smb.find(t) != smb.end()) {
				nd.out.push_back({r,l.second});
				node[r].in.push_back({n,l.second}); } };
		for (const auto& l : fn.in) {
			int t = l.first, r = ren[t]; if (t >= m) continue;
			if (smb.find(t) != smb.end()) {
				nd.in.push_back({r,l.second});
				node[r].out.push_back({n,l.second}); } };
		n++; };
	init();
}

template <class T>
void infoGraph<T>::init()
//! computes the flow inside the graph, counts iso nodes, normalizes edge weights
//! computes steady state distribution, and computes codelength.
{	niso = 0; T tw = T(0); logNode = T(0);
	for (const auto& nd : node) tw += nd.wt;
	for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n]; nd.wt /= tw;	// normalize teleportation weight
		if (nd.out.empty() && nd.self <= T(0)) { iso.push_back(n); niso++; }
		else {	T s = nd.self; for (auto& l : nd.out) s += l.second;
			nd.self /= s; for (auto& l : nd.out) l.second /= s; } }
	eigenvector();									// calculate steady state matrix
	for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n];		// update links to represent flow
		nd.self *= beta*nd.size;
		for (auto& l : nd.out) l.second *= beta*nd.size;
		for (auto& l : nd.out) {
			for (auto& li : node[l.first].in) {
				if (li.first == n) { li.second = l.second; break; } } };
		nd.niso = (nd.out.empty() && nd.self <= T(0))? nd.size: 0;
		nd.exit = nd.size-(alpha*nd.size+beta*nd.niso)*nd.wt-nd.self;
		logNode += plogp(nd.size); }
	calibrate();
}

template <class T>
void infoGraph<T>::eigenvector()
//! computes steady state distribution over the network
{	vecD<T> tmp(N); tmp = T(1)/N;
	for (unsigned int it = 0; it < 200; it++) { T d = T(0), s = T(0), sd = T(0);	// calculate dangling size
		for (unsigned int n = 0; n < niso; n++) d += tmp[iso[n]];
		for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n];	// flow from weight	  
 			nd.size = (alpha+beta*d)*nd.wt; }
		for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n];	// flow from network steps	  
			nd.size += beta*nd.self*tmp[n];
			for (const auto& l : nd.out) node[l.first].size += beta*l.second*tmp[n]; };
		for (unsigned int n = 0; n < N; n++) s += node[n].size;			// normalize 
		for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n];
			nd.size /= s; sd += std::abs(nd.size-tmp[n]); tmp[n] = nd.size; }
		if (sd < T(1e-10) && it >= 50) break; }
}

template <class T>
void infoGraph<T>::calibrate()
//! computes the code length of the network.
{	logExit = T(0); exitFlow = T(0); logSize = T(0);
	for (unsigned int n = 0; n < N; n++) { infoNode<T>& nd = node[n];
		logSize	+= plogp(nd.exit+nd.size); exitFlow += nd.exit;
		logExit += plogp(nd.exit); }
	exit = plogp(exitFlow);
}

//! implements an optimization context for the infomap community detection algorithm.

template <class T> struct Greedy {
	infoGraph<T>& g;			//!< the graph.
	unsigned int N;				//!< number of nodes in the graph.
	unsigned int E;				//!< number of empty (isolated) nodes (or communities).
	uvecD	mempty;				//!< list of empty (isolated) nodes (or communities).
	uvecD	nodeId;				//!< module number of each node
	vecD<T> mexit;	 			//!< exit flow per node (or community).
	vecD<T> msize;	 			//!< capacity per node (or community).
	vecD<T> mniso;				//!< list of isolated nodes.
	vecD<T> mwt;				//!< list of nodes (or communities) weight.
	uvecD	mmb;				//!< size of communities.
	T	logNode;			//!< plogp of total capacity.
	T	logExit;			//!< sum of plogp of exit flows.
	T	logSize;			//!< sum of plogp of capactity. 
	T	exitFlow;			//!< total exit flow.
	T	exit; 				//!< plogp of exit flow.
public:
	Greedy(infoGraph<T>& _g)
		: g(_g), N(g.N), E(0), mempty(N), nodeId(N), mexit(N), msize(N),
		mniso(N), mwt(N), mmb(N), logNode(g.logNode), logExit(g.logExit),
		logSize(g.logSize), exitFlow(g.exitFlow), exit(plogp(exitFlow))
		{ for (unsigned int n = 0; n < N; n++) { const infoNode<T>& nd = g.node[n];
	 		nodeId[n] = n; mexit[n] = nd.exit; msize[n] = nd.size;
	 		mniso[n] = nd.niso; mwt[n] = nd.wt;
	 		mmb[n] = nd.mb.size(); } }
	void	setMove(const uvecD& cls);
	bool	optimize();
	void	apply(const bool sort);
	T	codelen() const								//! returns MDL for this clustering.
	 	{ return exit-2.0*logExit+logSize-logNode; }
};

template <class T>
bool Greedy<T>::optimize()
//! moves node that yields best MDL compression.
{	const T cl = codelen(); bool moved = false; uniformDist<T> ud;
	unsigned int offset = 1; uvecD ro(N), red(N); red = 0;
	for (unsigned int n = 0; n < N; n++) ro[n] = n;					// generate random enumeration of nodes
	for (unsigned int n = 0; n < N-1; n++) std::swap(ro[n],ro[n+ud()*(N-n-1)]);
	for (unsigned int n = 0; n < N; n++) {						// pick nodes in random order
		const unsigned int f = ro[n], m = nodeId[f]; infoNode<T>& nd = g.node[f];
	 	if (offset > INT_MAX) { offset = 1; red = 0; }
		std::vector<std::pair<unsigned int,std::pair<T,T>>> flow;
	 	if (nd.out.size() == 0) {						// dangling node, add node to calculate flow below
	 	 	 red[m] = offset; flow.push_back({m,{T(0),T(0)}}); }
	 	else {	for (const auto& l : nd.out) {					// for all out links
				unsigned int nb = nodeId[l.first];
				if (red[nb] >= offset)
					flow[red[nb]-offset].second.first += l.second;
				else {	red[nb] = offset+flow.size();
		 			flow.push_back({nb,{l.second,T(0)}}); } } };
	 	for (const auto& l : nd.in) {	 	 				// for all in links
	 	 	 int nb = nodeId[l.first];
	 	 	 if (red[nb] >= offset)
				flow[red[nb]-offset].second.second += l.second;
	 	 	 else {	red[nb] = offset+flow.size();
		 		flow.push_back({nb,{l.second,T(0)}}); } };
		for (auto& l : flow) {							// for isolated nodes
	 	 	if (l.first == m) {
				l.second.first += (g.alpha*nd.size+g.beta*nd.niso)*(mwt[m]-nd.wt);
				l.second.second += (g.alpha*(msize[m]-nd.size)+g.beta*(mniso[m]-nd.niso))*nd.wt; }
			else {	l.second.first += (g.alpha*nd.size+g.beta*nd.niso)*mwt[l.first];
				l.second.second += (g.alpha*msize[l.first]+g.beta*mniso[l.first])*nd.wt; } };
		T mout = (g.alpha*nd.size+g.beta*nd.niso)*(mwt[m]-nd.wt);		// calculate flow to/from own module
		T min = (g.alpha*(msize[m]-nd.size)+g.beta*(mniso[m]-nd.niso))*nd.wt;
		if (red[m] >= offset) {
	 		mout = flow[red[m]-offset].second.first;
	 		min = flow[red[m]-offset].second.second; }
		if (mmb[m] > nd.mb.size() && E > 0)					// move to empty module (if node not already alone)
			flow.push_back({mempty[E-1],{T(0),T(0)}});
		const unsigned int nl = flow.size();
	 	for (unsigned int j = 0; j < nl-1; j++) {				// randomize link order for optimized search
	 	 	std::swap(flow[j],flow[j+ud()*(nl-j-1)]); }
	 	unsigned int b = m; T bout = T(0), bin = T(0), bd = T(0);
	 	for (const auto& l : flow) { if (l.first == m) continue; 	 	// find the move that minimizes the description length
	 		const T fout = l.second.first, fin = l.second.second;
			T delta_exit = plogp(exitFlow+mout+min-fout-fin)-exit;
			T delta_logExit = -plogp(mexit[m])-plogp(mexit[l.first]) 
				+plogp(mexit[m]-nd.exit+mout+min)+plogp(mexit[l.first]+nd.exit-fout-fin);
			T delta_logSize = -plogp(mexit[m]+msize[m])
				-plogp(mexit[l.first]+msize[l.first])
		 		+plogp(mexit[m]+msize[m]-nd.exit-nd.size+mout+min)
		 		+plogp(mexit[l.first]+msize[l.first]+nd.exit+nd.size-fout-fin);
			T d = delta_exit-T(2)*delta_logExit + delta_logSize;
			if (d-bd < T(-1e-10)) { b = l.first; bout = fout; bin = fin; bd = d; } };
	 	 if (b != m) {	 	 						// make best possible move
			if (mmb[b] == 0) E--;						// update empty module vector
			if (mmb[m] == nd.mb.size()) { mempty[E] = m; E++; }
	 	 	exitFlow -= mexit[m]+mexit[b];
	 	 	logExit -= plogp(mexit[m])+plogp(mexit[b]);
	 	 	logSize -= plogp(mexit[m]+msize[m])+plogp(mexit[b]+msize[b]);
	 	 	mexit[m] -= nd.exit-mout-min; msize[m] -= nd.size; mniso[m] -= nd.niso;
	 	 	mwt[m] -= nd.wt; mmb[m] -= nd.mb.size();
	 	 	mexit[b] += nd.exit-bout-bin; msize[b] += nd.size; mniso[b] += nd.niso;
	 	 	mwt[b] += nd.wt; mmb[b] += nd.mb.size();
	 	 	exitFlow += mexit[m]+mexit[b];
	 	 	logExit += plogp(mexit[m])+plogp(mexit[b]);			// update terms in map equation
	 	 	logSize += plogp(mexit[m]+msize[m])+plogp(mexit[b]+msize[b]);
	 	 	exit = plogp(exitFlow); nodeId[f] = b; moved = true; }
	 	 offset += N; }
	 return moved && codelen() < cl;
}

template <class T>
void Greedy<T>::apply(const bool sort)
//! build communities for next level.
{	uvec mod; uvecD n2m(N);					// contains ids of no-empty modules (nodes)
	if (sort) { std::multimap<T,unsigned int> map;
	 	 for (unsigned int n = 0; n < N; n++) {
	 	 	if (mmb[n]) map.insert({msize[n],n}); };
	 	 for (auto it = map.rbegin(); it != map.rend(); it++)
			mod.push_back(it->second); }
	else {	for (unsigned int n = 0; n < N; n++) { 
	 	 	 if (mmb[n]) mod.push_back(n); } };
	const unsigned int M = mod.size(); std::vector<infoNode<T>> tmp(M);		// create a new node list
	for (unsigned int m = 0; m < M; m++) {						// create new nodes
		infoNode<T>& nd = tmp[m]; nd.mb.clear();
		const unsigned int n = mod[m]; n2m[n] = m;
	 	nd.exit = mexit[n]; nd.size = msize[n];
	 	nd.niso = mniso[n]; nd.wt = mwt[n]; }	
	std::vector<std::map<unsigned int,T>> fout(M);					// calculate fout of links to different modules
	for (unsigned int n = 0; n < N; n++) { const infoNode<T>& nd = g.node[n];
	 	 const unsigned int m = n2m[nodeId[n]];					// final id of the module of the node i
	 	 std::copy(nd.mb.begin(), nd.mb.end(), std::back_inserter(tmp[m].mb) );
	 	 for (const auto& l : nd.out) {
	 	 	const unsigned int nm = n2m[nodeId[l.first]];
	 	 	if (l.first != n) { auto it = fout[m].find(nm);
				if (it != fout[m].end()) it->second += l.second;
			else fout[m].insert({nm,l.second}); } } };
	for (unsigned int m = 0; m < M; m++) { infoNode<T>& nd = tmp[m];		// create out at new level
	 	 for (auto& l : fout[m]) {
	 	 	 if (l.first != m) nd.out.push_back(l); } };
	std::vector<std::map<unsigned int,T>> fin(M);					// calculate outflow to different modules
	for (unsigned int n = 0; n < N; n++) { const infoNode<T>& nd = g.node[n];
	 	const unsigned int m = n2m[nodeId[n]];
		for (const auto& l : nd.in) {
	 	 	const unsigned int nm = n2m[nodeId[l.first]];
	 	 	if (l.first != n) { auto it = fin[m].find(nm);
				if (it != fin[m].end()) it->second += l.second;
			else fin[m].insert({nm,l.second}); } } };
	for (unsigned int m = 0; m < M; m++) { infoNode<T>& nd = tmp[m];		// create inflow from different modules
	 	 for (auto& l : fin[m]) { if (l.first != m) nd.in.push_back(l); } };
	mempty = ~0u; E = 0; g.node = tmp; g.N = M; N = M; g.calibrate();		// clear empty module
}

template <class T>
void Greedy<T>::setMove(const uvecD& cls)
//! computes the new code size if modules are merged as indicated by cls.
{	 for (unsigned int n = 0; n < N; n++) { const infoNode<T>& nd = g.node[n];
		const unsigned int t = cls(n); if (t == n) continue;
	 	T foutN = (g.alpha*nd.size+g.beta*nd.niso)*(mwt[n]-nd.wt);
	 	T finN = (g.alpha*(msize[n]-nd.size)+g.beta*(mniso[n]-nd.niso))*nd.wt;
	 	T foutT = (g.alpha*nd.size+g.beta*nd.niso)*mwt[t];
	 	T finT = (g.alpha*msize[t]+g.beta*mniso[t])*nd.wt;
	 	for (const auto& l : nd.out) { unsigned int nb = nodeId[l.first];	// for all out
			if (nb == n) foutN += l.second;
			else if (nb == t) foutT += l.second; };
	 	for (const auto& l : nd.in) { unsigned int nb = nodeId[l.first];	// for all in
			if (nb == n) finN += l.second;
			else if (nb == t) finT += l.second; };
	 	if (mmb[t] == 0) E--;							// update empty module vector
	 	if (mmb[n] == nd.mb.size()) { mempty[E] = n; E++; }
	 	exitFlow -= mexit[n] + mexit[t];
	 	logExit -= plogp(mexit[n])+plogp(mexit[t]);
	 	logSize -= plogp(mexit[n]+msize[n])+plogp(mexit[t]+msize[t]);
	 	mexit[n] -= nd.exit-foutN-finN; msize[n] -= nd.size;
	 	mniso[n] -= nd.niso; mwt[n] -= nd.wt; mmb[n] -= nd.mb.size();
	 	mexit[t] += nd.exit-foutT-finT; msize[t] += nd.size;
	 	mniso[t] += nd.niso; mwt[t] += nd.wt; mmb[t] += nd.mb.size();
	 	exitFlow += mexit[n] + mexit[t];
	 	logExit += plogp(mexit[n]) + plogp(mexit[t]);
	 	logSize += plogp(mexit[n] + msize[n]) + plogp(mexit[t] + msize[t]);
	 	exit = plogp(exitFlow); nodeId[n] = t; };
}

template <class T>
T infoGraph<T>::partition()
//! computes communities for this graph; returns MDL.
{	infoGraph<T> tmp(*this); uvecD cls(N), scls(N);
	bool cinit = false, sinit = false; T ocl = T(1000);
	for (unsigned int it = 0; it < 1000; it++) {
		if (it > 0) { cinit = true;						// indicates current clustering
	 	 	if (it % 2 == 0 && N > 1) { sinit = true; unsigned int si = 0;
				for (unsigned int i = 0; i < N; i++) { const auto& mb = node[i].mb;
		 			if (mb.size() <= 1) { scls[mb[0]] = si; cls[si++] = i; continue; };
		 	 		infoGraph<T> sf(tmp,mb); sf.partition();	// else extract subgraph and partition
					for (const infoNode<T>& nd : sf.node) {		// record membership changes
						for (const auto m : nd.mb) scls[mb[m]] = si;
						cls[si++] = i; } } }
			else {	for (unsigned int n = 0; n < N; n++)			// for each module
					for (const auto m : node[n].mb) cls[m] = n; }
			*this = tmp;
			if (sinit) { Greedy<T> opt(*this); opt.setMove(scls);		// use scls clustering
				opt.apply(false); sinit = false; } }			// flag it's used
		T ncl = T(0);
		while (true) {	Greedy<T> opt(*this);
			if (cinit) { cinit = false; opt.setMove(cls); }			// use cls clustering and flag it's used
			const T icl = opt.codelen();	 	 	 
			while (opt.optimize()) ;					// main greedy optimizing loop
	 	 	opt.apply(true); ncl = opt.codelen();
			if (std::abs(ncl-icl) < T(1e-10)) break; };
	 	if (ocl-ncl >= T(1e-10)) ocl = ncl; else return ncl; };
	throw rtException("Partitioning did not converge");
}

/*! \brief finds communities in a graph using the infomap algorithm.
 */

/*!
\param[in] A is a matrix representing a graph
\param[in] nt number of permutations
\return the mapping of nodes in g to communities
 */
//! \relates esscCommunities
template <class T>
uvecD infoComm(const matD<T>& A, const unsigned int nt, const bool verbose = false)
{	infoGraph<T> g(A); uvecD mem(A.M); T best = T(1000);
	for (unsigned int t = 0; t < nt; t++) {
		infoGraph<T> tmp(g); const T cl = tmp.partition(); if (cl >= best) continue;
		best = cl;
		if (verbose) { printf("[%4d] code %f for %u communities.\r", t, best, tmp.N);
				fflush(stdout); };
		for (unsigned int n = 0; n < tmp.N; n++) {
			for (const auto m : tmp.node[n].mb) mem[m] = n; } };
	return mem;
}
#undef allNodes

#endif

