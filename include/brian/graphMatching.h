#ifndef GRAPHMATCHING_H
#define GRAPHMATCHING_H

/*
 *
 * graphMatching.h: matching of uni- and bipartite graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: graphMatching.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (11/08/16): initial version FK
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief matching of uni- and bipartite graphs.

    \details For detailed information, refer to (\cite Gold96, <a href="ref59.pdf" target="_blank">PDF</a>;
		\cite Beckett16, <a href="ref60.pdf" target="_blank">PDF</a>).
*/

/*
 * based on code by:
 * Copyright (c) 2007 John Weaver
 * Copyright (c) 2015 Miroslav Krajicek
 *
 */

//! Provides a context for solving linear optimization using the Munkres algorithm.

template<typename T> class Munkres {
	enum class state { normal = 0, star = 1, prime = 2 };
	matD<T> m;				//!< distance matrix.
	matD<state> mask;			//!< cell state matrix.
	bvecD	rmask;				//!< row mask vector.
	bvecD	cmask;				//!< column mask vector.
	unsigned int saver;			//!< saved row index.
	unsigned int savec;			//!< saved column index.
	using rcpair = std::pair<unsigned int,unsigned int>;

	bool	findZero(unsigned int& r, unsigned int& c)				//! returns row and column of next zero element.
		{ for (r = 0 ; r < m.M; r++) { if (rmask[r]) continue;
			for (c = 0; c < m.N; c++) { if (cmask[c]) continue;
				if (isSmall<T>(m(r,c))) return true; } };
		  return false; }
	bool	inList(const rcpair& n, const std::list<rcpair>& l) const		//! returns true if pair n is in list l.
		{ for (const auto i : l) if (n == i) return true;
		  return false; }
	unsigned int step1()								//! stars the first element in each row/column that contains a zero.
		{ for (unsigned int r = 0 ; r < m.M; r++) {
			for (unsigned int c = 0; c < m.N; c++) {
				if (isSmall<T>(m(r,c))) {
		  			for (unsigned int i = 0; i < r; i++)
						if (state::star == mask(i,c)) goto next_col;
		  			mask(r,c) = state::star; goto next_row; }
			next_col: ; }
		  next_row: ; }
		  return 2; }
	unsigned int step2()								//! masks columns that contain a starred element.
		{ unsigned int n = 0;
		  for (unsigned int r = 0; r < m.M; r++)
	  		for (unsigned int c = 0; c < m.N; c++)
				if (state::star == mask(r,c)) { cmask[c] = true; n++; }
		  return n >= m.M? 0: 3; }						// if all are masked, we're done.
	unsigned int step3()								//! implements step 3 of the Munkres algorithm.
		{ if (findZero(saver,savec)) mask(saver,savec) = state::prime;		// find an uncovered zero in the matrix and prime it.
		  else return 5;							// if no such zero exists, go to Step 5
		  for (unsigned int c = 0; c < m.N; c++) {
	  		if (mask(saver,c) != state::star) continue;			// cover this row and uncover the column containing the starred zero
			rmask[saver] = true; cmask[c] = false; return 3; };		// return to Step 3.1 to find a new Z
		  return 4; }								// if No Z* exists in the r of the Z', go to Step 4.
	void	replaceState(const rcpair& p, const state is, const state as)		//! if pair p has state is, masks pair with state as.
		{ if (mask(p.first,p.second) == is) mask(p.first,p.second) = as; }
	unsigned int step4()								//! implements step 4 of the Munkres algorithm.
		{ std::list<rcpair> seq; seq.push_back({saver, savec});			// seq contains pairs of row/col values where we have found a star or a prime
		  unsigned int r,c = savec; bool madepair;
		  do { madepair = false;						// construct the alternating sequence of primed and starred zeros
	  		for (r = 0; r < m.M; r++) { if (mask(r,c) != state::star) continue;
				if (inList({r,c},seq)) continue;
		  		madepair = true; seq.push_back({r,c}); break; };
	  		if (!madepair) break;
			madepair = false;
			for (c = 0; c < m.N; c++) { if (mask(r,c) != state::prime) continue;
		  		if (inList({r,c},seq)) continue;
		  		madepair = true; seq.push_back({r,c}); break; };
		  } while (madepair);
		  for (const auto p : seq) {
			replaceState(p,state::star,state::normal);			// unstar each starred zero of the sequence
			replaceState(p,state::prime,state::star); } 			// star each primed zero of the sequence
		  for (unsigned int r = 0; r < mask.M; r++)				// erase all primes, uncover all cols and rows
	  		for (unsigned int c = 0; c < mask.N; c++)
				replaceState({r,c},state::prime,state::normal);
		  rmask = false; cmask = false; return 2; }
	unsigned int step5()								//! implements step 5 of the Munkres algorithm.
		{ T h = std::numeric_limits<T>::max();
		  for (unsigned int r = 0; r < m.M; r++) { if (rmask[r]) continue;	// let h be the smallest uncovered entry in m
			for (unsigned int c = 0; c < m.N; c++) { if (cmask[c]) continue;
				if (h > m(r,c) && isSmall<T>(m(r,c)) == false) h = m(r,c); } };
		  for (unsigned int r = 0; r < m.M; r++) { if (rmask[r] == false) continue;
			for (unsigned int c = 0 ; c < m.N; c++) m(r,c) += h; }; 	// add h to all covered rows
		  for (unsigned int c = 0; c < m.N; c++) { if (cmask[c]) continue;
			for (unsigned int r = 0; r < m.M; r++) m(r,c) -= h; }; 		// subtract h from all uncovered cols
		  return 3; }								// return to step 3, without altering stars, primes, or covers.
	void	minimizeDir(const bool oc)						//! minimizes entries along direction dir.
		{ const unsigned int os = oc? m.N: m.M, is = oc? m.M: m.N;
	  	  for (unsigned int i = 0; i < os; i++) {	  			// look for a minimum value to subtract from all values along the "outer" direction.
			T min = oc? m(0,i): m(i,0);
			for (unsigned int j = 1; j < is && min > T(0); j++)		// if the current minimum is greater than zero, keep looking
		  		min = std::min(min, oc? m(j,i): m(i,j));
			if (min > T(0)) {
		  		for (unsigned int j = 0; j < is; j++)
					if (oc) m(j,i) -= min; else m(i,j) -= min; } } }
public:
	Munkres<T>(const matD<T>& _m)							//! Allocates a context for linear optimization in matrix m.
		: m(_m), mask(m.M,m.N), rmask(m.M), cmask(m.N), saver(0), savec(0)
		{ assert(m.M == m.N); rmask = false; cmask = false; }
	uvecD	solve()									//! solves the optimization problem; returns row vector of matches.
		{ minimizeDir(true); minimizeDir(false);				// prepare the matrix
		  for (unsigned int step = 1; step; ) {					// follow the steps
			switch (step) {
			case 1: step = step1(); break;					// step is always 2
			case 2: step = step2(); break;					// step is always either 0 or 3
			case 3: step = step3(); break;					// step in [3, 4, 5]
			case 4: step = step4(); break;					// step is always 2
			case 5: step = step5(); break; } }				// step is always 3
		  for (unsigned int r = 0; r < m.M; r++)				// mark results
			for (unsigned int c = 0; c < m.N; c++) {
				if (mask(r,c) == state::star) m(r,c) = T(0);
				else m(r,c) = T(-1); }
		  uvecD res(m.M);
		  for (unsigned int r = 0; r < m.M; r++)				// res[r] contains best matching column c
			for (unsigned int c = 0; c < m.N; c++)
				if (m(r,c) == T(0)) { res[r] = c; break; };
		  return res; }
};

//! Implements an algorithm for matching two undirected graphs.

template<typename T> class gradAssignment {
	const T betai = T(0.1);			//!< initial temperature
	const T betaf = T(10);			//!< final temperature
	const T betar = T(1.5);			//!< temperature increase
	const unsigned int maxit = 1000;	//!< max number of iterations per temperature level
	const T eps = T(1e-4);			//!< convergence limit
	const matS<T>& A;			//!< points to affinity matrix
	const unsigned int n1;			//!< number of nodes in first graph
	const unsigned int n2;			//!< number of nodes in second graph

	matD<T>	computeGradients(const matD<T>& M, const T beta) const			//! returns gradient of energy functional.
		{ matD<T> Mt = M; 
		  for (unsigned int r = 0; r < A.M; r++) { T q{0};
			const unsigned int i = r % M.M, j = r / M.M;
			for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)
				q += A.x[t]*M.x[A.ci[t]];
			Mt(i,j) = std::exp(std::min(q*beta,T(50.0))); };
		  return Mt; }
public:
	gradAssignment(const matS<T>& _A, const unsigned int _n1, const unsigned int _n2)
	//! allocates a context for graph matching using graduated assignment.
		: A(_A), n1(_n1), n2(_n2) { }
	uvecD	work(const bool verbose = false) const					//! returns a vector of label assignments.
		{ matD<T> M(n2,n1); M = 1.0/(n1*n2);
		  for (T beta = betai; beta < betaf; beta *= betar) {
			for (unsigned it = 0; it < maxit; it++) {
				const matD<T> Mt = bistNormalize(computeGradients(M,beta));
				const T e = absdiff(Mt,M); M = Mt;
				if (verbose) { printf("[%3d] %6.4e %6.4e\r", it, beta, e); fflush(stdout); }
				if (e < eps) break; };
			if (verbose) printf("\n"); };
		  return greedyMax(M); }
	matD<T>	assign() const
		{ const uvecD u = work(); assert(u.N == n1); fmatD M(n1,n2); M = T(0);
		  for (unsigned int i = 0; i < u.N; i++) M(i,u(i)) = T(1);
		  return M; }
};

/*! \brief finds matches in two graphs.
 */

/*!
\param[in] A points to the (node,edge) affinity matrix of both graphs
\param[in] n1 corresponds to the number of nodes in the first graph
\param[in] n2 corresponds to the number of nodes in the second graph
\return a vector of correspondencies of nodes in graph 1 to nodes in graph 2.
 */
//! \relates gradAssignment
template <typename T> 
uvecD graphMatch(const matS<T>& A, const unsigned int n1, const unsigned int n2)
{	return gradAssignment<T>(A,n1,n2).work(false);
}

//! Implements an algorithm for finding communities in a bipartite graph.

template<typename T> class bipCommunities {
	using uvec = std::vector<unsigned int>;

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
		{ uvecD u = unique(v); matD<T> L(u.N,v.N); L = 0.0f;
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
		{ T q = checkDiv(div[i],div[j]); if (q <= qbest) return false;
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
for (unsigned int i = 0; i < red.N; i++) assert(red[i] < red.N);
}

//! \relates bipCommunities
template <typename T> 
vecD<T> bipCommunities<T>::initBD() const
//! initializes blue label degrees.
{	vecD<T> bd(blue.N); bd = 0.0f; if (maxLabel(blue) == 0) return bd;
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
		if (qn > 0.0f && qn <= qbest) {						// so do not jump out right here
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
\param[in] A points to the affinity matrix of both graphs
\param[in] mini corresponds to the minimal number of modules to consider
\param[in] reps corresponds to the number of iterations per additional module
\return the best modularity value
 */
//! \relates bipCommunities
template <typename T> 
T bipComm(uvecD& red, uvecD& blue, const matD<T>& A, const unsigned int mini = 0,
	const unsigned int reps = 10)
{	const bool flip = A.M > A.N, verbose = true;
	bipCommunities<T> lp(flip? trp(A): A); T qbest = lp.work(red,blue,0,verbose);
	if (mini) { const uvecD u = unique(red); const unsigned int nm = u.N;
		for (unsigned int m = 1; m <= mini; m++) {
			for (unsigned int it = 0; it < reps; it++) {
				uvecD r = red, b = blue; T q = lp.work(r,b,m+nm,false);
				if (q > qbest) { red = r; blue = b; qbest = q; } };
			if (verbose) { printf("[+%3d] %6.4f\r", m, qbest); fflush(stdout); } } };
	if (flip) std::swap(red,blue);
	if (verbose) { printf("final q %6.4f for %u communities.\n", qbest, unique(red).N); };
	return qbest;
}


/*
 * refer to:
 * Campigotto R, Cespedes PC, Guillaume JL (2014)
 * A Generalized and Adaptive Method for Community Detection
 * arVix:1406.2518v1
 *
 */

#define allNodes(n)	(unsigned int n = 0; n < N; n++)

//! Implements an algorithm for finding communities in a graph.

template<typename T> class uniCommunities {
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
	uniCommunities(matD<T>& _A)
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
unsigned int uniCommunities<T>::moveNodes(const uvecD& rn)
//! move nodes between communities; returns number of moves.
{	vecD<T> dw(N); dw = T(-1); uvecD pos(N); pos = 0; unsigned int nc = 0, nb = N;
	for (unsigned int r = 0; r < rn.N; r++) { const unsigned int i = rn(r);		// for all nodes i
		for (unsigned int j = 0; j < nb; j++) dw[pos[j]] = T(-1);		// clear weights
		pos[0] = n2c[i]; dw[pos[0]] = 0; nb = 1;				// init first node of this community
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
matD<T> uniCommunities<T>::partition(uvecD& map) const
//! generates a new matrix of communities as computed by aggregate().
{	uvecD rn(N); rn = ~0u; for allNodes(i) rn[n2c(i)]++;				// collect communities
	unsigned int C = 0; for allNodes(i) if (rn[i] < N) rn[i] = C++;			// determine number of communities C
	std::vector<std::vector<unsigned int>> cm(C); 
	for allNodes(i) cm[rn[n2c(i)]].push_back(i);					// put nodes into new communities
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
bool uniCommunities<T>::aggregate(const T eps)
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
void uniCommunities<T>::detMove(const uvecD& map)
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

#undef allNodes

/*! \brief finds communities in a graphs.
 */

/*!
\param[in] A is a matrix representing a graph
\param[in] eps a small number used for convergence
\return the mapping of nodes in g to communities
 */
//! \relates uniCommunities
template <typename T>
uvecD uniComm(const matD<T>& _A, const unsigned int nit = 10)
//! finds communities in matrix A; returns mapping of nodes onto communities.
{	uvecD map; T qbest{0}; unsigned int nc = _A.M;
	for (unsigned int it = 0; it < nit; it++) {					// for nit tries
		matD<T> A = _A; uvecD m(A.M);						// copy matrix
		for (unsigned int i = 0; i < m.N; i++) m[i] = i;			// initialize mapping
		for (bool imp = true; imp; ) {						// while improvements made...
			uniCommunities<T> lv(A); imp = lv.aggregate();			// aggregate nodes
			if (imp) A = lv.partition(m);					// partition into new graph
			else {	const T q = lv.modularity();				// check modularity
				if (q > qbest) { qbest = q; map = m; nc = A.M;		// keep best mapping so far
					printf("[%4d] q %6.4f for %u communities.\r",
						it, qbest, nc); fflush(stdout); } } } };
	printf("\n"); return map;							// return best mapping
}

/*! \brief finds communities in a graphs.
 */

/*!
\param[in] A is a matrix representing a graph
\param[in] map a vector of nodes mapping to communities
\return the modularity of this mapping
 */
//! \relates uniCommunities
template <typename T>
T uniComm(const matD<T>& _A, const uvecD& map)
//! builds communities in matrix A according to mapping map; returns modularity.
{	matD<T> A = _A; uniCommunities<T> lv(A); lv.detMove(map); return lv.modularity();
}

/*
 * refer to:
 * Yan J, Cho M, Zha H, Yang X (2016)
 * Multi-Graph Matching via Affinity Optimization with Graduated Consistency Regularization
 * IEEE TPAMI 38, 1228-1242
 *
 * Algorithm 2 CAO-C on p.1233
 *
 */

//! Implements an algorithm for finding a consistent matches in a set of graphs.

#define allGraphs(c) (unsigned int c = 0; c < N; c++)
#define allNodes(c) (unsigned int c = 0; c < n; c++)

template<typename T> class CAOc {
	const T lInit = T(0.2);			//!< initial consistency weight
	const T lStep = T(1.1);			//!< weight increase
	const T lMax = T(1.0);			//!< maximum weight
	std::vector<matD<T>> X;			//!< vector of initial assignment matrices
	const std::vector<matS<T>>& K;		//!< vector of affinity matrices
	const unsigned int N;			//!< number of graphs
	const unsigned int n;			//!< number of nodes
	const unsigned int inCnt;		//!< inlier count
	const T den;				//!< normalization constant
	matD<T>	mask;				//!< inlier node consistency mask

	std::vector<size_t> sortIndex(const std::vector<T>& v) const			//! returns index of elements sorted in descending order
		{ std::vector<size_t> u(v.size());
		  for (size_t i = 0; i < u.size(); i++) u[i] = i;
		  std::sort(u.begin(),u.end(),
			[&v](size_t a, size_t b) { return v[a] > v[b]; });
		  return u; }
	T	rowDiff(const matD<T>& P, const matD<T>& Q, unsigned int i) const	//! returns sum of absolute differences in row i
		{ T d{0}; for allNodes(j) d += std::abs(P(i,j)-Q(i,j));
		  return d; }
	T	maskedDiff(const matD<T>& P, const matD<T>& Q, unsigned int r) const	//! returns row-wise masked absolute difference between P and Q
		{ T d{0}; unsigned int s = 0;
		  for allNodes(i) { if (mask(i,r) > T(0)) { d += rowDiff(P,Q,i); s++; } };
		  return s? d/T(s): d; }
	matD<T> consistencyMask() const							//! computes node-wise consistency of each node in each graph, 
		{ matD<T> mask(n,N); mask = T(0);					// see Def 4 on p.1230
		  for allGraphs(r) { std::vector<T> c(n); for allNodes(i) c[i] = T(0);
			for (unsigned int i = 0; i < N-1; i++) {
				for (unsigned int j = i+1; j < N; j++) {
					const matD<T> P = X[i*N+r]*X[j*N+i], &Q = X[j*N+r];
		  			for allNodes(a) c[a] += T(1)-T(0.5)*rowDiff(P,Q,a); } };
			for allNodes(i) c[i] /= T(0.5)*T(N*(N-1));			// normalize the summation value
			std::vector<size_t> u = sortIndex(c);				// sort indices by descending values
			for (unsigned int i = 0; i < inCnt; i++) mask(u[i],r) = T(1); };
		  return mask; }
	T	inlierScore(unsigned int inCnt) const					//! returns normalization factor for findPath()
		{ T dmax{0}; const range r(0,inCnt);					// see footnote on p.1233
		  for allGraphs(i) { for (unsigned int j = i+1; j < N; j++) {
			const matD<T>& J = X[j*N+i]; vecD<T> v = asVector(J(r,r));
		  	for (unsigned int i = 0; i < v.N; i++) v[i] = std::round(v(i));
			const T d = dot(v,K[j*N+i]*v); dmax = std::max(d,dmax); } };
		  assert(dmax > T(0)); return dmax; }
	T	spConsistency(const matD<T>& C, unsigned int i, unsigned int j) const	//! computes single pair consistency of mapping between graphs i and j.
		{ T e{0}; for allGraphs(r) { if (r == i || r == j) continue;		// see Def 1 on p.1230
			const matD<T> P = X[i+r*N]*X[r+j*N]; e += maskedDiff(C,P,r); };
		  e = T(1)-e/(T(2)*N); assert(e >= T(0) && e <= T(1)); return e; }
	T	spAffinity(const matD<T>& P, const matS<T>& Q, unsigned int r) const	//! computes single pair affinity of mapping between graphs i and j.
		{ vecD<T> p = asVector(P); unsigned int k = 0;
		  for allNodes(i) { for allNodes(j) p[k++] *= mask(j,r); };
		  T a = dot(p,Q*p)/den; assert(a >= T(0) && a <= T(1)); return a;  }
	matD<T>	pConsistency(const std::vector<matD<T>>& X) const			//! computes pair consistency of mapping X.
		{ matD<T> C(N,N); C = T(0); for allGraphs(i) { C(i,i) = T(1);		// see Def 2 on p.1230
			for (unsigned int j = i+1; j < N; j++) { T e{0};
				for allGraphs(c) { const matD<T> P = X[c*N+i]*X[j*N+c];
					e += maskedDiff(P,X[j*N+i],c); };
				C(i,j) = T(1)-e/(T(2)*N); C(j,i) = C(i,j); } };
		  return C; }
	matD<T>	findPath(unsigned int i, unsigned int j, T lambda, bool ex) const	//! finds matching between graphs i and j.
		{ vecD<T> cn(N), af(N); cn = T(0); af = T(0);
		  for allGraphs(c) { const matD<T> P = X[c*N+i]*X[j*N+c];  		// for each candidate solution c
			if (ex) cn[c] = spConsistency(P,i,j);				// compute exact consistency and affinity
			af[c] = spAffinity(P,K[j*N+i],c); };
		  const vecD<T> ft = ex? (T(1)-lambda)*af+lambda*cn: af;		// compute fitness (Eq 7)
		  const unsigned int c = ft.imax(); return X[c*N+i]*X[j*N+c]; };	// return solution with maximum fitness
	bool	spanPath(std::vector<unsigned int>& seq, bvecD& vis, const matD<T>& A,
			unsigned int s, unsigned int e)					//! returns path seq from s to e on the span tree A
		{ bool succ = false; const unsigned int l = seq.size();
		  for allGraphs(c) { if (vis(c)) continue;				// for each unvisited graph c
			if (A(s,c) > T(0)) { seq.push_back(c); vis(c) = true;		// if edge (i,c) exists, add to path
				if (c == e) return true;				// return true if end point found
				succ |= spanPath(seq,vis,A,c,e); } };			// else continue from e
		  if (succ == false && l > 1) seq.pop_back();
		  return succ; }
	std::vector<matD<T>> mstMatch(const std::vector<matD<T>>& X)			//! post-process solution
		{ const matD<T> C = pConsistency(X), A = minimumSpanningTree(C);	// build minimum spanning tree on consistency matrix
		  std::vector<matD<T>> M(N*N); for allGraphs(i) { M[i*N+i] = X[i*N+i];
			for (unsigned int j = i+1; j < N; j++) {			// for each graph pair
        			std::vector<unsigned int> seq; seq.push_back(i);
        			bvecD vis(N); vis = false; vis(i) = true;
        			spanPath(seq,vis,A,i,j); matD<T> P = X[0];		// compute the path on tree
				for (unsigned int s = 0; s < seq.size()-1; s++)		// follow path on tree
					P = P*X[seq[s+1]*N+seq[s]];
        			M[j*N+i] = P; M[i*N+j] = trp(P); } };			// set results
		  return M; }
public:
	CAOc(const std::vector<matD<T>>& _X, const std::vector<matS<T>>& _K, unsigned int _in)
		: X(_X), K(_K), N(std::sqrt(X.size())), n(X[0].M), inCnt(_in),
		den(inlierScore(inCnt)), mask(consistencyMask())
		{ assert(N*N == X.size()); assert(X.size() == K.size()); }		// check that X is square and matches K
	std::vector<matD<T>> work(const unsigned int maxit, const T eps = T(1e-6), const bool verbose = false)	
		{ T lambda = lInit;
		  for (unsigned int it = 0; it < maxit; it++) {				// see Algorithm 2, line 3
		  	const bool ex = it >= 2; T d = T(0);
		  	for allGraphs(i) {
			for (unsigned int j = i+1; j < N; j++) {			// for each graph pair (line 4)
				if (verbose) { printf("pair %02d-%02d\r", i,j); fflush(stdout); };
            			const matD<T> P = findPath(i,j,lambda,ex);		// evaluate the fitness function (line 5)
				d += absdiff(X[j*N+i],P); 				// compute update
				X[j*N+i] = P; X[i*N+j] = trp(P); } };			// save matching
			d *= T(2)/T(N*N*n);
			if (verbose) { printf("[%2d] lambda %6.4f delta %6.4f\n", it,lambda,d); fflush(stdout); };
			if (d < eps) break; 						// check convergence
			mask = consistencyMask();					// recompute mask
			if (ex) lambda = std::min(lMax,lStep*lambda); };		// update weight (line 7)
		  return mstMatch(X); }							// post-processing (line 9)
};

#undef allGraphs
#undef allNodes

/*! \brief finds consistent matches in a set of graphs.
 */

/*!
\param[in] X references a vector of (n,n) matrices representing all pairings of n graphs
\param[in] K references to a vector of (n,n) affinity matrices of between all pairings of n graphs
\param[in] n is the number of inlier nodes
\return a vector of (n,n) matrices representing final pairings of n graphs.
 */
//! \relates CAOc
template <typename T> 
std::vector<matD<T>> graphSetMatch(std::vector<matD<T>>& X, const std::vector<matS<T>>& K,
	const unsigned int n)
{	return CAOc<T>(X,K,n).work(100,1e-6,true);
}

#endif

