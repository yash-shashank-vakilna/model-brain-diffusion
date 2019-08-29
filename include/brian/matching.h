#ifndef MATCHING_H
#define MATCHING_H

/*
 *
 * matching.h: matching of uni- and bipartite graphs
 * BRIAN Software Package Version 3.0
 *
 * $Id: matching.h 455 2016-12-12 15:08:30Z frithjof $
 *
 * 0.10 (11/08/16): initial version FK
 * v406 (28/09/16): bumped to version 3.0
 *
 */

//! Provides a context for solving linear optimization using the Munkres algorithm.

/*
 * based on code by:
 * Copyright (c) 2007 John Weaver
 * Copyright (c) 2015 Miroslav Krajicek
 *
 */

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

//! Implements an algorithm for matching two undirected graphs via graduated assignment.

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


//! Implements an algorithm for finding a consistent matches in a set of graphs.

/*
 * refer to:
 * Yan J, Cho M, Zha H, Yang X (2016)
 * Multi-Graph Matching via Affinity Optimization with Graduated Consistency Regularization
 * IEEE TPAMI 38, 1228-1242
 *
 * Algorithm 2 CAO-C on p.1233
 *
 */

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

