#ifndef PRECONDITIONER
#define PRECONDITIONER

/*
 *
 * preconditioner.h: preconditioners for sparse linear systems
 * BRIAN Software Package Version 3.0
 *
 * $Id: preconditioner.h 508 2017-03-26 20:13:21Z frithjof $
 *
 * 0.10 (17/06/11): preconditioners introduced
 * 0.20 (23/09/12): split from solver
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.50 (26/11/14): block Jacobi preconditioner added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements preconditioners for sparse linear systems.
*/

#define allRows(A, r)	(unsigned int r = 0; r < A.M; r++)
#define allCols(A, t)	(unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)

//! Abstract base class of all preconditioners for iterative sparse solvers.

template<typename T, typename U = T> class precS {					// base class for all preconditioners
public:
	precS()										//! constructs a preconditioner.
		{ clearFpException(); }							// note: clears fp errors
        virtual ~precS() { }
        virtual vecD<U> operator()(const vecD<U>& b) const				//! applies preconditioner to vector b.
		{ return b; }								// returns identity
};

//! Implements a Jacobi preconditioner.

template<typename T, typename U = T> class precJacobi : public precS<T,U> {
	vecD<U>	d;				//!< inverse of diagonal
public:
	precJacobi(const matS<U>& A)							//! constructs a Jacobi preconditioner.
		: precS<T,U>(), d(A.idia()) { checkFpException(); }
        vecD<U>	operator()(const vecD<U>& b) const					//! applies Jacobi preconditioner to vector b.
		{ return d*b; }
};

//! Implements a block Jacobi preconditioner.

template<typename T, typename U = T> class precBlockJacobi : public precS<T,U> {
	const unsigned int bs;			//!< block size
	std::vector<matD<U>> Ai;		//!< local inverses
	matD<U>	getBlock(const matS<U>& A, const unsigned int s, const unsigned int l)	//! returns inverse sub-matrix A(s:s+l,s:s+l).
		{ matD<U> S(l,l); S = U(0);						// allocate sub-matrix
		  for (unsigned int r = s; r < s+l; r++) { 				// fill sub-matrix
			for (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++) {
				const unsigned int c = A.ci[t];
				if (c >= s && c < s+l) S(r-s,c-s) = A.x[t]; } };
		  return inv(S); }							// return inverse
	void	init(const matS<U>& A)							//! computes block inverses of size bs.
		{ for (unsigned int i = 0; i < Ai.size(); i++) {
			const unsigned int s = i*bs, l = std::min((i+1)*bs,A.M)-s;
			assert(l > 0); Ai[i] = getBlock(A,s,l); } }
public:
	precBlockJacobi(const matS<U>& A, const unsigned int _bs)			//! constructs a block Jacobi preconditioner.
		: precS<T,U>(), bs(_bs), Ai((A.M-1)/bs+1)
		{ init(A); checkFpException(); }
        vecD<U>	operator()(const vecD<U>& b) const					//! applies block Jacobi preconditioner to vector b.
		{ vecD<U> x(b.N);
		  for (unsigned int i = 0; i < Ai.size(); i++) {
			const unsigned int s = i*bs, l = Ai[i].M; vecD<U> v(l); v = U(0);
		 	for (unsigned int r = 0; r < l; r++) v[r] = b(r+s);
			const vecD<U> xi = Ai[i]*v;
		 	for (unsigned int r = 0; r < l; r++) x[r+s] = xi(r); };
		  return x; };
};

//! Helper class for approximate inverse preconditioner: represents a mapped vector.

template<typename T, typename U = T>
class airow {
public:
	struct map_entry {
		U	v;
		unsigned int id;

		map_entry(const U _v, unsigned int _id) : v(_v), id(_id) { }
		map_entry() { }
	};
	struct heap_entry {
		U	v;
		typename std::map<unsigned int,map_entry>::iterator it;

		heap_entry(const U& v, const typename std::map<unsigned int,map_entry>::iterator& i)
			: v(v), it(i) { }
		heap_entry() { }
	};
	std::map<unsigned int,map_entry> m;	// sorted by index
	std::vector<heap_entry> h;		// sorted by decreasing absolute value
private:
	void	swap(const unsigned int i, const unsigned int j)
	  	{ std::swap(h[i],h[j]);
		  h[i].it->second.id = i; h[j].it->second.id = j; }
	void	downheap(unsigned int i)
		{ unsigned int c0 = (i+1)*2-1, c1 = (i+1)*2;
		  while (c0 < h.size() || c1 < h.size()) { unsigned int mc = c0;
			if (c1 < h.size() && std::abs(h[c1].v) < std::abs(h[c0].v)) mc = c1;
			if (std::abs(h[c0].v) < std::abs(h[i].v) ||
				(c1 < h.size() && std::abs(h[c1].v < std::abs(h[i].v)))) swap(i,mc);
			else break;
			i = mc; c0 = (i+1)*2-1; c1 = (i+1)*2; } }
	void	upheap(unsigned int i)
		{ unsigned int p = (i-1)/2; while (i != 0) {
			if (std::abs(h[i].v) >= std::abs(h[p].v)) return;
			swap(i,p); i = p; p = (i-1)/2; } }
	void	replace(const unsigned int i, const U v)
		{ const U mv = h.empty()? U(0): std::abs(h.begin()->v);
		  if (std::abs(v) < mv) return;
		  if (h.empty() == false) { const auto it = h.begin()->it;
		  	swap(0,h.size()-1); h.pop_back(); downheap(0); m.erase(it); };
		  insert(i,v); }
	void	add(const unsigned int i, const U v)
		{ const auto it = m.find(i); it->second.v += v;
		  const unsigned int j = it->second.id;
		  const U ov = h[j].v; h[j].v = it->second.v;
		  if (std::abs(h[j].v) < std::abs(ov)) upheap(j); else downheap(j); }
	bool	find(const unsigned int i) const
		{ return m.count(i) != 0; }
public:
	airow()
		: m(), h() { }
	unsigned int size() const
		{ return m.size(); }
	std::map<unsigned int,U> mv(const matS<U>& A) const
		{ std::map<unsigned int,U> b;
		  for (const auto& e : m) {
			const U v = e.second.v; const unsigned int r = e.first;
			for allCols(A,t) { auto it = b.find(A.ci[t]);
				const U p = A.x[t]*v;
				if (it == b.end()) b[A.ci[t]] = p;
				else it->second += p; } };
		  return b; }
	T	dot(const std::map<unsigned int,U>& b) const
		{ auto ai = m.begin(); auto bi = b.begin(); T s{0};
		  while (ai != m.end() && bi != b.end()) {
			const unsigned int a_ind = ai->first, b_ind = bi->first;
			if (a_ind == b_ind) { s += ai->second.v * bi->second; ++ai; ++bi; }
			else if (a_ind < b_ind) ai++;
			else bi++; }
		  return s; }
	void	drop(airow<T,U>& x, const U v, const T tol, const unsigned int nz) const
		{ for (auto& e : m) { const unsigned int i = e.first;
			const U t = v*e.second.v; if (std::abs(t) < tol) continue;
			if (x.find(i)) x.add(i,t);
			else if (nz == ~0u || x.size() < nz) x.insert(i,t);
			else x.replace(i,t); } }
	void	insert(const unsigned int i, const U v)
		{ map_entry me(v,~0u); const auto it = m.insert(std::make_pair(i,me)).first;
		  heap_entry he(v,it); h.push_back(he);
		  he.it->second.id = h.size()-1; upheap(h.size()-1); }
};

template<typename T, typename U = T>
matS<U> asMatrix(const std::vector<airow<T,U>>& src)
{	unsigned int nz = 0; for (const auto& s : src) nz += s.size();
	matS<U> A(src.size(),src.size(),nz); unsigned int r = 0, t = 0; A.ri[r] = 0; 
	for (const auto& s : src) {
		for (const auto& e : s.m) { A.ci[t] = e.first; A.x[t++] = e.second.v; }
		A.ri[++r] = t; }
	return A;
}

template <typename T, typename U = T> class precAinv : public precS<T,U> {       
	matS<U> W;				//!< stores approximate inverse matrix
	vecD<U> d;				//!< stores inverse diagonal
public:
	precAinv(const matS<U>& A, const T tol, const unsigned int nz = 5,
		const bool lin = false, const unsigned int lp = 0)
		: precS<T,U>(), W(), d(A.M)
		{ std::vector<airow<T,U>> w(A.M);
		  for allRows(A,r) w[r].insert(r,T(1)); 
		  for allRows(A,r) { const auto u = w[r].mv(A); const T p = w[r].dot(u);
			if (isSmall<U>(p)) throw rtException("precAinv: small diagonal value");
			d[r] = U(1.0/p);
			for (auto it = u.upper_bound(r); it != u.end(); it++) {
				unsigned int i = it->first, rc = nz;
				if (lin) { rc = lp+A.ri[i+1]-A.ri[i]; if (rc < 1) rc = 1; }
				w[r].drop(w[i],-it->second/p,tol,rc); } }
		  W = asMatrix(w); }
	vecD<U>	operator()(const vecD<U>& x) const
		{ const vecD<U> t = (W*x)*d; return W.tm(t); }
};

template <typename T, typename U = T> class precAinvNS : public precS<T,U> {       
	matS<U> W;				//!< stores approximate inverse matrix
	matS<U> Z;				//!< stores approximate inverse matrix
	vecD<U> d;				//!< stores inverse diagonal
public:
	precAinvNS(const matS<U>& A, const T tol, const unsigned int nz = 5,
		const bool lin = false, const unsigned int lp = 0)
		: precS<T,U>(), W(), Z(), d(A.M)
		{ std::vector<airow<T,U>> wt(A.M), z(A.M); const matS<T> At = trp(A);
		  for allRows(A,r) { wt[r].insert(r,T(1)); z[r].insert(r,T(1)); };
		  for allRows(A,r) {
			const auto u = wt[r].mv(At), l = z[r].mv(A); const T p = wt[r].dot(l);
			if (isSmall<U>(p)) throw rtException("precAinvNS: small diagonal value");
			d[r] = U(1.0/p);
			for (auto it = u.upper_bound(r); it != u.end(); it++) {
				unsigned int i = it->first, rc = nz;
				if (lin) { rc = lp+A.ri[i+1]-A.ri[i]; if (rc < 1) rc = 1; }
				z[r].drop(z[i],-it->second/p,tol,rc); }
			for (auto it = l.upper_bound(r); it != l.end(); it++) {
				unsigned int i = it->first, rc = nz;
				if (lin) { rc = lp+A.ri[i+1]-A.ri[i]; if (rc < 1) rc = 1; }
				wt[r].drop(wt[i],-it->second/p,tol,rc); } }
		  Z = asMatrix(z); W = asMatrix(wt); }
	vecD<U>	operator()(const vecD<U>& x) const
		{ const vecD<U> t = (Z*x)*d; return W.tm(t); }
};

enum class amgSmootherType { Jacobi, GaussSeidel, polynomial };

//! Stores a level for an algebraic multigrid preconditioner.

template <typename T, typename U> class level {

	class smoother {								//! abstract base class for AMG smoothers.
	public:
		smoother() { }
        	virtual ~smoother() { }
		virtual smoother* clone() const	= 0;
		virtual void pre(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const = 0;
		virtual void post(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	= 0;
	};

	//! Jacobi smoother for AMG preconditioner.

	class Jacobi : public smoother {
		const vecD<U> inv;		//!< stores diagonal
	public:
		Jacobi(const matS<U>& A, const T omega = T(1))				//! constructs a Jacobi smoother for sparse matrix A.
			: smoother(), inv(A.idia()*omega) { }
		Jacobi(const Jacobi& b)							//! copies Jacobi smoother from b.
			: smoother(), inv(b.inv) { }
		Jacobi& operator=(const Jacobi& b)					//! assigns Jacobi smoother from b.
			{ if (this != &b) inv = b.inv; return *this; }
 		Jacobi& operator=(Jacobi&& b)						//! move assigns Jacobi smoother from b.
			{ assert(this != &b); inv = b.inv; return *this; }
		Jacobi*	clone() const							//! clones a Jacobi instance.
			{ return new Jacobi(*this); }
		void	pre(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies pre-multiplication x = A * b.
			{ (void)A; x = inv*b; }
		void	post(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies post-multiplication x = A * b.
			{ const vecD<U> y = b-A*x; x += inv*y; }
	};

	//! Gauss-Seidel smoother for AMG preconditioner.

	class GaussSeidel : public smoother {
		const unsigned int nit;		//!< number of GS iterations

		void	sweep(vecD<U>& x, const matS<U>& A, const vecD<U>& b,
				const int s, const int e, const int i) const		//! applies a Gauss-Seidel sweep x = A * b.
			{ if (x.N == 0) { x.resize(A.M); x = U(0); };
			  for (int r = s; r != e; r += i) { U d{0}, rs = U{0};
				for allCols(A,t) { int c = A.ci[t]; U v = A.x[t];
					if (c == r) d = v; else rs += v*x[c]; };
				if (d != U(0)) { x[r] = (b(r)-rs)/d; } } }
		void	cycle(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies sweep cycles.
			{ for (unsigned int it = 0; it < nit; it++) {
				sweep(x,A,b,0,A.M,1); sweep(x,A,b,A.M-1,-1,-1); } }
	public:
		GaussSeidel(const unsigned int n = 0)					//! constructs an empty Gauss-Seidel smoother.
			: smoother(), nit(n) { }
		GaussSeidel(const GaussSeidel& b)					//! copies Gauss-Seidel smoother from b.
			: smoother(), nit(b.nit) { }
		GaussSeidel& operator=(const GaussSeidel& b)				//! assigns Gauss-Seidel smoother from b.
			{ if (this != &b) nit = b.nit; return *this; }
 		GaussSeidel& operator=(GaussSeidel&& b)					//! move assigns Gauss-Seidel smoother from b.
			{ assert(this != &b); nit = b.nit; return *this; }
		GaussSeidel* clone() const						//! clones a Gauss-Seidel instance.
			{ return new GaussSeidel(*this); }
		void	pre(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies pre-multiplication x = A * b.
			{ cycle(x,A,b); }
		void	post(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies post-multiplication x = A * b.
			{ cycle(x,A,b); }
	};

	//! polynomial smoother for AMG preconditioner.

	class polynomial : public smoother {
		const vecD<T> cf;		//!< stores coefficients of the polynomial 

		vecD<T>	setCoefs(const T rho, const T lo = T(1.0/30.0),
				const T up = T(1.1)) const				//! determines coefficients of the polynomial.
			{ const unsigned int n = 3; vecD<T> d(n), rt(n);		// Chebyshev roots for the interval [-1,1]
			  const T x0 = lo*rho, x1 = up*rho;
			  for (unsigned int i = 0; i < n; i++)
				rt[i] = T(std::cos(M_PI*(i+0.5)/n));
			  for (unsigned int i = 0; i < n; i++)
				rt[i] = T(0.5)*(x1-x0)*(T(1)+rt[i])+x0;			// Chebyshev roots for the interval [x0,x1]
			  d[0] = T(1); d[1] = -(rt[0]+rt[1]+rt[2]);
			  d[2] = rt[0]*rt[1]+rt[1]*rt[2]+rt[2]*rt[0]; 
			  const T f = rt[0]*rt[1]*rt[2]; return d/f; }
	public:
		polynomial(const matS<U>& A)						//! constructs a polynomial smoother for eigenvalue rho.
			: smoother(), cf(setCoefs(real(A.spectralRadius()))) { }
		polynomial(const polynomial& b)						//! copies polynomial smoother from b.
			: smoother(), cf(b.cf) { }
		polynomial& operator=(const polynomial& b)				//! assigns polynomial smoother from b.
			{ if (this != &b) cf = b.cf; return *this; }
 		polynomial& operator=(polynomial&& b)					//! move assigns polynomial smoother from b.
			{ assert(this != &b); cf = b.cf; return *this; }
		polynomial* clone() const						//! clones a polynomial instance.
			{ return new polynomial(*this); }
		void	pre(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies pre-multiplication x = A * b.
			{ x = b*U(cf(0));
			  for (unsigned int i = 1; i < cf.N; i++) x += A*x+b*U(cf(i)); }
		void	post(vecD<U>& x, const matS<U>& A, const vecD<U>& b) const	//! applies post-multiplication x = A * b.
			{ vecD<U> r = b-A*x, h = r*U(cf(0));
			  for (unsigned int i = 1; i < cf.N; i++) {
				h = A*h+r*U(cf(i)); x += h; } }
	};

       	const matS<U> A;			//!< LHS
       	vecD<U> b;				//!< RHS
       	matS<U> P;				//!< prolongation operator
       	matS<U> R;				//!< restriction operator
       	smoother* sm;				//!< smoothing operator

	smoother* initSmoother(const amgSmootherType t) const
		{ switch (t) {
		  case amgSmootherType::Jacobi:
			return new Jacobi(A);
		  case amgSmootherType::GaussSeidel:
			return new GaussSeidel(10);
		  case amgSmootherType::polynomial:
			return new polynomial(A);
		  default: throw optException("unknown smoother type"); } }
	vecD<U>	setOperators(const T th, const U omega);
public:
	level(const matS<U>& _A, const amgSmootherType t)				//! constructs a level for matrix A.
		: A(_A), b(_A.M), P(), R(), sm(initSmoother(t)) { b = U(1e-3); }
	level(const matS<U>& _A, const vecD<U>& _b, const amgSmootherType t)		//! construct a level for matrix A and RHS b.
		: A(_A), b(_b), P(), R(), sm(initSmoother(t)) { }
        level(const level& l)								//! copy constructs from level l.
		: A(l.A), b(l.b), P(l.P), R(l.R), sm(l.sm->clone()) { }
        level(level&& l)								//! move-constructs from level l.
		: A(l.A), b(l.b), P(l.P), R(l.R), sm(l.sm) { l.sm = nullptr; }
        level& operator=(const level& l)						//! assigns from level l.
		{ if (this != &l) { A = l.A; b = l.b; R = l.P, R = l.R;
			sm = l.sm->clone(); }; return *this; }
	~level() { delete sm; }
	unsigned int dim() const							//! returns number of rows.
		{ return A.m(); }
	vecD<U>	coarse(vecD<U>& xc, const vecD<U>& bc) const				//! reduce solution.
		{ sm->pre(xc,A,bc); return R*(bc-A*xc); }
	void	correct(vecD<U>& x, const vecD<U>& xc, const vecD<U>& bc) const		//! correct solution.
		{ x += P*xc; sm->post(x,A,bc); }
	vecD<U>	solveDense(const vecD<U>& bc) const					//! dense solve of A*x = b on the coarsest grid.
		{ return solve(asDense(A),bc); }
	void	init(const T th, matS<U>& An, vecD<U>& bn);
};

template <typename T, typename U>
void level<T,U>::init(const T th, matS<U>& An, vecD<U>& bn)
//! initialize P and R at this level.
{	const vecD<U> d = A.d(); matS<U> S(A);						// estimate omega
 	for allRows(A,r) for allCols(A,t) { S.x[t] = A.x[t]/d(r); };
	const U omega = U(4.0/3.0)/S.spectralRadius();
	bn = setOperators(th,omega); An = R*(A*P);
}

template <typename T, typename U>
vecD<U> level<T,U>::setOperators(const T th, const U omega)
//! computes prolongation P and restriction R.
{	const ivecD agg = A.strength(th).aggregate();
	unsigned int na = ITOU(agg(agg.amax())+1), nz = agg.N, nu = 0;
	for (unsigned int r = 0; r < agg.N; r++) if (agg(r) < 0) nu++;			// count isolated nodes
	matS<U> Q(nz,na,nz-nu); unsigned int t = 0;
	for (unsigned int r = 0; r < agg.N; r++) {
		Q.ri[r] = t; if (agg(r) < 0) continue;					// do not use isolated nodes
		Q.ci[t] = ITOU(agg(r)); Q.x[t] = b(r); t++; };
	assert(t == Q.nz); Q.ri[agg.N] = Q.nz; vecD<U> bn(na); bn = U(0); 
	for allRows(Q,r) for allCols(Q,t) { bn[Q.ci[t]] += norm(Q.x[t]); };		// compute norm over each aggregate
	for (unsigned int i = 0; i < bn.N; i++) { bn[i] = std::sqrt(bn[i]); };
	for allRows(Q,r) for allCols(Q,t) { Q.x[t] /= bn[Q.ci[t]]; };			// rescale columns of Q
	const vecD<U> d = A.d(); matS<U> S(A);
	for allRows(S,r) for allCols(S,t) { const unsigned int c = S.ci[t];
		S.x[t] = -omega*S.x[t]/d(c); if (r == c) S.x[t] += U(1); };		// I-omega*D^-1*A
	P = S*Q; R = trp(P); return bn;							// set prologation and restriction operators
}

//! Algebraic multigrid preconditioner.
/*!
 * For detailed information, refer to (\cite Vanek96, <a href="ref32.pdf" target="_blank">PDF</a>)
*/

template <typename T, typename U = T> class precAMG : public precS<T,U> {
	std::vector<level<T,U>> levels;		//!< vector of resolution levels
	void	solve(vecD<U>& x, const vecD<U>& b, const unsigned int i = 0) const	//! computes approximate solution.
		{ if (i+1 == levels.size()) x = levels[i].solveDense(b);		// dense solve on coarse grid
		  else { const vecD<U> bc = levels[i].coarse(x,b);			// restrict to coarse grid
			vecD<U> xc; solve(xc,bc,i+1);					// compute coarse grid solution
			levels[i].correct(x,xc,b); } }					// apply coarse grid correction
public:
	precAMG(const matS<U>& A, const T th = T(0), const amgSmootherType t = amgSmootherType::polynomial)
		//! constructs an AMG preconditioner with drop level tol.
		 : precS<T,U>(), levels()
		{ levels.reserve(20); levels.push_back({A,t});				// avoid reallocations which force matrix copies
		  while (levels.back().dim() > 50) {					// for all levels with matrices greater than 50 rows
			matS<U> An; vecD<U> bn; levels.back().init(th,An,bn);		// init level and generate sub-sampled matrices
			levels.push_back({An,bn,t}); };					// allocate new level
		  checkFpException();
	 }
	vecD<U>	operator()(const vecD<U>& b) const					//! applies AMG preconditioner using vector b.
		{ vecD<U> x(b.N); x = U(0); solve(x,b); return x; }
};

//! Incomplete LU preconditioner.
/*!
 * For detailed information, refer to Saad Y (1992)
   ILUT: A dual threshold incomplete LU factorization, Algorithm 3.2 
*/

template <typename T, typename U = T> class precILU : public precS<T,U> {       
	using ivp = std::pair<unsigned int,U>;	// stores a matrix element
	using row = std::vector<ivp>;		// stores a sparse row

	std::vector<row> L;			//!< left (lower) elements
	std::vector<row> R;			//!< right (upper) elements
	vecD<U> D;				//!< inverse diagonal elements

	row	mergeRows(const row& w, const row& u, const U p) const			//! combines current row w with right row u and pivot p.
		{ unsigned int wi = 0, ui = 0, nw = w.size(), nu = u.size(); row z;
		  while (true) {
			if (wi < nw && ui < nu) {
				if (w[wi].first < u[ui].first) {
					z.push_back(w[wi]); wi++; }
				else if (w[wi].first == u[ui].first) {
					z.push_back({w[wi].first,w[wi].second-p*u[ui].second});
					wi++; ui++; }
				else { z.push_back({u[ui].first,-p*u[ui].second});
					ui++; } }
			else if (wi == nw && ui < nu) {
				z.push_back({u[ui].first,-p*u[ui].second}); ui++; }
			else if (wi < nw && ui == nu) { z.push_back(w[wi]); wi++; }
			else return z; } }
	void	truncRow(row& w, const unsigned int nf) const				//! keeps nf absolute largest elements in w.
		{ if (w.size() > nf) w.erase(w.begin()+nf,w.end()); 			// drop small elements
		  sort(w.begin(),w.end()); }	 					// sort by index
	void	insertSorted(row& w, const ivp& e)
		{ const T v = std::abs(e.second); if (v == T(0)) return;
		  unsigned int f = 0; while (f < w.size() && std::abs(w[f].second) > v) ++f;
		  ivp t = e; for (unsigned int j = f; j < w.size(); j++) std::swap(w[j], t);
		  w.push_back(t); }
	void	solveRow(vecD<U>& x, const row& w, const unsigned int r) const		// solve row by back substitution
		{ for (const auto& e: w) x[r] -= e.second*x[e.first]; }
public:
	precILU(const matS<U>& A, const unsigned int nf, const T tol)
		//! constructs an ILU preconditioner with up to nf elements above drop level tol.
		 : precS<T,U>(), L(A.M), R(A.M), D(A.M)
		{ D = U(0); for allRows(A,r) { T tau{0}; row w,l,u;			// collect elements of row r in w
			for allCols(A,t) { const U v = A.x[t]; tau += SQR(v);
				w.push_back({A.ci[t],v}); }
			tau = std::sqrt(tau)*tol; sort(w.begin(),w.end());		// set drop level and sort elements in w by index
			for (unsigned int k = 0; k < w.size(); k++) {
				const unsigned int c = w[k].first; if (c >= r) break;	// iterate over lower diagonal parts of A:
				const U p = w[k].second/D[c]; assert(std::isfinite(p));	// divide by diagonal
				if (std::abs(p) > tau) { w[k].second = p;		// merge in elements larger than tau
					if (R[c].size()) w = mergeRows(w,R[c],p); }
				else { w.erase(w.begin()+k,w.begin()+k+1); k--; } };
			for (const auto& wi : w) { if (wi.second == U(0)) continue;	// split remaining entries into left, right and diagonal
				if (wi.first < r) insertSorted(l,wi);
				else if (wi.first == r) { D[r] = wi.second;		// save inverse of diagonal element
					if (isSmall<T>(wi.second)) throw rtException("precILU: zero diagonal at %u", int(r)); }
				else insertSorted(u,wi); };
			truncRow(l,nf); L[r] = l; truncRow(u,nf); R[r] = u; } }
	vecD<U>	operator()(const vecD<U>& b) const					//! apply ILU preconditioner to vector b.
		{ vecD<U> x = b;
		  for (unsigned int r = 0; r < D.N; r++) solveRow(x,L[r],r); 		// forward block L solve
		  for (unsigned int r = D.N-1; r != ~0u; r--) {				// backward block R solve
			solveRow(x,R[r],r); x[r] /= D(r); }	
		  return x; }
};

#undef allRows
#undef allCols

#endif
