#ifndef PSOLVER_H
#define PSOLVER_H

/*
 *
 * psolver.h: partitioned vector, matrix and solver classes
 * BRIAN Software Package Version 3.0
 *
 * $Id: psolver.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (08/12/11): for BRIAN2.2 by FK
 * 0.40 (16/12/13): documented
 * 0.50 (17/11/14): released for BRIAN2.7
 * 0.60 (22/11/14): pprecILU corrected
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements partitioned preconditioner and solver classes.
*/

#define allRows(A, r)	(unsigned int r = 0; r < A.p.len(); r++)
#define sallRows(A, r)	(unsigned int r = 0; r < A.M; r++)
#define allCols(A, t)	(unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)

//! Abstract base class of all partitioned preconditioners.

template<typename T> class pprecS {
public:
	pprecS()									//! constructs a partitioned preconditioner.
		{ }
	virtual ~pprecS()
		{ }
        virtual pvecD<T> operator()(const pvecD<T>& b) const				//! applies preconditioner to partitioned vector b.
		{ return b; };								// returns identity
};

//! Implements a partitioned Jacobi preconditioner.

template<typename T> class pprecJacobi : public pprecS<T> {
	const pvecD<T> d;			//!< stores inverse of diagonal
public:
	pprecJacobi(const pmatS<T> &A)							//! constructs a partitioned Jacobi preconditioner.
		: pprecS<T>(), d(A.idia()) { };
        pvecD<T> operator()(const pvecD<T>& b) const					//! applies Jacobi preconditioner to vector b.
		{ return d*b; };
};
     
//! Implements a partitioned block Jacobi preconditioner.

template<typename T> class pprecBlockJacobi : public pprecS<T> {
	const unsigned int bs;			//!< block size
	std::vector<matD<T>> bm;		//!< local inverses
	matD<T>	getBlockInv(const pmatS<T>& A, const unsigned int b) const
		{ matD<T> Ai(bs,bs); Ai = T(0); const unsigned int s = b*bs;
		  for (unsigned int r = 0; r < bs; r++) { 
			for (unsigned int t = A.ri[s+r]; t < A.ri[s+r+1]; t++) {
				const unsigned int c = A.ci[t];
				if (c >= s && c < s+bs) Ai(r,c-s) = A.x[t]; } };
		  return inv(Ai); };
public:
	pprecBlockJacobi(const pmatS<T>& A, const unsigned int _bs)			//! constructs a block Jacobi preconditioner.
		 : pprecS<T>(), bs(_bs), bm(A.p.len()/bs)
		{ for (unsigned int i = 0; i < bm.size(); i++) bm[i] = getBlockInv(A,i);
		  checkFpException(); }
        pvecD<T> operator()(const pvecD<T>& b) const					//! applies block Jacobi preconditioner to vector b.
		{ pvecD<T> x(b.p); vecD<T> bi(bs);
		  for (unsigned int i = 0; i < bm.size(); i++) {
		 	for (unsigned int j = 0; j < bs; j++) bi[j] = b(i*bs+j);
			vecD<T> xi = bm[i]*bi;
		 	for (unsigned int j = 0; j < bs; j++) x[i*bs+j] = xi[j]; };
		  return x; };
};

//! partitioned Jacobi smoother for pprecAMG.

template <typename T> class pJacobi {
	const pvecD<T> inv;			//!< stores inverse of diagonal
public:
	pJacobi<T>(const pmatS<T>& A, const T omega)					//! construct a partitioned Jacobi smoother for sparse matrix A.
		: inv(A.idia()*omega) { }
	void	pre(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies pre-multiplication x = A * b.
		{ (void)A; x = inv*b; }
	void	post(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies post-multiplication x = A * b.
		{ pvecD<T> y = b-A*x; x -= inv*y; }
};

//! partitioned Gauss-Seidel smoother for pprecAMG.

template <typename T> class pGaussSeidel {
	const unsigned int nit;			//!< number of GS iterations

	void	sweep(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b,
			const int s, const int e, const int i) const			//! applies a Gauss-Seidel sweep x = A * b.
		{ x.repartition(A.p); x.sync();
		  for (int r = s; r != e; r += i) { T d{0}, rs{0};
			for allCols(A,t) { int c = A.ci[t]; T v = A.x[t];
				if (c == r) d = v; else rs += v*x[c]; };
			if (d != T(0)) { x[r] = (b(r)-rs)/d; } } }
	void	forward(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const	//! applies one forward Gauss-Seidel sweep
		{ sweep(x,A,b,0,A.p.len(),1); }
	void	backward(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const	//! applies one backward Gauss-Seidel sweep
		{ sweep(x,A,b,A.p.len()-1,-1,-1); }
	void	cycle(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies one sweep cycle.
		{ for (unsigned int it = 0; it < nit; it++) {
			forward(x,A,b); backward(x,A,b); } }
public:
	pGaussSeidel(const unsigned int n = 1)						//! constructs an empty Gauss-Seidel smoother.
		: nit(n) { }
	void	pre(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies pre-multiplication x = A * b.
		{ cycle(x,A,b); }
	void	post(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies post-multiplication x = A * b.
		{ cycle(x,A,b); }
};

//! partitioned polynomial smoother for pprecAMG.

template <typename T> class ppolynomial {
	const vecD<T> def;			//!< stores coefficients of the polynomial 
	vecD<T>	setCoefs(const T rho, const T lo = T(1.0/30.0), const T up = T(1.1)) const //! determines coefficients of the polynomial.
		{ const unsigned int n = 3; vecD<T> _def(n);
		  const T x0 = lo*rho, x1 = up*rho; vecD<T> rt(n);			// Chebyshev roots for the interval [-1,1]
		  for (unsigned int i = 0; i < n; i++)
			rt[i] = T(cos(M_PI*(T(i)+T(0.5))/n));
		  for (unsigned int i = 0; i < n; i++)					// Chebyshev roots for the interval [x0,x1]
			rt[i] = T(0.5*(x1-x0)*(T(1)+rt[i])+x0);
		  _def[0] = T(-1.0); _def[1] = rt[0]+rt[1]+rt[2];
		  _def[2] = -(rt[0]*rt[1]+rt[1]*rt[2]+rt[2]*rt[0]); 
		  T f = -rt[0]*rt[1]*rt[2]; _def /= f; return _def; }
public:
	ppolynomial<T>(const T rho)							//! constructs a polynomial smoother for eigenvalue rho.
		: def(setCoefs(rho)) { }
	void	pre(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies pre-multiplication x = A * b.
		{ x = b*def(0); for (unsigned int i = 1; i < def.N; i++)
			x += A*asDense(x)+b*def(i); }
	void	post(pvecD<T>& x, const pmatS<T>& A, const pvecD<T>& b) const		//! applies post-multiplication x = A * b.
		{ pvecD<T> r = b-A*x, h = r*def(0); h.repartition(A.p);
		  for (unsigned int i = 1; i < def.N; i++) {
			h = A*h+r*def(i); x += h; } }
};

//! Stores a level for a partitioned algebraic multigrid preconditioner.

template <typename T>
class plevel {
       	const pmatS<T> A;			//!< sparse matrix
	const T	rho;				//!< spectral radius of A
      	ppolynomial<T> sm;			//!< smoothing operator
       	pmatS<T> R;				//!< restriction operator
       	pmatS<T> P;				//!< prolongation operator
       	pvecD<T> b;				//!< right-hand side
 
	void	prolongatorSmoother(pvecD<T>& bn, const T th, const T omega);
public:
	plevel()									//! constructs an empty level.
		: A(), rho(T(0)), sm(), R(), P(), b() { }
	plevel(const pmatS<T>& _A)							//! constructs a level for matrix A.
		: A(_A), rho(A.spectralRadius()),sm(ppolynomial<T>(rho)), R(), P(), b(A.p)
		{ b = 1; }
	plevel(const pmatS<T>& _A, const pvecD<T>& _b)					//! constructs a level for matrix A and RHS b.
		 : A(_A), rho(A.spectralRadius()),sm(ppolynomial<T>(rho)),
		R(), P(), b(_b) { }
	unsigned int dim() const							//! returns number of rows.
		{ return A.m(); }
	unsigned int nz() const 							//! returns number of non-zero elements in A.
		{ return A.nz; }
	pvecD<T> coarse(pvecD<T>& x, const pvecD<T>& b) const				//! reduces solution.
		{ sm.pre(x,A,b); pvecD<T> r = b-A*x; return R*asDense(r); }
	void	correct(pvecD<T>& x, const pvecD<T>& xc, const pvecD<T>& b) const	//! corrects solution.
		{ x += P*asDense(xc); sm.post(x,A,b); }
	void	solveDense(pvecD<T>& x, const pvecD<T>& b) const			//! solves of A*x = b on the coarsest grid (at root).
		{ vecD<T> r = solve(asDense(A),b.gather()); x.scatter(r); }
	void	init(const T th, pmatS<T>& An, pvecD<T>& bn);
};

template <typename T>
void plevel<T>::init(const T th, pmatS<T>& An, pvecD<T>& bn)
//! initialize this level.
{	T omega; { pvecD<T> d = A.d(); pmatS<T> DiA(A);					// estimate omega
 	  for allRows(A,r) for allCols(A,t) DiA.x[t] = A.x[t]/d[r];
	  omega = T(4.0/3.0)/DiA.spectralRadius(); };					// delete DiA after estimation
	prolongatorSmoother(bn, th, omega); 
	R = trp(P); An = R*(A*P); bn.repartition(An.p);
}

template <typename T>
void plevel<T>::prolongatorSmoother(pvecD<T>& bn, const T th, const T omega)
//! computes prolongation and restriction.
{	// takes the input aggregate vector agg and the RHS b
	// on output, Q contains columns of vertices per aggregate, scaled by columnwise norm
	pmatS<T> St = A.strength(th); pivecD agg = St.aggregate();
	unsigned int na = ITOU(agg.max()+1), nu = 0, l = A.p.len();
	for (unsigned int i = 0; i < l; i++) if (agg(i) < 0) nu++;
	pmatS<T> Q(A.M,na,l-nu,A.p); unsigned int t = 0; b.sync();
	for (unsigned int r = 0; r < l; r++) {
		Q.ri[r] = t; if (agg[r] < 0) continue;
		Q.ci[t] = ITOU(agg(r)); Q.x[t] = b.x[r]; t++; };
	Q.ri[l] = Q.nz; Q.assemble();
	pmatS<T> Qt = trp(Q); bn.repartition(Qt.p); bn = T(0);
	for allRows(Qt,r) for allCols(Qt,t) bn[r] += SQR(Qt.x[t]); 			// compute norm over each aggregate
	for (unsigned int i = 0; i < Qt.p.len(); i++) bn[i] = std::sqrt(bn[i]);
	vecD<T> tbn = asDense(bn);
	for allRows(Q,r) for allCols(Q,t) Q.x[t] /= tbn[Q.l2g(t)];			// rescale columns of Q
	pmatS<T> S(A); vecD<T> d = asDense(S.d());
	const unsigned int s = S.p.st();
	for allRows(S,r) for allCols(S,t) { const unsigned int c = S.l2g(t); 
		T v = omega*S.x[t]/d[c], I = (r+s == c)? T(1): T(0); S.x[t] = I-v; };	// I-omega*D^-1*A
	S.assemble(); P = S*Q;
}

//! Algebraic multigrid preconditioner.
/*!
 * For detailed information, refer to (\cite Vanek96, <a href="ref32.pdf" target="_blank">PDF</a>)
*/

template <typename T> class pprecAMG : public pprecS<T> {
private:
	std::vector<plevel<T>> levels;		//!< vector of resolution levels
	void	solve(pvecD<T>& x, const pvecD<T>& b, const unsigned int i) const	//! compute approximate solution.
		{ if (i+1 == levels.size()) { levels[i].solveDense(x,b); }		// dense solve on coarse grid
		  else { pvecD<T> bc = levels[i].coarse(x,b);				// restrict to coarse grid
			pvecD<T> xc(bc.p); solve(xc,bc,i+1);				// compute coarse grid solution
			levels[i].correct(x,xc,b); } };					// apply coarse grid correction
public:
	pprecAMG(const pmatS<T>& _A, const T th = 0, const bool verbose = false)	//! construct a partitioned AMG preconditioner with drop level tol.
		 : pprecS<T>(),levels() 
		{ levels.reserve(20); levels.push_back({_A});				// avoid reallocations which force matrix copies
		  while (levels.back().dim() > 10*_A.p.nprocs) {
			pmatS<T> An; pvecD<T> bn;
			if (verbose && _A.p.isRoot())
				printf("level %zd %u %u\n", levels.size(), levels.back().dim(), levels.back().nz());
			levels.back().init(th, An, bn);
			levels.push_back({An, bn}); };
		  if (verbose && _A.p.isRoot())
				printf("level %zd %u\n", levels.size(), levels.back().dim()); };
	pvecD<T> operator()(const pvecD<T>& b) const					//! applies AMG preconditioner to vector x.
		{ pvecD<T> x(b.p); x = 0; solve(x,b,0); return x; };
};

//! Partitioned incomplete LUT preconditioner.

template <typename T> class pprecILU : public pprecS<T> {       
	const unsigned int N;			//!< number of columns
	std::vector<vecS<T>> L;			//!< left elements
	std::vector<vecS<T>> R;			//!< right elements
	pvecD<T> D;				//!< diagonal elements
	void	qswap(vecD<T>& a, uvecD& v, const unsigned int i,
			const unsigned int j) const   					//! swaps elements i and j in vector (a,v).
		{ std::swap(a[i],a[j]); std::swap(v[i],v[j]); }
	void	qsplit(vecD<T>& a, uvecD& v, const unsigned int n,
			const unsigned int t) const  					//! sorts vector (a,v) for the t (of n) most important elements.
		{ unsigned int f = 0, l = n-1; if (n == 0 || t < f || t > l) return;
		  while (true) { unsigned int m = f; T key = a[m];			// outer loop -- while m != t
			for (unsigned int i = f+1; i <= l; i++)
				if (a[i] > key) { m++; qswap(a,v,m,i); };
			qswap(a,v,m,f); if (m == t) return;				// test for while loop
			if (m > t) l = m-1; else f = m+1; } }
	vecS<T> allocateRow(const vecD<T>& x, const uvecD& u, const unsigned int l,
		const unsigned int n, const unsigned int o) const 			//! allocates sparse row from vector (x,u) with max n elements.
		{ vecD<T> xi(l); uvecD ci(l); vecS<T> row(n);
		  for (unsigned int i = 0; i < l; i++) {
			xi[i] = std::abs(x(i+o)); ci[i] = i+o; };
		  qsplit(xi, ci, l, n);
		  for (unsigned int i = 0; i < n; i++) {
			row.u[i] = u(ci[i]); assert(u(ci[i]) != ~0u);
			row.x[i] = x(ci[i]); }
		  return row; }
	void	factorize(const pmatS<T>& A, const unsigned int nf, const T tol,
		const std::vector<vecS<T>>& RL, const vecD<T>& DL)			//! factorizes this partition of A into L,D,R.
		{ uvecD u(A.N), col(A.N); vecD<T> x(A.N); u = ~0u;
#define fillIn(a) u[c] = a; col[a] = c; x[a] = xi; 
		  const unsigned int s = A.p.st();	
		  for allRows(A,r) { const unsigned int rs = r+s; 			// set col for all rows...
			T rn = 0; unsigned int lu = 0, ll = 0;
			col[rs] = rs; u[rs] = rs; x[rs] = T(0);				// set diagonal element
			for allCols(A,t) { unsigned int c = A.l2g(t); T xi = A.x[t]; // fill L- and R-part of row r
				rn += std::abs(xi);
				if (c == rs) x[rs] = xi;
				else if (c < rs) { fillIn(ll); ll++; }
				else { fillIn(rs+lu+1); lu++; } }
			if (isSmall<T>(rn)) throw rtException("Zero row %u", int(rs));
			T nz = A.ri[r+1]-A.ri[r], lim = tol*rn/nz;
			for (unsigned int i = 0; i < ll; i++) {				// eliminate previous rows
				unsigned int ri = col[i], t = i;
				for (unsigned int k = i+1; k < ll; k++)			// determine smallest column index
					if (col[k] < ri) { ri = col[k]; t = k; };
				if (t != i) { 						// swaps
					unsigned int c = col[i]; col[i] = col[t];
					col[t] = c; u[ri] = i; u[c] = t;
					T xi = x[i]; x[i] = x[t]; x[t] = xi; }
				x[i] *= ri >= s? D[ri-s]: DL(ri); u[ri] = ~0u;
				const vecS<T>& row = ri >= s? R[ri-s]: RL[ri];		// get the multiplier
				for (unsigned int k = 0; k < row.N; k++) {
					unsigned int c = row.u[k]; T xi = -x[i]*row.x[k];
					if (std::abs(xi) < lim && u[c] == ~0u) continue; // disregard small fill-in element
					if (u[c] != ~0u) x[u[c]] += xi;
					else if (c < rs) { fillIn(ll); ll++; }		// fill in L
					else { fillIn(rs+lu+1); lu++; } } }		// fill in R
			u[rs] = ~0u;
			for (unsigned int i = 0; i < lu; i++) u[col[rs+i+1]] = ~0u;	// restore u
			if (isSmall<T>(std::abs(x[rs])))
				throw rtException("Zero diagonal %u", int(rs));
			D[r] = T(1.0)/x[rs]; 
			L[r] = allocateRow(x, col, ll, std::min(ll, nf),0);		// update/store row of L
			R[r] = allocateRow(x, col, lu, std::min(lu, nf),rs+1); } }	// update/store row of R
#undef fillIn
	void	recvFactorization(const int p, std::vector<vecS<T>>& RL,
			vecD<T>& d) const 						//! receives a factorization from process proc.
		{ unsigned int len = 0, i = 0; unsigned char* buf; vecS<T> e; 
		  buf = reinterpret_cast<unsigned char *>(mpiRecvBuf(len,p));		// receive rows R
		  while (i < len) { i += e.unpack(buf+i); RL.push_back(e); }; delete [] buf; // unpack each row and store in RL
		  buf = reinterpret_cast<unsigned char *>(mpiRecvBuf(len,p));		// receive diagonal d
		  memcpy(d.x,buf,len); delete [] buf; }					// copy buffer to global diagonal
	void	sendFactorization(const int p, const std::vector<vecS<T>>& RL,
			vecD<T>& d) const						//! sends a factorization to process proc.
		{ unsigned int len = 0, t = 0, n = ITOU(RL.size());			// compute message length
		  for (unsigned int i = 0; i < n; i++) len += RL[i].size();		// sum up size of each row
		  unsigned char *buf = new unsigned char [len]; 			// allocate buffer
		  for (unsigned int i = 0; i < n; i++) t += RL[i].pack(buf+t);		// pack rows into buffer
		  mpiSendBuf(buf,len,p); delete [] buf;					// send RL to next process
		  mpiSendBuf(d.x,d.N*sizeof(T),p); }					// send d to next process
	void	solveRow(pvecD<T>& v, const vecS<T>& r, const unsigned int j) const 	//! solves in row (r,j) using vector v.
		{ for (unsigned int i = 0; i < r.N; i++) {
			unsigned int c = r.u[i]; if (c == ~0u) continue;
			v[j] -= v[c]*r.x[i]; } }
public:
	pprecILU(const pmatS<T>& A, const unsigned int nf, const T tol)			//! constructs an ILU preconditioner up to nf places or drop level tol.
		 : pprecS<T>(), N(A.p.len()), L(N), R(N), D(A.p)
		{ { std::vector<vecS<T>> RL; vecD<T> DL(A.M);
		    if (A.p.myid > 0) recvFactorization(A.p.myid-1,RL,DL);		// receive factorization from process on the left
		    factorize(A,nf,tol,RL,DL);						// factorize this partition
		    if (A.p.myid < A.p.nprocs-1) {					// if there is a process on the right
			RL.insert(RL.end(),R.begin(),R.end()); 				// merge factorization
		 	const unsigned int s = D.p.st(), l = D.p.len();	
		  	for (unsigned int i = 0; i < l; i++) DL[i+s] = D(i);		// merge diagonal
			sendFactorization(A.p.myid+1,RL,DL); } };			// send RL and d to next process
		  uvecD map = A.g2lmap(); const unsigned int l = A.p.len();
		  for (unsigned int r = 0; r < l; r++) {				// map all L,R through g2lmap via u
			unsigned int *ul = L[r].u, *ur = R[r].u;
			for (unsigned int i = 0; i < L[r].N; i++) ul[i] = map(ul[i]);
			for (unsigned int i = 0; i < R[r].N; i++) ur[i] = map(ur[i]); }
		  checkFpException(); }
	pvecD<T> operator()(const pvecD<T>& b) const					//! applies ILU preconditioner to vector x.
		{ pvecD<T> x = b; x.sync();
		  for (unsigned int r = 0; r < N; r++) solveRow(x,L[r],r);		// block L solve
		  x.sync(); for (unsigned int r = N-1; r != ~0u; r--) {
			solveRow(x,R[r],r); x[r] *= D(r); };				// block R solve
		  return x; }
};

#undef allRows
#undef sallRows
#undef allCols

//! Abstract base class of all iterative sparse parallel solvers.

template<typename T> class psolver {
protected:
	static const unsigned int maxit = 10000;//!< max number of iterations
	const pmatS<T>&		A;		//!< ref to sparse matrix
	const pprecS<T>&	M;		//!< ref to preconditioner
	const T			eps;		//!< final residual tolerance
public:
	psolver(const pmatS<T>& _A, const pprecS<T>& _M = pprecS<T>(),
		const T _eps = T(1e-6))							//! construct a solver for matrix A and preconditioner M.
		 : A(_A), M(_M), eps(_eps) { clearFpException(); }
	virtual ~psolver() { }
	virtual void solve(pvecD<T>& x, const pvecD<T>& b, const bool verbose = false)	//! solve A x = b up to tolerance eps.
		 = 0;
	virtual pvecD<T> mv(const pvecD<T>& x) const					//! sparse matrix-vector product.
		{ return A*x; }
	virtual pvecD<T> prec(const pvecD<T>& x) const					//! apply preconditioner to vector x.
		{ return M(x); }
};

#define smat	psolver<T>::A			// to make code a bit more readable

//! Partitioned stabilized bi-conjugate gradient solver.

template<typename T> class pbicgstab : public psolver<T> {
public:
	pbicgstab(const pmatS<T>& _A, const pprecS<T>& _M, const T _eps = T(1e-6))	//! construct bicgstab solver for matrix A and preconditioner M.
		: psolver<T>(_A,_M,_eps) { }
	void	solve(pvecD<T>& x, const pvecD<T>& b, const bool verbose = false)	//! solve A x = b up to tolerance eps.
		{ x = T(0); T bn = norm(b); if (isSmall<T>(bn)) bn = T(1.0);
		  T alpha = T(1), omega = T(1), rho2 = T(1);
		  pvecD<T> r = b-this->mv(x), rt = r, p, v; T cv = norm(r)/bn;
		  for (unsigned int it = 0; it < this->maxit && cv > this->eps; it++) {
			T rho1 = dot(rt,r);
			if (isSmall<T>(rho1)) throw rtException("rho breakdown in bicgstab");
			if (it == 0) p = r;
			else { T t = (rho1/rho2)*(alpha/omega); p = r+(p-v*omega)*t; };
			pvecD<T> phat = this->prec(p); v = this->mv(phat);
			T rtv = dot(rt,v);
			if (isSmall<T>(rtv)) rtException("rtv breakdown in bicgstab");
			alpha = rho1/rtv; pvecD<T> s = r-v*alpha; x += alpha*phat;
			cv = norm(s)/bn; if (cv < this->eps) break;
			pvecD<T> shat = this->prec(s); pvecD<T> t = this->mv(shat);
			omega = dot(t,s)/dot(t,t); x += omega*shat; r = s-omega*t;
			cv = norm(r)/bn; rho2 = rho1;
			if (verbose && smat.p.isRoot()) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); } };
		  checkFpException();
		  if (cv > this->eps) throw rtException("no convergence in bicgstab"); }
};

//! Partitioned generalized minimum residual (GMRES) solver.
/*!
Implements the algorithm in [Figure 2.6](http://www.netlib.org/linalg/html_templates/node29.html) of:\n
Barrett et al., Templates for the Solution of Linear Systems.
Siam Book 1994
*/ 

template<typename T> class pgmres : public psolver<T> {
	unsigned int m;				//!< subspace dimension

	void	applyRotation(T& dx, T& dy, const T cs, const T sn) const		//! applies a Givens rotation.
		{ T t = cs*dx+conj(sn)*dy; dy = -sn*dx+cs*dy; dx = t; }
	void	generateRotation(const T dx, const T dy, T& cs, T& sn)  const		//! computes a Givens rotation.
		{ if (isSmall<T>(dy)) { cs = T(1); sn = T(0); }
		  else { T t1 = norm(dx), t2 = norm(dy), nm = std::sqrt(t1+t2);
			cs = dx/nm; sn = dy/nm; } }
	void	upperTriSolve(vecD<T>& s, const matD<T>& H, const unsigned int n) const	//! solves the final upper triangular system.
		{ for (unsigned int i = n-1; i != std::numeric_limits<unsigned int>::max(); i--) { T si = s[i];
			for (unsigned int j = i+1; j < n; j++) si -= H(i,j)*s[j];
  		  s[i] = si/H(i,i); } }
public:
	pgmres(const pmatS<T>& _A, const pprecS<T>& _M, const unsigned int _m = 10,
			const T _eps = T(1e-6))						//! construct GMRES solver for matrix A and preconditioner M.
		: psolver<T>(_A,_M,_eps), m(_m) { }
	void	solve(pvecD<T>& x, const pvecD<T>& b, const bool verbose = false)	//! solve A x = b up to tolerance eps.
		{ x = T(0); T bn = norm(b); if (isSmall<T>(bn)) bn = T(1);
		  matD<T> H(m+1, m); H = T(0); vecD<T> cs(m), sn(m), s(m+1);
		  cs = T(1); sn = T(0); s = T(0);
		  for (unsigned int it = 0; it < this->maxit; it++) {
			pvecD<T> r = this->prec(b-this->mv(x)); T rn = norm(r), nm = rn;
			if (verbose && smat.p.isRoot()) { printf("[%4d] %6.4e\n", it, rn/bn); fflush(stdout); };
			if (rn/bn <= this->eps) { checkFpException(); return; };
			std::vector<pvecD<T>> V(m+1); V[0] = r/T(rn); s[0] = T(rn);
			for (unsigned int i = 0; i < m; i++) {
				if (nm/bn <= this->eps) { m = i; break; };
				pvecD<T> v = this->prec(this->mv(V[i]));
				for (unsigned int j = 0; j <= i; j++) {
					H(j,i) = dot(V[j],v); v -= H(j,i)*V[j]; };
				H(i+1,i) = norm(v); V[i+1] = v/H(i+1,i);
				for (unsigned int j = 0; j < i; j++)
					applyRotation(H(j,i), H(j+1,i), cs[j], sn[j]);
				generateRotation(H(i,i), H(i+1,i), cs[i], sn[i]); 
				applyRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
				s[i+1] = -sn[i]*s[i]; s[i] = cs[i]*s[i];
				nm = std::abs(s[i+1]); };
			upperTriSolve(s,H,m);
			for (unsigned int i = 0; i < m; i++) x += V[i]*s[i]; };
		     checkFpException(); throw rtException("no convergence in pgmres"); }
};

//! Partitioned induced dimension reduction (IDR) solver.
/*!
 * For detailed information, refer to (\cite Gijzen11, <a href="ref37.pdf" target="_blank">PDF</a>)
*/
/*!
Implements the algorithm in:\n
van Gijzen & Sonneveld (2011)\n
Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits Bi-orthogonality Properties.\n
[ACM ToMS (38) 1-19](http://ta.twi.tudelft.nl/nw/users/gijzen/idrs_toms.pdf).
*/

template<typename T> class pidrs : public psolver<T> {
	const unsigned int m;			//!< subspace dimension

	vecD<T>	orthogonalize(const std::vector<pvecD<T>>& V, pvecD<T>& v) const	//! orthogonalize vector v against basis V
		{ unsigned int n = ITOU(V.size()); vecD<T> s(n); s = T(0);		// return something useful if n == 0
		  for (unsigned int i = 0; i < n; i++) s[i] = dot(V[i],v);
		  for (unsigned int i = 0; i < n; i++) v -= V[i]*s[i];
		  return s; }
	void	randomInit(pvecD<T>& v)
		{ uniformDist<T> ud; for (unsigned int i = 0; i < v.N; i++) v[i] = ud(); 
		  T n = norm(v); v /= T(n); }
public:
	pidrs(const pmatS<T>& _A, const pprecS<T>& _M, const unsigned int _m = 4,
		const T _eps = T(1e-6))							//! construct IDR solver for matrix A and preconditioner M.
		 : psolver<T>(_A,_M,_eps), m(_m) { }
	void	solve(pvecD<T>& x, const pvecD<T>& b, const bool verbose = false)	//! solve A x = b up to tolerance eps.
		{ x = T(0); T bn = norm(b); if (isSmall<T>(bn)) bn = T(1.0);
		  T omega = T(1.0); T kappa = T(0.7);					// leave kappa fixed at setting found in the paper
		  pvecD<T> r = b-this->mv(x); T cv = norm(r)/bn;
		  std::vector<pvecD<T>> V; pvecD<T> v(b.p); vecD<T> f(m);
		  randomInit(v); v /= T(norm(v)); V.push_back(v);			// init test space
		  for (unsigned int i = 1; i < m; i++) {
			randomInit(v); orthogonalize(V,v); v /= T(norm(v)); V.push_back(v); };
		  matD<T> M(m,m); M.id(); std::vector<pvecD<T>> G, H; pvecD<T> a(b.p); a = 0;
		  for (unsigned int i = 0; i < m; i++) { G.push_back(a); H.push_back(a); };
		  for (unsigned int it = 0; it < this->maxit && cv > this->eps; it++) {
			for (unsigned int i = 0; i < m; i++) f[i] = dot(V[i],r);	// generate rhs for small system
			for (unsigned int i = 0; i < m; i++) {
				vecD<T> c = ::solve(M,f); v = r;			// solve small system
				for (unsigned int j = i; j < m; j++) v -= c[j]*G[j];	// compute basis vectors g,h
				pvecD<T> h = this->prec(v)*omega;
				for (unsigned int j = i; j < m; j++) h += c[j]*H[j];
				pvecD<T> g = this->mv(h);
				for (unsigned int j = 0; j < i; j++) {			// bi-orthogonalize the basis vectors
					T t = dot(V[j],g)/M(j,j); g -= t*G[j]; h -= t*H[j]; };
				H[i] = h; G[i] = g;					// and save them
				for (unsigned int j = i; j < m; j++) 			// update column of M = (V'*G)' 
					M(j,i) = dot(V[j],g);
				T t = f[i]/M(i,i); r -= t*g; x += t*h;			// make r orthogonal to M
				cv = norm(r)/bn; if (cv < this->eps) break;
				for (unsigned int j = 0; j <= i; j++) f[j] = T(0);	// f = M'*r
				for (unsigned int j = i+1; j < m; j++) f[j] -= t*M(j,i); };
			pvecD<T> t = this->prec(r); v = this->mv(t);			// entering G_j+1
			T vr = dot(r,v); omega = vr/dot(v,v);				// update omega
			T rho = std::abs(vr/(norm(v)*norm(r))); 
			if (rho < kappa) omega *= kappa/rho;
			r -= omega*v; x += omega*t; cv = norm(r)/bn;			// update solution & check residual
			if (verbose && smat.p.isRoot()) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); } };
		  checkFpException();
		  if (cv > this->eps) throw rtException("no convergence in pidrs"); }
};

//! Partitioned least-squares minimum residual (LSMR) solver.
/*!
 * For detailed information, refer to (\cite Fong11, <a href="ref35.pdf" target="_blank">PDF</a>)
*/
template <typename T> class plsmr : public psolver<T> {
	const T	damp;			//!< damping factor

	T	sign(const T a) { return a < T(0)? T(-1): T(1); }			//! returns sign of a as int.
	void	symOrtho(const T a, const T b, T& c, T& s, T& r)			//! constructs a rotation.
		{ if (isSmall<T>(b)) { c = sign(a); s = 0; r = std::abs(a); }
		  else if (isSmall<T>(a)) { c = 0; s = sign(b); r = std::abs(a); }
		  else if (std::abs(b) > std::abs(a)) {
			T tau = a/b; s = sign(b)/T(std::sqrt(1+SQR(tau))); c = s*tau; r = b/s; }
		  else { T tau = b/a; c = sign(a)/T(std::sqrt(1+SQR(tau))); s = c*tau; r = a/c; } }
public:
	plsmr(const pmatS<T>& _A, const T d = 0, const T _eps = T(1e-6))			//! construct LSMR solver for matrix A.
		 : psolver<T>(_A,pprecS<T>(),_eps), damp(d) { }
	void	solve(pvecD<T>& x, const pvecD<T>& b, const bool verbose = false)	//! solve A x = b up to tolerance eps.
		{ // see Fong & Saunders, SIAM J Sci Comp March 2011.
		  const pmatS<T> At = trp(this->A); const unsigned int n = this->A.N; 
		  pvecD<T> u = b; x.repartition(this->A.p); x = T(0); 
		  pvecD<T> v(n); v = T(0); T beta = norm(u), alpha = T(0);		// 1- section 5.2
		  if (beta > T(0)) { u /= T(beta); v = At*u; alpha = norm(v); };
		  if (alpha > T(0)) v /= T(alpha);
		  T zetabar = alpha*beta; 
		  T alphabar = alpha, normA2 = alpha*alpha, betadd = beta, betad = T(0);
		  T cbar = T(1), sbar = T(0), rho = T(1), rhobar = T(1), rhodold = T(1);
		  T tautildeold = T(0), thetatilde = T(0), zeta = T(0), d = T(0);
		  pvecD<T> h = v; pvecD<T> hbar(n); hbar = T(0); 
		  if (zetabar == T(0)) throw rtException("Zeta is zero");
		  for (unsigned int it = 0; it < this->maxit; it++) {			// 2- repeat forever
	        	u = this->mv(v)-u*T(alpha); beta = norm(u);			// 3- bi-diagonalization
        		if (beta > T(0)) { u /= T(beta); v = At*u-v*T(beta); };
 			alpha = norm(v); if (alpha > 0) v /= T(alpha);
			T chat, shat, alphahat, rhoold = rho, c, s;
			symOrtho(alphabar, damp, chat, shat, alphahat); 		// 4- construct rotation phat
			symOrtho(alphahat, beta, c, s, rho);				// 5- construct and apply rotation p
			T thetanew = s*alpha; alphabar = c*alpha;
			T rhobarold = rhobar, zetaold = zeta, thetabar = sbar*rho;
			symOrtho(cbar*rho, thetanew, cbar, sbar, rhobar);		// 6- construct and apply rotation pbar
			zeta = cbar*zetabar; zetabar = -sbar*zetabar;
			hbar = h-T(thetabar*rho/(rhoold*rhobarold))*hbar;		// 7- update h, x, hbar
			x += T(zeta/(rho*rhobar))*hbar; h = v-T(thetanew/rho)*h;
			T betaacute = chat*betadd, betacheck = -shat*betadd;		// 8- apply rotation phat, p
			T betahat = c*betaacute; betadd = -s*betaacute;
			T ctildeold, stildeold, rhotildeold, thetatildeold = thetatilde;
        		symOrtho(rhodold, thetabar, ctildeold, stildeold, rhotildeold);	// 9- construct and apply rotation ptilde
 			thetatilde = stildeold*rhobar; rhodold = ctildeold*rhobar;
			betad = -stildeold*betad+ctildeold*betahat;
			tautildeold = (zetaold-thetatildeold*tautildeold)/rhotildeold;	// 10- update tau
			T taud = (zeta-thetatilde*tautildeold)/rhodold; d += SQR(betacheck); 
			T normr = T(std::sqrt(d+SQR(betad-taud)+SQR(betadd)));		// 11- compute residuals
			if (normr == 0) throw rtException("System is singular");
			normA2 += SQR(beta)+SQR(alpha);
			T rtol = T(std::abs(zetabar)/(std::sqrt(normA2)*normr));		// 12- check termination criterion
			if (verbose) { printf("[%4d] %6.4e\r", it, rtol); fflush(stdout); };
			if (rtol < this->eps) { checkFpException(); return; } };
		  checkFpException(); throw rtException("No convergence in lsmr"); }
};

#endif

