#ifndef SOLVER_H
#define SOLVER_H

/*
 *
 * solver.h: iterative solvers for sparse linear systems
 * BRIAN Software Package Version 3.0
 *
 * $Id: solver.h 508 2017-03-26 20:13:21Z frithjof $
 *
 * 0.10 (26/05/11): for BRIAN2 by FK
 * 0.11 (17/06/11): preconditioners introduced
 * 0.12 (03/12/11): switched to CSR
 * 0.13 (05/12/11): replace LU code with dense BLAS call via matD
 * 0.14 (10/12/11): reimplementation of ainv classes
 * 0.20 (23/09/12): IDRs implemented, revised Ainv, split preconditioner, base class enhanced
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.41 (25/09/14): checked LDL&  lsmr, minor corrections
 * 0.50 (24/11/14): cg added, idrs modified to vectors
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Provides iterative solvers for sparse linear systems.
*/

//! Abstract base class of all iterative sparse solvers.

template<typename T, typename U = T> class solver {
protected:
	static const unsigned int maxit = 10000; //!< max number of iterations
	const matS<U>&		A;		//!< ref to sparse matrix
	const precS<T,U>&	M;		//!< ref to preconditioner
	const T			eps;		//!< final residual tolerance
public:
	solver(const matS<U>& _A, const precS<T,U>& _M = precS<T,U>(), const T _eps = T(1e-6))	//! constructs a solver for matrix A and preconditioner M.
		 : A(_A), M(_M), eps(_eps) { clearFpException(); }
	virtual ~solver() { }
	virtual vecD<U> solve(const vecD<U>&, const bool) = 0;				//! solves A x = b up to tolerance eps.
	virtual vecD<U> mv(const vecD<U>& x) const					//! computes a sparse matrix-vector product.
		{ return A*x; }
	virtual vecD<U> tm(const vecD<U>& x) const					//! computes a sparse transpose matrix-vector product.
		{ return A.tm(x); }
	virtual vecD<U> prec(const vecD<U>& x) const					//! applies preconditioner to vector x.
		{ return M(x); }
};

//! Conjugate gradient solver.
/*!
Implements the algorithm in [Figure 2.5](http://www.netlib.org/linalg/html_templates/node29.html) of:\n
Barrett et al., Templates for the Solution of Linear Systems.
Siam Book 1994
*/ 

template<typename T, typename U = T> class cg : public solver<T,U> {
public:
	cg(const matS<U>& _A, const precS<T,U>& _M, const T _eps = T(1e-6))		//! constructs cg solver for matrix A and preconditioner M.
		: solver<T,U>(_A,_M,_eps) { }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false)			//! solves A x = b up to tolerance eps.
		{ vecD<U> x(b.N); x = U(0);
		  const T bn = norm(b); if (bn < solver<T,U>::eps) return x;
		  vecD<U> r = b, z = solver<T,U>::prec(r), p = z; T rn = dot(r,z); 
		  for (unsigned int it = 0; it < solver<T,U>::maxit; it++) {
			const vecD<U> t = solver<T,U>::mv(p); const T alpha = rn/dot(t,p);
			x += alpha*p; r -= alpha*t; z = solver<T,U>::prec(r);
			const T new_rn = dot(r,z);
			if (it % 100 == 0) { const T cv = norm(solver<T,U>::mv(x)-b)/bn;
				if (verbose) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); };
				if (std::abs(cv) < solver<T,U>::eps) return x; }
			const T beta = new_rn/rn; rn = new_rn; p = z+beta*p; }
		  printf("cg: no convergence.\n"); return x; }
};

//! Stabilized bi-conjugate gradient solver.
/*!
Implements the algorithm in [Figure 2.10](http://www.netlib.org/linalg/html_templates/node32.html#941) of:\n
Barrett et al. (1994) Templates for the Solution of Linear Systems.
Siam Book.
*/ 

template<typename T, typename U = T> class bicgstab : public solver<T,U> {
public:
	bicgstab(const matS<U>& _A, const precS<T,U>& _M, const T _eps = T(1e-6))	//! constructs bicgstab solver for matrix A and preconditioner M.
		: solver<T,U>(_A,_M,_eps) { }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false)			//! solves A x = b up to tolerance eps.
		{ vecD<U> x(b.N), p; x = U(0); T beta{0}, omega{1}, alpha{0}, rho1{0};
		  const T bn = norm(b); if (bn < solver<T,U>::eps) return x;
		  vecD<U> r = b-solver<T,U>::mv(x); vecD<U> rt = r,v;
		  T cv = norm(r)/bn; if (cv < solver<T,U>::eps) return x;
		  for (unsigned int it = 0; it < solver<T,U>::maxit; it++) {
			const T rho = dot(rt,r);
			if (rho == T(0)) throw rtException("bicgstab: rho breakdown");
			if (it) { beta = (rho/rho1)*(alpha/omega); p = r+beta*(p-omega*v); }
			else p = r;
			const vecD<U> ph = solver<T,U>::prec(p);
			v = solver<T,U>::mv(ph); alpha = rho/dot(rt,v);
			const vecD<U> s = r-alpha*v;
			if (norm(s) < solver<T,U>::eps) { x += alpha*ph; return x; }
			const vecD<U> sh = solver<T,U>::prec(s);
			const vecD<U> t = solver<T,U>::mv(sh); omega = dot(t,s)/dot(t,t);
			x += alpha*ph+omega*sh; r = s-omega*t;
			cv = norm(r)/bn; rho1 = rho;
			if (verbose) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); };
			if (cv < solver<T,U>::eps) return x;
			if (omega == T(0)) throw rtException("bicgstab: omega breakdown"); };
		  printf("bicgstab: no convergence.\n"); return x; }
};


//! Generalized minimum residual (GMRES) solver.
/*!
Implements the algorithm in [Figure 2.6] (http://www.netlib.org/linalg/html_templates/node29.html) of:\n
Barrett et al., Templates for the Solution of Linear Systems.
Siam Book 1994
*/ 

template<typename T, typename U = T> class gmres : public solver<T,U> {
	unsigned int m;				//!< subspace dimension

	void	applyRotation(T& dx, T& dy, const T cs, const T sn) const		//! applies a Givens rotation.
		{ const T t = cs*dx+sn*dy; dy = -sn*dx+cs*dy; dx = t; }
	void	generateRotation(const T dx, const T dy, T& cs, T& sn) const		//! computes a Givens rotation.
		{ if (dy == T(0)) { cs = T(1); sn = T(0); }
		  else if (std::abs(dy) > std::abs(dx)) { const T t = dx/dy;
			sn = T(1)/std::sqrt(T(1)+SQR(t)); cs = t*sn; }
		  else { const T t = dy/dx;
			cs = T(1)/std::sqrt(T(1)+SQR(t)); sn = t*cs; } }
	void	upperTriSolve(vecD<T>& s, const matD<T>& H) const			//! solves the final upper triangular system.
		{ for (unsigned int i = H.M-1; i != ~0u; i--) { T si = s[i];
			for (unsigned int j = i+1; j < H.N; j++) si -= H(i,j)*s[j];
  		  	s[i] = si/H(i,i); } }
public:
	gmres(const matS<U>& _A, const precS<T,U>& _M, const unsigned int _m = 10,
			const T _eps = T(1e-6))						//! constructs GMRES solver for matrix A and preconditioner M.
		: solver<T,U>(_A,_M,_eps), m(_m) { }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false)			//! solves A x = b up to tolerance eps.
		{ vecD<U> x(b.N); x = U(0);
		  const T bn = norm(b); if (bn < solver<T,U>::eps) return x;
		  matD<T> H(m+1,m); vecD<T> cs(m),sn(m),s(m+1);
		  H = T(0); cs = T(1); sn = T(0);
		  for (unsigned int it = 0; it < solver<T,U>::maxit; it++) {
			vecD<U> v = b-solver<T,U>::mv(x); const T cv = norm(v)/bn;
			const vecD<U> r = solver<T,U>::prec(v); const T rn = norm(r);
			if (verbose) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); };
			if (cv < solver<T,U>::eps) return x;
			std::vector<vecD<U>> V(m+1); V[0] = r/rn; s = T(0); s[0] = rn;
			for (unsigned int i = 0; i < m; i++) {
				v = solver<T,U>::prec(solver<T,U>::mv(V[i]));
				for (unsigned int j = 0; j <= i; j++) {
					H(j,i) = dot(V[j],v); v -= H(j,i)*V[j]; };
				H(i+1,i) = norm(v); V[i+1] = v/H(i+1,i);
				for (unsigned int j = 0; j < i; j++)
					applyRotation(H(j,i),H(j+1,i),cs(j),sn(j));
				generateRotation(H(i,i),H(i+1,i),cs[i],sn[i]); 
				applyRotation(H(i,i),H(i+1,i),cs(i),sn(i));
				applyRotation(s[i],s[i+1],cs(i),sn(i)); };
			upperTriSolve(s,H);
			for (unsigned int i = 0; i < m; i++) x += V[i]*s[i]; };
		  printf("gmres: no convergence\n"); return x; }
};


//! Least-squares minimum residual (LSMR) solver.
/*!
 * For detailed information, refer to (\cite Fong11, <a href="ref35.pdf" target="_blank">PDF</a>)
*/
template <typename T, typename U = T> class lsmr : public solver<T,U> {
	const T	damp;				//!< damping factor

	T	sign(const T a) { return a < T(0)? T(-1): T(1); }			//! returns sign of a as int.
	void	symOrtho(const T a, const T b, T& c, T& s, T& r)			//! constructs a rotation.
		{ if (isSmall<T>(b)) { c = sign(a); s = T(0); r = std::abs(a); }
		  else if (isSmall<T>(a)) { c = 0; s = sign(b); r = std::abs(a); }
		  else if (std::abs(b) > std::abs(a)) {
			const T tau = a/b; s = sign(b)/T(std::sqrt(T(1)+SQR(tau)));
			c = s*tau; r = b/s; }
		  else { const T tau = b/a; c = sign(a)/T(std::sqrt(T(1)+SQR(tau)));
			s = c*tau; r = a/c; } }
public:
	lsmr(const matS<U>& _A, const T d = 0, const T _eps = T(1e-6))			//! constructs LSMR solver for matrix A.
		 : solver<T,U>(_A,precS<T,U>(),_eps), damp(d) { }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false)			//! solves A x = b up to tolerance eps.
		{ unsigned int n = solver<T,U>::A.N; vecD<U> x(b.N), u = b; x = U(0);
		  const T bn = norm(b); if (bn < solver<T,U>::eps) return x;
		  vecD<U> v(n); v = U(0); T beta = norm(u), alpha = T(0);		// section 5.2, step 1
		  if (beta > T(0)) { u /= U(beta); v = solver<T,U>::tm(u); alpha = norm(v); };
		  if (alpha > T(0)) v /= U(alpha);
		  T zetabar = alpha*beta; if (zetabar == T(0)) throw rtException("lsmr: zeta is zero");
		  T alphabar = alpha, normA2 = alpha*alpha, betadd = beta, betad = T(0);
		  T cbar = T(1), sbar = T(0), rho = T(1), rhobar = T(1), rhodold = T(1);
		  T tautildeold = T(0), thetatilde = T(0), zeta = T(0), d = T(0);
		  vecD<U> h = v; vecD<U> hbar(n); hbar = U(0); 
		  if (zetabar == T(0)) throw rtException("lsmr: zeta is zero");
		  for (unsigned int it = 0; it < 50*solver<T,U>::maxit; it++) {		// 2- repeat forever
	        	u = solver<T,U>::mv(v)-u*U(alpha); beta = norm(u);		// 3- bidiagnoalization
        		if (beta > T(0)) { u /= U(beta); v = solver<T,U>::tm(u)-v*U(beta); };
 			alpha = norm(v); if (alpha > 0) v /= U(alpha);
			T chat, shat, alphahat, rhoold = rho, c, s;
			symOrtho(alphabar, damp, chat, shat, alphahat); 		// 4- construct rotation phat
			symOrtho(alphahat, beta, c, s, rho);				// 5- construct and apply rotation p
			T thetanew = s*alpha; alphabar = c*alpha;
			T rhobarold = rhobar, zetaold = zeta, thetabar = sbar*rho;
			symOrtho(cbar*rho, thetanew, cbar, sbar, rhobar);		// 6- construct and apply rotation pbar
			zeta = cbar*zetabar; zetabar = -sbar*zetabar;
			hbar = h-U(thetabar*rho/(rhoold*rhobarold))*hbar;		// 7- update h, x, hbar
			x += U(zeta/(rho*rhobar))*hbar; h = v-U(thetanew/rho)*h;
			T betaacute = chat* betadd, betacheck = -shat* betadd;		// 8- apply rotation phat, p
			T betahat = c*betaacute; betadd = -s*betaacute;
			T ctildeold, stildeold, rhotildeold, thetatildeold = thetatilde;
        		symOrtho(rhodold, thetabar, ctildeold, stildeold, rhotildeold);	// 9- construct and apply rotation ptilde
 			thetatilde = stildeold*rhobar; rhodold = ctildeold* rhobar;
			betad = -stildeold*betad+ctildeold*betahat;
			tautildeold = (zetaold-thetatildeold*tautildeold)/rhotildeold;	// 10- update tau
			T taud = (zeta - thetatilde*tautildeold)/rhodold; d += SQR(betacheck); 
			T normr = T(std::sqrt(d+SQR(betad-taud)+SQR(betadd)));		// 11- compute residuals
			if (normr == 0) throw rtException("lsmr: system is singular");
			normA2 += SQR(beta)+SQR(alpha);
//			T rtol = T(std::abs(zetabar)/(std::sqrt(normA2)*normr));		// 12- check termination criterion
			if (it % 100 == 0) { const T cv = norm(b-solver<T,U>::mv(x))/bn;
				if (verbose) { printf("[%4d] %6.4e\r", it, cv); fflush(stdout); };
				if (cv < solver<T,U>::eps) return x; } };
		  printf("lsmr: no convergence.\n"); return x; }
};
#undef allCols
#undef allRows
#undef allARows
#undef allACols

//! Induced dimension reduction (IDR) solver.
/*!
 * For detailed information, refer to (\cite Gijzen11, <a href="ref37.pdf" target="_blank">PDF</a>)
*/
/*!
Implements the algorithm in:\n
van Gijzen&  Sonneveld (2011)\n
Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits Bi-orthogonality Properties.\n
[ACM ToMS (38) 1-19](http://ta.twi.tudelft.nl/nw/users/gijzen/idrs_toms.pdf).
ported from Matlab code published on the website.
*/

template<typename T, typename U = T> class idrs : public solver<T,U> {
	const T	mp = 1e3*std::numeric_limits<T>::epsilon();	// number close to machine precision:
	const T	angle{0.7};			//!< angle limit for correcting omega
	const bool smooth{false};		//!< switches smoothing on
	const bool replace{false}; 		//!< switches residual replacement on
	const unsigned int m;			//!< subspace dimension

	vecD<U>	orthogonalize(const std::vector<vecD<U> >& V, vecD<U>& v) const		//! orthogonalizes vector v against basis V
		{ unsigned int n = ITOU(V.size()); vecD<U> s(n); s = U(0);		// return something useful if n == 0
		  for (unsigned int i = 0; i < n; i++) s[i] = dot(V[i],v);
		  for (unsigned int i = 0; i < n; i++) v -= V[i]*s[i];
		  return s; }
	void	randomInit(vecD<U>& v) const
		{ random(v); T n = norm(v); v /= U(n); }
	T	omega(const vecD<T>& t, const vecD<T>& s, const T angle) const
		{ const T ns = norm(s), nt = norm(t), ts = dot(t,s);
		  const T rho = std::abs(ts/(nt*ns));
		  T om = ts/(nt*nt); if (rho < angle) om *= angle/rho;
		  return om; }
public:
	idrs(const matS<U>& _A, const precS<T,U>& _M, const unsigned int _m = 4,
		const T _eps = T(1e-6))							//! constructs IDR solver for matrix A and preconditioner M.
		 : solver<T,U>(_A,_M,_eps), m(_m) { }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false)			//! solves A x = b up to tolerance eps.
		{ vecD<U> x(b.N), f(m), xs, rs; x = U(0); 
		  const T bn = norm(b); if (bn < solver<T,U>::eps) return x;
		  bool update = false; T om = T(1); const underscore _;
		  std::vector<vecD<U>> P, G, H; matD<T> M(m,m); M.id();
		  vecD<U> v(b.N); randomInit(v); P.push_back(v);			// init P
		  for (unsigned int i = 1; i < m; i++) { randomInit(v); 
			orthogonalize(P,v); v /= U(norm(v)); P.push_back(v); };
		  for (unsigned int i = 0; i < m; i++) { G.push_back(x); H.push_back(x); };
		  vecD<T> r = b-solver<T,U>::mv(x); T rn = norm(r);
		  if (smooth) { xs = x; rs = r; };
		  for (unsigned int it = 0; it < solver<T,U>::maxit; it++) {
			if (rn/bn < solver<T,U>::eps) return x;
			for (unsigned int k = 0; k < m; k++) f[k] = dot(r,P[k]);	// new RHS for small system
			for (unsigned int k = 0; k < m; k++) { v = r;
				const vecD<T> c = ::solve(M(_(k,m),_(k,m)),f(_(k,m)));	// solve small system
				for (unsigned int j = k; j < m; j++) v -= G[j]*c(j-k);	// make v orthogonal to P
				v = om*solver<T,U>::prec(v);				// preconditioning
				for (unsigned int j = k; j < m; j++) v += H[j]*c(j-k);	// make v orthogonal to P
				H[k] = v; G[k] = solver<T,U>::mv(v);			// compute new U(:,k) and G(:,k), G(:,k) is in space G_j
				for (unsigned int j = 0; j < k; j++) {			// bi-orthogonalise the new basis vectors
					const T alpha = dot(P[j],G[k])/M(j,j);
					G[k] -= alpha*G[j]; H[k] -= alpha*H[j]; };
				for (unsigned int j = k; j < m; j++)			// new column of M = P'*G  (first k-1 entries are zero)
					M(j,k) = dot(G[k],P[j]);
      				if (M(k,k) == T(0)) throw rtException("idrs: zero diagonal");
				const T bt = f(k)/M(k,k); x += bt*H[k];			// make r orthogonal to q_i, i = 1..k 
				r -= bt*G[k]; rn = norm(r); 
				if (replace && rn > bn*solver<T,U>::eps/mp) update = true;
				if (smooth) { const vecD<T> t = rs-r;
					const T gm = dot(t,rs)/dot(t,t);
					rs -= gm*t; xs -= gm*(xs-x); rn = norm(rs); };
				if (rn/bn < solver<T,U>::eps) return smooth? xs: x;
				for (unsigned int j = k+1; j < m; j++) f[j] -= bt*M(j,k); }; // new f = P'*r (first k components are zero)
			v = solver<T,U>::prec(r);					// preconditioning
			const vecD<U> t = solver<T,U>::mv(v);				// computation of a new omega
			om = omega(t,r,angle); if (om == T(0)) throw rtException("idrs: omega breakdown");
     			r -= om*t; x += om*v; 
			if (replace && rn > bn*solver<T,U>::eps/mp) update = true;
			if (update && rn < bn) { r = b-solver<T,U>::mv(x); update = false; };
			rn = norm(r);
			if (verbose) { printf("[%4d] %6.4e\r", it, rn/bn); fflush(stdout); };
				if (rn/bn < solver<T,U>::eps) return smooth? xs: x;
			if (smooth) { const vecD<T> t = rs-r;
				const T gm = dot(t,rs)/dot(t,t);
				rs -= gm*t; xs -= gm*(xs-x); rn = norm(rs); } }
		  printf("idrs: no convergence.\n"); return smooth? xs: x; }
};

//! A direct solver for semi-positive definite sparse matrices

template <typename T, typename U = T> class LDL {
#define allRows(r) (unsigned int r = 0; r < L.M; r++)
#define allCols(t,r) (unsigned int t = L.ri[r]; t < L.ri[r+1]; t++)
#define allACols(t,r) (unsigned int t = A.ri[r]; t < A.ri[r+1]; t++)
	matS<U> L;				//!< lower triangular matrix
	vecD<U>	D;				//!< diagonal vector

	void	decompose(const matS<U>& A, const uvecD& par)				//! performs the numeric decomposition phase
                { unsigned int m = L.M, len; uvecD flag(m), pat(m), nc(m);		// note: uses A's upper triangle only
		  vecD<U> y(L.N); y = U(0); nc = 0;
                  for allRows(r) { unsigned int top = m; flag[r] = r;
                        for allACols(t,r) { unsigned int c = A.ci[t]; if (c > r) continue;
                                y[c] += A.x[t];
                                for (len = 0; flag[c] != r; c = par(c)) { pat[len++] = c; flag[c] = r; }
                                while (len > 0) pat[--top] = pat[--len]; };
                        D[r] = y[r]; y[r] = U(0);
                        for ( ; top < m; top++) {
                                unsigned int c = pat[top], t = L.ri[c]; U yc = y[c]; y[c] = U(0);
                                for ( ; t < L.ri[c]+nc[c]; t++) y[L.ci[t]] -= L.x[t]*yc;
                                U v = yc/D[c]; D[r] -= v*yc; L.ci[t] = r; L.x[t] = v; nc[c]++; };
                        if (D[r] == U(0)) throw rtException("LDL: zero diagonal"); } }
	void	assign(const LDL<T,U>& b)
		{ L = b.L, D = b.D; }
public:
	LDL()	: L(), D() { }
	LDL(const matS<U>& A)								//! constructs direct LDL solver for matrix A.
		: L(), D(A.M)
		{ // allocate LL^T and perform the symbolic decomposition phase - uses A's upper triangle only
		  unsigned int n = D.N; uvecD f(n), nc(n), par(A.M);
		  for (unsigned int r = 0; r < n; r++) {
			f[r] = r; nc[r] = 0; par[r] = ~0u;
			for allACols(t,r) { unsigned int c = A.ci[t]; if (c >= r) continue;
				for ( ; f[c] != r; c = par[c]) {
					if (par[c] == ~0u) par[c] = r;
					nc[c]++; f[c] = r; } } };
		  const unsigned int nz = nc.sum(); L.resize(n,n,nz);
		  for allRows(r) L.ri[r+1] = L.ri[r]+nc[r];
		  decompose(A,par); }
	LDL(const LDL<T,U>& b)								//! copies from matrix b.
		 : L(), D() { assign(b); }
	vecD<U>	solve(const vecD<U>& b, const bool verbose = false) const		//! solves A x = b.
		{ (void)verbose; vecD<U> x(b); clearFpException();
		  for allRows(r) { for allCols(t,r) x[L.ci[t]] -= L.x[t]*x[r]; };
		  for allRows(r) x[r] /= D(r);
		  for (unsigned int r = L.M-1; r != ~0u; r--) {
			for allCols(t,r) x[r] -= L.x[t]*x[L.ci[t]]; };
		  checkFpException(); return x; }
	LDL<T,U>& operator=(const LDL<T,U>& b)						//! assigns from matrix b.
		{ if (this != &b) assign(b); return *this; }
	LDL<T,U>& operator=(LDL<T,U>&& b)						//! assigns from matrix b.
		{ assert(this != &b); L = std::move(b.L); D = std::move(b.D); return *this; }
#undef allCols
#undef allRows
#undef allARows
};

#endif		// solver.h

