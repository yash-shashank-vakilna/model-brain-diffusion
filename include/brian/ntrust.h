#ifndef NTRUST_H
#define NTRUST_H

/*
 *
 * ntrust.h: Newton trust-region solver
 * BRIAN Software Package Version 3.0
 *
 * $Id: ntrust.h 425 2016-10-24 02:16:19Z frithjof $
 *
 * 0.10 (06/03/16): first implementation
 * v406 (28/09/16): bumped to version 3.0
 *
 * based on MATLAB code by C.T. Kelley, 1997
 *
 * see:
 * Steihaug T (1983)
 * The conjugate gradient method and trust regions in large scale optimization.
 * SIAM J Numer Anal 20, 626-637.
 *
 */

/*! \file
    \brief Newton trust-region solver

    \details based on MATLAB code by C.T. Kelley, 1997.
*/

//! Implements a context for optimizing nonlinear problems using a Newton trust region method.

template<typename T> class ntrust {
	const T	tol;				//!< final tolerance of solution.
	const unsigned int maxit;		//!< maximum number of iterations.
	const T	eps;				//!< a small number.

	virtual T funcAt(const vecD<T>& x) const = 0;					//! evaluates function at x (nb. implement in subclass).
	virtual vecD<T>	gradAt(const vecD<T>& x) const					//! evaluates gradient at x.
		{ const T y = funcAt(x); const unsigned int n = x.N;
		  vecD<T> g(n); g = T(0); for (unsigned int i = 0; i < n; i++) {
			vecD<T> w(n); w = T(0); w[i] = eps;
			g[i] = (funcAt(x+w)-y)/eps; };
		  return g; }
	virtual vecD<T>	precAt(const vecD<T>& x) const					//! applies preconditioner to x.
		{ return x; }
	virtual vecD<T>	dirDeriv(const vecD<T>& x, const vecD<T>& s, const vecD<T>& g) const //! returns a finite difference directional derivative.
		{ const T ns = norm(s); 
		  if (ns > T(0)) { const T d = eps/ns;
			return (gradAt(x+d*s)-g)/d; }
		  else { vecD<T> z(x.N); z = T(0); return z; } }
	matD<T>	HessianAt(const vecD<T>& x, const vecD<T>& g)				//! returns a forward difference Hessian f''(x). Required for dogleg().
		{ const unsigned int n = x.N; matD<T> H(n,n); H = T(0);
		  for (unsigned int i = 0; i < n; i++) {
			vecD<T> s(n); s = T(0); s[i] = T(1.0);
			H.setCol(i, dirDeriv(x,s,g)); };
		  return symmetrize(H); }
	T	solveEq(const vecD<T>& x, const vecD<T>& p, const T d) const		//! solves linear equation (x,p,d).
		{ const T a = dot(p,p), b = T(2.0)*dot(x,p), c = dot(x,x)-SQR(d);
		  return (-b+std::sqrt(SQR(b)-T(4.0)*a*c))/(T(2.0)*a); }
	bool	getPoint(vecD<T>& xn, const matD<T>& H, const vecD<T>& g, const vecD<T>& xc,
			const T rc) const;
	vecD<T>	solveCG(std::vector<vecD<T>>& dir, const vecD<T>& x, const vecD<T>& g, const T d,
			const unsigned int nit = 20) const;
	vecD<T>	getStep(std::vector<vecD<T>>& dir, const T r) const;
public:
	ntrust(const T _tol = T(1e-6), const unsigned int _maxit = 100,
		const T _eps = T(1e-6))							//! Allocates a context for nonlinear optimization.
		: tol(_tol), maxit(_maxit), eps(_eps)
		{ }
	virtual ~ntrust()
		{ }
	vecD<T>	dogleg(const vecD<T>& x0, const bool verbose = false);
	vecD<T>	cg(const vecD<T>& x0, const bool verbose = false);
};

template <typename T>
bool ntrust<T>::getPoint(vecD<T>& xn, const matD<T>& H, const vecD<T>& g, const vecD<T>& xc,
	const T rc) const
//! returns new point xn based on Hessian H and gradient g at xc and trust region radius rc.
{	T mu = dot(g,H*g), sigma = T(0), gn = norm(g); bool bnd = false;
	if (mu > T(0)) { sigma = dot(g,g)/mu; if (sigma*gn > rc) sigma = rc/gn; }	// compute Cauchy point
	else { bnd = true; sigma = rc/gn; }						// if steepest descent is a direction of negative curvature, leap to TR boundary
	const vecD<T> p = xc-sigma*g;
	if (bnd) xn = p;								// if p is on the TR boundary, that's the trial point.
	else {	const vecD<T> ds = g*T(-1.0), dn = inv(H)*ds, xi = xc+dn;		// else p is in the interior and the steepest descent is a direction of positive curvature
		if (dot(ds,dn) <= T(0)) xn = p;						// if the Newton direction goes uphill, revert to p.
		else if (norm(dn) <= rc) xn = xi;					// else if the Newton point is inside, accept it.
		else {	const vecD<T> d1 = sigma*g; const T t = solveEq(d1,d1+dn,rc);	// else find the intersection of the dog leg path with TR boundary.
			xn = p+t*(xi-p); bnd = true; } };
	return bnd;
}

template <typename T>
vecD<T> ntrust<T>::dogleg(const vecD<T>& x0, const bool verbose)
//! implements a dog leg trust region, Newton model, dense algorithm. 
{	T y = funcAt(x0); vecD<T> x = x0, xo = x, xn = x, g = gradAt(x); 		// current, old and new position
	T r = std::min(norm(g),T(10.0)), ro = r, rn = r;				// current, old and new trust radius
	unsigned int step = 1; matD<T> H;
	for (unsigned int it = 0; it < maxit; it++) { if (norm(g) <= tol) break;
		if (step == 1) H = HessianAt(x,g);					// we're good, so compute new Hessian
		const bool bnd = getPoint(xn,H,g,x,r); const vecD<T> dx = xn-x;		// now adjust the TR radius using the trial point
		const T t = dot(g,dx)+T(0.5)*dot(dx,H*dx), q = (funcAt(xn)-y)/t;
		if (q < T(0.25)) { xn = x; r = norm(dx)*T(0.5);				// reduce radius
				step = step == 3? 4: 2; }
		else if (q > T(0.75) && bnd) { r *= T(2.0); step = 3; }			// expand radius
		else step = 1;
		rn = r;									// set x and r for next round
		if (step == 3) { xo = xn; ro = r; }
		else if (step == 4) { xn = xo; rn = ro; step = 1; };
		x = xn; r = rn;
		if (verbose) { printf("%e %e %e %e\n", y, x(0), x(1), norm(g)); };
		if (step == 1) { y = funcAt(x); g = gradAt(x); } };			// we're good, so compute new function value and gradient
	return x;
}

template <typename T>
vecD<T> ntrust<T>::solveCG(std::vector<vecD<T>>& dir, const vecD<T>& x, const vecD<T>& g, const T rc,
	const unsigned int nit) const
//! solves the trust region problem with preconditioned conjugate-gradient.
{	vecD<T> xn(x.N); xn = T(0); vecD<T> b = g*T(-1.0), r = b-dirDeriv(x,xn,g);	// set initial guess and residuals
	vecD<T> z = precAt(r), p = z; T rho = dot(z,r), rn = norm(r);
	const T cv = eps*norm(b), hdel = rc*(T(1.0)-eps); T rhoold = T(0);
	for (unsigned int it = 0; it < nit; it++) {
		if (rn < cv || norm(xn) > hdel) break;					// check for convergence
		p = it? z+(rho/rhoold)*p: z; vecD<T> w = dirDeriv(x,p,g);
		T alpha = dot(p,w);
		if (alpha <= T(0)) alpha = solveEq(xn,p,rc);				// alpha < 0: head to the TR boundary 
		else {	alpha = rho/alpha;						// else set new alpha
			if (norm(xn+alpha*p) > rc) alpha = solveEq(xn,p,rc); };
		xn += alpha*p; dir.push_back(alpha*p); r -= alpha*w;			// get new gradient and residuals
		rn = norm(r); rhoold = rho; z = precAt(r); rho = dot(z,r); };
	return xn;
}

template <typename T>
vecD<T> ntrust<T>::getStep(std::vector<vecD<T>>& dir, const T r) const
//! finds the point of intersetion of the trust region boundary and the PL path.
{	vecD<T> s = dir[0]; bool inside = true; unsigned int k = 0;
	if (norm(s) > r || dir.size() == 1) { dir.resize(1); return s*(r/norm(s)); };
	for (unsigned int i = 1; i < dir.size(); i++) {
		if (norm(s+dir[i]) > r && inside) {
			k = i; const vecD<T> p = dir[i];
			s += solveEq(s,p,r)*p;  inside = false; }
		else s += dir[i]; };
	dir.resize(k+1); return s;
}

template <typename T>
vecD<T> ntrust<T>::cg(const vecD<T>& x0, const bool verbose)
//! implements a Steihaug Newton-CG-Trust region algorithm.
{	T r = norm(x0), y = funcAt(x0); vecD<T> g = gradAt(x0), x = x0; 
	for (unsigned int it = 0; it < maxit; it++) { // if (norm(g) <= tol) break; 
		std::vector<vecD<T>> dir; const vecD<T> s = solveCG(dir,x,g,r);  	// get a new step
		vecD<T> xn = x+s, w = dirDeriv(x,s,g); T yn = funcAt(xn); 		// now adjust the TR radius using the trial point
		T pred = dot(g,s)+T(0.5)*dot(w,s), q = (yn-y)/pred;
		if (q > T(0.25)) { x += s; y = funcAt(x); g = gradAt(x);
			if (q < T(0.25)) r *= T(0.5); }					// reduce radius
			else if (q > T(0.75) && norm(s) > r-eps) r *= T(2.0);		// expand radius
		else {	for (unsigned int t = 0; q <= T(0.25); t++) { if (t >= 30) return x;
				r = T(0.5)*std::min(r,norm(s));
				const vecD<T> s = getStep(dir,r); xn = x+s; yn = funcAt(xn); 
				if (yn < y) { w = dirDeriv(x,s,g);
					pred = dot(g,s)+T(0.5)*dot(w,s); };
				q = (yn-y)/pred; };
			x = xn; y = funcAt(x); g = gradAt(x); };
		if (verbose) { printf("[%3d] %e %e\n", it, y, norm(g)); } };
	return x;
}

#endif

