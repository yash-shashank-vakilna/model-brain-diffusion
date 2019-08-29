#ifndef GRADOPTIMIZER_H
#define GRADOPTIMIZER_H

/*
 *
 * gradOptimizer.h: gradient-based optimization
 * BRIAN Software Package Version 3.0
 *
 * $Id: gradOptimizer.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (28/08/15): initial release
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions for optimization classes.
*/


//! Implements the base class for problems using gradient-based optimizers.

template <typename T>
class cfo {
public:
	cfo()
		{ }
	virtual ~cfo()
		{ }
	virtual T	func(const vecD<T>& p) = 0;	//!< returns function value at parameter vector p.
	virtual vecD<T> grad(const vecD<T>& p) = 0;	//!< returns gradient at parameter vector p.
	virtual matS<T> hess(const vecD<T>& p) = 0;	//!< returns Hessian at parameter vector p.
};

//! Implements an optimizer based on scaled conjugate gradients.

template <typename T>
T CG(cfo<T>& obj, vecD<T>& p, const T gtol = T(1e-6), const unsigned int maxit = 500u, const bool verbose = false)
//! Given cost function object cfo and initial parameter vector p, returns optimized function value and parameter set (in p).
{	T f = obj.func(p); vecD<T> r = obj.grad(p)*T(-1.0), d = r;			// function value, negative gradient and search direction
	bool success = true; vecD<T> s;           					// used as approximation to H*p in loop below
	T sigma = T(1.0e-2), lambda = T(0.1), lbar = T(0), delta = T(0);		// update for lambda
	for (unsigned int it = 0; it < maxit; it++) { const T d2 = dot(d,d);
		if (success == true) {                                  		// if last step led to reduction of cost-function
			const T sigma_k = sigma/std::sqrt(d2);                      	// normalized step-length when estimating H*p
			s = (obj.grad(p+sigma_k*d)+r)/sigma_k;				// approximation to H*p
			delta = dot(d,s); };						// approximation to p'*H*p
		s += (lambda-lbar)*d; delta += (lambda-lbar)*d2;			// if <0 then H+(l-lb)*I not positive definite
		if (delta <= T(0)) {							// if it H is not positive definite
      			s += (lambda-T(2.0)*(delta/d2))*d;            			// make H more diagonal dominant to ensure pos def
      			lbar = T(2.0)*(lambda-delta/d2);
      			delta = lambda*d2-delta; lambda = lbar; }
		const T mu = dot(d,r), alpha = mu/delta, nf = obj.func(p+alpha*d);	// value of cost-function at attempted new point
		const T Delta = T(2.0)*delta*(f-nf)/SQR(mu);				// >0 means attempted step reduced cost-function
		if (Delta >= T(0)) { f = nf; p += alpha*d; lbar = T(0); success = true;
			if (it % p.N == 0) { r = obj.grad(p)*T(-1.0); d = r; }
			else { vecD<T> oldr = r; r = obj.grad(p)*T(-1.0);
				T beta = (dot(r,r)-dot(oldr,r))/mu; d = r+beta*d; }
			if (Delta > T(0.75)) lambda *= T(0.5); }
		else { lbar = lambda; success = false; }
		if (Delta < T(0.25)) lambda *= T(4.0);
		T test = T(0);
		for (unsigned int i = 0; i < p.N; i++)
			test = std::max(test,std::abs(r(i))*std::max(std::abs(p(i)),T(1.0)));
		test /= std::max(f,T(1.0));
		if (verbose) { printf("[%u] %e %e\r", it, nf, test); fflush(stdout); };
		if (test < gtol) return f; };
	return 0; // throw rtException("CG did not converge");
}

//! Implements an optimizer based on the Levenberg-Marquardt method.

template <typename T>
T LM(cfo<T>& obj, vecD<T>& p, const T tol = T(1e-6), const unsigned int maxit = 500u, const bool verbose = false)
//! Given cost function object cfo and initial parameter vector p, returns optimized function value and parameter set (in p).
{	T f = obj.func(p); bool success = true;						// calculate initial values
	T lambda = T(0.1), olambda = T(0); vecD<T> g; matS<T> H;
	for (unsigned int it = 0; it < maxit; it++) {
		if (success) { g = obj.grad(p); H = obj.hess(p); }
		const T lf = (T(1.0)+lambda)/(T(1.0)+olambda);
		for (unsigned int i = 0; i < p.N; i++) H(i,i) *= lf;
		const precAinv<T> pre(H,T(0.01)); bicgstab<T> sv(H,pre); 
		vecD<T> step; bool solveFailed = false;
		try { step = sv.solve(step); }
		catch(...) { solveFailed = true; };
		const T nf = obj.func(p-step);
		if (T(2.0)*std::abs(f-nf) <
			tol*(std::abs(f)+std::abs(nf)+std::numeric_limits<T>::epsilon())) return nf;
		success = nf <= f;
		if (verbose) { printf("[%2d] %e %e %e %e %u\r", it, nf, f, f-nf, lambda, success); fflush(stdout); };
		if (success && solveFailed == false) {
			olambda = 0; p -= step; lambda *= T(0.1); f = nf; }
		else { olambda = lambda; lambda *= T(10.0); } };
	return 0; // throw rtException("LM did not converge");
}
#endif

