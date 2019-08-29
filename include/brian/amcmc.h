#ifndef AMCMC_H
#define AMCMC_H

/*
 *
 * amcmc.h: adaptive Markov chain Monte Carlo sampler
 * BRIAN Software Package Version 3.0
 *
 * $Id: amcmc.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (24/05/10): first implementation
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 * For detailed information, refer to (\cite Rosenthal07, <a href="ref1.pdf" target="_blank">PDF</a>)
 */

/*! \file
    \brief Adaptive Markov chain Monte Carlo sampler.

    \details For detailed information, refer to (\cite Rosenthal07, <a href="ref1.pdf" target="_blank">PDF</a>).
*/

#define allDim(d)	(unsigned int d = 0; d < D; d++)

//! Implements an adaptive Markov chain Monte Carlo sampler

template<typename T> class amcmc {
	const float bfrac = T(0.2);
	const unsigned int D;			//!< problem dimensionality
	const unsigned int N;			//!< number of blocks
	const unsigned int L;			//!< block length
	const unsigned int B;			//!< burn in period
	const unsigned int S;			//!< subsampling factor
	const T ar;				//!< acceptance ratio for variance scaling
	vecD<T>	xo;				//!< old location
	vecD<T>	xn;				//!< next location
	vecD<T>	ls;				//!< move variance scaling
	std::vector<vecD<T>> x;			//!< jumps
	uniformDist<T> ud;

	vecD<T>	acf(const unsigned int d, const unsigned int _M) const;			// autocorrelation function for variable d
	vecD<T>	quantiles(const unsigned int d) const;					// median of variable d
        virtual T density(const vecD<T>& xn) = 0;					// sample density, x_i \in [0,1]
        virtual T functional(const vecD<T>& xn) = 0;					// sample functional
public:
	amcmc(const unsigned int _D, const unsigned int _N, const unsigned int _L,
		const unsigned int _S = 1, const T iv = 1.0)				//! allocates a MCMC sampler.
		/*! Parameters: D the dimensionality of the problem, N the number of blocks,
		    L the block length, S the subsampling factor, and iv the initial location.
		    For detailed information, refer to (\cite Rosenthal07, <a href="ref1.pdf" target="_blank">PDF</a>).
		 */
		 : D(_D), N(_N), L(_L), B(static_cast<unsigned int>(N*L*bfrac)), S(_S), 
		ar(T(L*std::max(0.234,std::min(0.44,2.38/std::sqrt(_D))))), xo(D), xn(D), ls(D), x((N*L-B)/S), ud()
		{ xo = iv; xn = xo; ls = 0.0; }
	virtual ~amcmc() { }
	std::vector<vecD<T>> acf(const unsigned int M) const				//! returns ACF for all dimensions.
		{ std::vector<vecD<T>> a(D);
		  for allDim(d) { a[d] = acf(d,M); }; return a; }
	std::vector<vecD<T>> quantiles() const						//! returns quantiles for all dimensions.
		{ std::vector<vecD<T>> a(D);
		  for allDim(d) { a[d] = quantiles(d); }; return a; }
	vecD<T> median() const								//! returns median for all dimensions.
		{ vecD<T> a(D);
		  for allDim(d) { vecD<T> m = quantiles(d); a[d] = m(3); }; return a; }
	T	work(bool verbose = false);
};

template <typename T> 
T amcmc<T>::work(bool verbose)
//! sample over N blocks of length L using variance scaling.
{	unsigned int it = 0, ns = 0, ss = 0; bool fr = false;
	T fv = T(0.0), fe = T(0.0), fs = T(0.0), po = density(xo); ls = 0;		// re-init variance scaling
	for (unsigned int n = 1; n <= N; n++) { ivecD ac(D); ac = 0;			// for all blocks
		for (unsigned int l = 0; l < L; l++) { it++; 				// for all steps in this block
			for allDim(d) {							// for all parameters of this model
				xn[d] = xo[d]+T(std::exp(ls[d]))*normalDist<T>::sample();	// propose move
				T pn = density(xn);					// compute posterior probability
				if (ud() < std::exp(pn-po)) {				// accept move
					xo[d] = xn[d]; ac[d]++; po = pn; fr = false; }
				else xn[d] = xo[d]; }					// reject move
			if (it > B && ++ss == S) { ss = 0;				// after burnin, evaluate proposal
				if (fr == false) { fv = functional(xo); fr = true; };
				if (ns < x.size()) { x[ns] = xo; ns++; fs += fv; } } };	// keep result
        	fe = ns? fs/ns: 0.0;							// at the end of a block: average std::log likelihood
		T a = T(std::min(1.0/std::sqrt(n), 0.01));					// adapt variances
		for allDim(d) ls[d] += (ac[d] > ar)? a: -a;				// change std::log sigma based on acceptance
		if (verbose) {
			printf("[%4d] sq %e, f %e ", n, norm(xo), fe);
			printf("ls "); print(ls); printf("\n"); fflush(stdout); } };
	if (verbose) printf("\n");
	return fe;
}

template <typename T> 
vecD<T> amcmc<T>::acf(const unsigned int d, const unsigned int _m) const
//! compute autocorrelation function for variable d.
{	unsigned int n = x.size(), m = std::min(_m,n); vecD<T> tr(n), a(m); T s = 0;		// allocate local variables
	for (unsigned int t = 0; t < n; t++) { T xt = x[t](d); s += xt; tr[t] = xt; };	// copy to vector
	tr -= s/n;									// subtract mean
	for (unsigned int i = 0; i < m; i++) { s = 0;
		for (unsigned int t = 0; t < n-m; t++) s += tr[t]*tr[t+i];
		a[i] = s; };
	return a(0)? a/a(0): T(1.0);
}

template <typename T> 
vecD<T> amcmc<T>::quantiles(const unsigned int d) const
//! compute median of variable d.
{	unsigned int n = ITOU(x.size()); vecD<T> tr(n);
	for (unsigned int t = 0; t < n; t++) tr[t] = x[t](d);				// copy to vector
	return tr.quantiles();
}

#undef allDim

#endif

