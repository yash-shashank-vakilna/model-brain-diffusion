#ifndef OPTIMIZER_H
#define OPTIMIZER_H

/*
 *
 * optimizer.h: linear & nonlinear minimization
 * BRIAN Software Package Version 2.5
 *
 * $Id: optimizer.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.30 (18/11/13): initial release
 * 0.40 (16/12/13): documented
 * 0.50 (22/08/15): gradient-based optimizers added
 * 0.60 (16/01/16): simplex: use start position in simplex
 *
 */

/*! \file
    \brief Definitions for optimization classes.
*/

//! Implements minimization in one dimension

template<class T> class linemin {
	const T	eps;					//!< residual tolerance
	virtual T goal(const T p) = 0;							//! interface to cost function.
public:
	linemin<T>()									//! allocates an empty one-dimensional optimizer.
		 : eps(T(std::sqrt(std::numeric_limits<T>::epsilon()))) { }
	virtual ~linemin() { }
	T	optimize(T& pos, const T min, const T max, const T th);
};

template <class T>
T linemin<T>::optimize(T& x, const T min, const T max, const T t)
//! optimizes between min and max, returns minimum in pos.
{	// seeks a local minimum of a function F(X) in an interval [min,max].
	// see Brent (2002) Algorithms for Minimization Without Derivatives, Dover.
	T c = T(0.5*(3.0-std::sqrt(5.0))), sa = min, sb = max; x = sa+c*(max-min);		// c is the square of the inverse of the golden ratio
	T w = x, v = w, e = 0.0, fx = goal(x), fw = fx, fv = fw, d = 0, u = 0;
	while (true) {	T m = T(0.5*(sa+sb)), tol = eps*std::abs(x)+t, t2 = T(2.0*tol);
		if (std::abs(x-m) <= T(t2-0.5*(sb-sa))) break;				// check the stopping criterion
		T r = 0.0, q = r,p = q;							// fit a parabola
		if (tol < std::abs(e)) {
			r = (x-w)*(fx-fv); q = (x-v)*(fx-fw);
			p = (x-v)*q-(x-w)*r; q = T(2.0*(q-r)); if (0.0 < q) p = -p;
			q = std::abs(q); r = e; e = d; }
		if (std::abs(p) < std::abs(0.5*q*r) && q*(sa-x) < p && p < q*(sb-x)) {
			d = p/q; u = x+d;						// take the parabolic interpolation step
			if (u-sa < t2 || sb-u < t2) d = x < m? tol: -tol; }		// f must not be evaluated too close to a or b
		else { e = x < m? sb-x: sa-x; d = c*e; };				// a golden-section step
		if (tol <= std::abs(d)) u = x+d;
		else if (0.0 < d) u = x+tol; else u = x-tol;				// f must not be evaluated too close to x
		T fu = goal(u);
		if (fu <= fx) {	if (u < x) sb = x; else sa = x;				// update a, b, v, w, and x
			v = w; fv = fw; w = x; fw = fx; x = u; fx = fu; }
		else {	if (u < x) sa = u; else sb = u;
			if (fu <= fw || almostEqual(w,x)) { v = w; fv = fw; w = u; fw = fu; }
			else if (fu <= fv || almostEqual(v,x) || almostEqual(v,w)) { v = u; fv = fu; } } };
	return fx;
}

//! Implements the Levenberg-Marquardt algorithm for nonlinear optimization

template <class T> class LevenbergMarquardt {
	static const unsigned int maxit = 100;
	const T	eps = T(1e-8);
	const T	cv = T(1e-2);
	const T	tau = T(1e-3);
	virtual vecD<T> adder(const vecD<T>& x, const vecD<T>& h) const			// default is Euclidean update, subclass may override.
		{ return x+h; }
	virtual matD<T> estimateJacobian(const vecD<T>& x, const vecD<T>& inc) const	//! estimate the Jacobian of the goal function above
		{ unsigned int m = 0, n = x.N; vecD<T> xm = x; matD<T> ret;
		  for (unsigned int j = 0; j < n; j++) { const T f = T(0.5/inc(j));
			xm[j] = x(j)+inc(j); vecD<T> df = goal(xm);			// create the modified "x" vector:
			xm[j] = x(j)-inc(j); df -= goal(xm); xm[j] = x(j);
			if (j == 0) { m = df.N; ret.resize(m,n); }; ret.setCol(j,f*df); };
		  return ret; }
	virtual vecD<T> goal(const vecD<T>& p) const = 0;				//! interface to cost function.
public:
	LevenbergMarquardt(const T _eps = T(1e-8)) : eps(_eps)
		{ }
	virtual ~LevenbergMarquardt() { }
	vecD<T> optimize(vecD<T>& x, const vecD<T>& inc, const bool verbose = false);
};

template <typename T>
vecD<T> LevenbergMarquardt<T>::optimize(vecD<T>& x, const vecD<T>& inc, const bool verbose)
{	assert(inc.N == x.N);
	matD<T> J = estimateJacobian(x,inc), H = trp(J)*J;				// compute the Jacobian and Hessian of goal()
	vecD<T> f = goal(x), g = trp(J)*f;						// compute the gradient
	T sw = tau*H.maxD(), v = T(2.0), F = dot(f,f);
	for (unsigned int it = 0; norm(g) > eps && it < maxit; it++) {
		matD<T> K = H; K.addToD(sw); vecD<T> h = inv(K)*g*(T(-1.0));		// h = -( H + step I )^-1 * g
		if (verbose) printf("[%u] %e\n", it, norm(g));
//		if (norm(h) < cv*(norm(x)+cv)) break;					// check convergence
		vecD<T>	xn = adder(x,h), fn = goal(xn);					// increment x and update f
		const T Fn = dot(fn,fn), l = (F-Fn)/dot(h,sw*h-g);			// denom = h^t * (step * h - g)
		if (l > 0) { x = xn; f = fn; F = Fn;					// there is an improvement, so accept new x
			J = estimateJacobian(x,inc); H = trp(J)*J; g = trp(J)*f;	// update Jacobian, Hessian, gradient
			sw *= T(std::max(0.33,1.0-pow(2.0*l-1.0,3.0))); v = T(2.0); }	// and step width
		else {	sw *= v; v *= T(2.0); } }					// no improvement
	return f;
}

//! Helper class for simplex algorithm, represents vertex position & function value

template<class T> struct svtx {								// local, represents vertex position & function value
	vecD<T>	pos;					//!< current position
	T	val;					//!< function value
	svtx()										//! allocates an empty svtx record.
		 : pos(), val(T(0)) { }
	bool operator< (const svtx<T>& rhs) const					//! comparison by increasing function value.
		{ return (this->val < rhs.val); }
};

//! Implements non-linear minimization in N dimensions using the simplex algorithm

template<class T> class simplex {
	static const unsigned int maxit = 100;		//!< maximum number of optimization iterations
	const unsigned int np;				//!< number of points in simplex
	T	eps;					//!< residual tolerance
	std::vector<svtx<T>> v;				//!< vertices and function values
	vecD<T> ce;					//!< center
	vecD<T> pos;					//!< current position
	vecD<T> range;					//!< range

	T	eval(const unsigned int i, const T fac);	
	void	center()								//! returns center of simplex.
		{ ce = T(0); for (unsigned int i = 0; i < np; i++) ce += v[i].pos; }
	T	conv() const								//! computes convergence criterion.
		{ return std::abs(v[0].val-v[np-1].val)/(std::abs(v[0].val)+std::abs(v[np-1].val)+eps); }
	virtual T goal(const vecD<T>& p) = 0;						// interface to cost function.
public:
	simplex<T>(const unsigned int dim, T _eps = T(1e-8))				//! allocates a simplex optimizer for dim dimensions.
		 : np(dim+1), eps(_eps), v(np), ce(dim), pos(), range()
		{ vecD<T> t(dim); t = T(0);
		  for (unsigned int i = 0; i < dim; i++) { v[i].pos = t; v[i].pos[i] = T(1.0); };
		  v[dim].pos = t; center(); }
	virtual ~simplex() { }
	T	optimize(vecD<T>& pos, const vecD<T>& r, const bool verbose = false);
};

template <typename T>
T simplex<T>::eval(const unsigned int i, const T fac)
//! evaluate cost function for vertex i and factor fac.
{	T fac1 = (T(1.0)-fac)/T(np); T fac2 = fac1-fac;
	vecD<T> dp = ce*fac1-v[i].pos*fac2, p = pos+range*dp; T val = goal(p);
	if (val < v[i].val) { v[i].pos = dp; v[i].val = val; }; return val;
}
	
template <typename T>
T simplex<T>::optimize(vecD<T>& p, const vecD<T>& r, const bool verbose)
//! optimize starting at p with range r and final tolerance eps, returns minimum in p.
{	pos = p; range = r; assert(pos.N == np-1); assert(range.N == np-1);
	for (unsigned int i = 0; i < np; i++) v[i].val = eval(i, T(1.0));
	for (unsigned int it = 0; it < maxit; it++) {
		std::sort(v.begin(), v.end()); center();				// sort by increasing function value
		if (conv() < eps) break;						// if converged, return best point
		T val = eval(np-1, T(-1.0));						// reflect worst point
		if (val < v[0].val) val = eval(np-1, T(2.0));				// try std::expanding the simplex in this direction
		else if (val >= v[np-1].val) {						// else if worse
			T sv = v[np-1].val; val = eval(np-1, T(0.5));			// try contracting the simplex
			if (val <= sv) continue;					// if still worse
			for (unsigned int i = 0; i < np; i++) {				// contract the simplex around the best point
				v[i].pos = (v[i].pos+v[0].pos)*T(0.5);
				vecD<T> q = pos+range*v[i].pos; v[i].val = goal(q); } };
		if (verbose) { printf("[%2d] %f\r", it, v[0].val); fflush(stdout); } };
	if (verbose) printf("\n");
	p += range*v[0].pos; return v[0].val;
}

//! Implements non-linear minimization in N dimensions using Powell's algorithm

template<class T> class powell {
	const unsigned int maxit = 100;			//!< maximum number of optimization iterations
	const unsigned int dims;			//!< number of dimensions
	T		eps;				//!< residual tolerance
	unsigned int parts;				//!< number of partitions
	vecD<T>		vmin;				//!< lower bounds of search space
	vecD<T>		vmax;				//!< upper bounds of search space
	vecD<T>		start;				//!< start point
	T		best;				//!< best function value

	T	eval(const vecD<T>& pos, const T t, const vecD<T>& dir)			//! evaluate cost function in direction dir, from position pos.
		{ vecD<T> p = pos+t*dir; return goal(p); }
	void	minimizeOnInterval(vecD<T>& pos, const vecD<T>& dir);
	virtual T goal(const vecD<T>& p) = 0;						// interface to cost function.
public:
	powell<T>(const unsigned int n, T _eps = T(1e-8))				//! allocates a Powell optimizer for n dimensions.
		 : dims(n), eps(_eps), parts(0), vmin(), vmax(), start(), best(T(0)) { }
	virtual ~powell() { }
	T	optimize(vecD<T>& pos, const vecD<T>& min, const vecD<T>& max, 
			const unsigned int p = 16, const bool verbose = false);
};

template <typename T>
void powell<T>::minimizeOnInterval(vecD<T>& pos, const vecD<T>& dir)
//! find minimum from position pos in direction dir.
{	T tmax = std::numeric_limits<T>::max(), tmin = -tmax; 
	for (unsigned int i = 0; i < dims; i++) {
		if (pos[i] < vmin[i]) pos[i] = vmin[i];
		if (pos[i] > vmax[i]) pos[i] = vmax[i];
		const T b0 = vmin[i]-pos[i]; T b1 = vmax[i]-pos[i];
		if (dir(i) > eps) { tmin = std::max(tmin,b0/dir(i)); tmax = std::min(tmax,b1/dir(i)); }
		else if (dir(i) < -1.0*eps) { tmax = std::min(tmax,b0/dir(i)); tmin = std::max(tmin,b1/dir(i)); } };
	T tran = tmax-tmin, t = 0;
	T fmin = std::numeric_limits<T>::max(); unsigned int imin = 0;
	for (unsigned int i = 0; i <= parts; i++) {					// bracket a minimum on the interval [pos+tmin*dir,pos+tmax*dir]
		const T f = eval(pos, tmin+T(i)*tran/T(parts), dir);
		if (f < fmin) { imin = i; fmin = f; } };
	if (imin == 0) { t = tmin; imin = 1; }
	else if (imin == parts) { t = tmax; imin = parts-1; }
	else { t = tmin+imin*tran/parts; };
	T f = fmin; T t0 = tmin+(imin-1)*tran/parts; T f0 = eval(pos, t0, dir);
	T t1 = tmin+(imin+1)*tran/parts; T f1 = eval(pos, t1, dir);
	for (unsigned int i = 0; i <= parts; i++) {					// use inverse parabolic interpolation to find the minimum	
		if (std::abs(t1-t0) <= eps*std::abs(t)+eps) break;			// test for convergence (do not change these parameters)
		const T dt0 = t0-t, dt1 = t1-t, df0 = f0-f, df1 = f1-f;			// compute vertex of interpolating parabola
		const T temp0 = dt0*df1, temp1 = dt1*df0, delta = temp1-temp0;
		if (std::abs(delta) < eps) break;
		const T tmid = t+T(0.5)*(dt1*temp1-dt0*temp0)/delta;			// update bracket
		const T fmid = eval(pos, tmid, dir);
		if (tmid < t) { if (fmid <= f) { t1 = t; f1 = f; t = tmid; f = fmid; }
		else { if (isSmall<T>(t0) && f0 <= fmin) { t = t0; break; } else { t0 = tmid; f0 = fmid; } } }
		else if (tmid > t) { if (fmid <= f) { t0 = t; f0 = f; t = tmid; f = fmid; }
		else { if (isSmall<T>(t1-1.0) && f1 <= fmin) { t = t1; break; } else { t1 = tmid; f1 = fmid; } } }
		else break; };
	f = eval(pos, t, dir);
	if (f <= best) { pos += t*dir; best = f; };
}

template <typename T>
T powell<T>::optimize(vecD<T>& pos, const vecD<T>& pmin, const vecD<T>& pmax,
	const unsigned int p, const bool verbose)
//! starting at pos, find local optimum within bounds (pmin,pmax), on intervals p, returns optimum in pos.
{	std::vector<vecD<T>> dir(dims);
	for (unsigned int i = 0; i < dims; i++) { dir[i] = pos; dir[i] = 0; dir[i](i) = T(1.0); };
	vmin = pmin; vmax = pmax; start = pos; parts = p;
	vecD<T> savePos = pos; best = goal(pos);
	for (unsigned int it = 0; it < maxit; it++) {
		for (unsigned int i = 0; i < dims; i++)
			minimizeOnInterval(pos, dir[i]);				// find minimum in specified directions
		vecD<T> conj = pos-savePos; T l = norm(conj);				// estimate a conjugate direction
		if (l < eps) break;
		conj /= l; minimizeOnInterval(pos, conj);				// minimize in conjugate direction
		for (unsigned int i = 0; i < dims-1; i++) dir[i] = dir[i+1];		// cycle the directions
		dir[dims-1] = conj;							// add conjugate direction to list
		savePos = pos;								// set parameters for next pass
		if (verbose) { printf("[%2d] %f\r", it, best); fflush(stdout); } };
	if (verbose) printf("\n");
	return best;
}

//! Helper class for genetic optimizer: operations on bit strings

struct bitstring {
	unsigned int N;					//!< length of bit string
	bool	*x;					//!< pointer to data

	void	alloc(const unsigned int _N)						//! allocate space for N bits.
		{ N = _N; x = new bool [N]; }
	void	clear()									//! deallocate space.
		{ delete [] x; x = nullptr; N = 0; }
	void	assign(const bitstring &b)						//! assign from b.
		{ if (N != b.N) { clear(); alloc(b.N); };
		  for (unsigned int i = 0; i < N; i++) x[i] = b.x[i]; }
public:
	bitstring(const unsigned int _N = 0)						//! allocates an empty bitstring.
		 : N(0), x(nullptr) { alloc(_N); }
	bitstring(const bitstring& b)							//! copies from b.
		 : N(0), x(nullptr) { assign(b); }
	bitstring(bitstring&& b)							//! moves from b.
		 : N(b.N), x(b.x) { b.x = nullptr; }
	~bitstring() { clear(); }
	bitstring& operator=(const bitstring& b)					//! assigns from b.
		{ if (this != &b) assign(b); return *this; }
	bitstring& operator=(bitstring&& b)						//! move assigns from b.
		{ assert(this != &b); delete [] x; x = b.x; b.x = nullptr;
		  N = b.N; return *this; }
	bool	get(const unsigned int s, const unsigned int e, unsigned int& b) const	//! returns bits [s,e] in b.
		{ if (s > e || s >= N || e > N || e-s > 32) return false;
	 	  b = 0; for (unsigned int i = s; i < e; i++) { 
			if (i != s) { b <<= 1; }; b |= x[i]; };
		  return true; }
	bool	set(const unsigned int s, const unsigned int e, unsigned int b)		//! sets bits [s,e] from b.
		{ if (s > e || s >= N || e > N || e-s > 32) return false;
	  	  for (unsigned int i = e-1; i >= s; i--) { 
			x[i] = b & 0x1; b >>= 1; if (i == 0) break; };
		  return true; }
	bool	swap(const unsigned int s, const unsigned int e, bitstring& b)		//! swaps bits [s,e] with b.
		{ if (s > e || s >= N || e > N) return false;
	  	  for (unsigned int i = s; i < e; i++) std::swap(x[i],b.x[i]);
		  return true; }
};

#define RANDINT(l,h)	FTOU((l)+(h-l+1)*drand48())
#define RANDMUT(p)	FTOU(std::ceil(std::log(drand48())/std::log(1.0-p)))

//! Helper class for genetic optimizer: implements a member of a population

template<class T> struct member {
	unsigned int 	ng;				//!< number of genes
	unsigned int 	len;				//!< length of gene in bits
	bitstring	genome;				//!< bitstring of genes
	T 		fitness;			//!< fitness (see evalFitness)
	bool 		needsEval;			//!< to avoid unnecessary recalculations

	member<T>(const unsigned int n, const unsigned int l)				//! allocates a member with n genes of length l.
		 : ng(n), len(l), genome(n*l), fitness(std::numeric_limits<T>::max()), needsEval(true)
		{ const unsigned int m = mask();
	  	  for (unsigned int i = 0; i < ng; i++) setGene(i,RANDINT(0,m)); }
	member(const member& m)								//! copies from member m.
		 : ng(m.ng), len(m.len), genome(m.genome), fitness(m.fitness), needsEval(m.needsEval)
		{ }
	~member() { }
	member&	operator=(const member& b)						//! assigns from member m.
		{ if (this != &b) { ng = b.ng; len = b.len; genome = b.genome;
			fitness = b.fitness; needsEval = b.needsEval; };
		  return *this; }
	member&	operator=(member&& b)							//! move assigns from member m.
		{ assert(this != &b); ng = b.ng; len = b.len; fitness = b.fitness; 
		  genome = std::move(b.genome); needsEval = b.needsEval;
		  return *this; }
	unsigned int mask() const							//! constructs a bit mask of length len.
		{ unsigned int m = 0;
		  for (unsigned int i = 0; i < len; i++) { 
			if (i != 0) { m <<= 1u; }; m |= 1u; };
		  return m; }
	void	mutate()								//! makes a random mutation. 
		{ const unsigned int pos = RANDINT(0,ng*len-1), b = RANDINT(0,1);
	 	  genome.set(pos, pos+1, b); needsEval = true; }
	void	crossover(member& m)							//! cross over genes with member m. 
		{ unsigned int s = RANDINT(0,ng*len-2), e = RANDINT(0,ng*len-1);
		  if (e < s) std::swap(s,e); else if (e == s) e++;
		  if (genome.swap(s,e,m.genome)) needsEval = true; }
	void	exchange(member& m)							//! exchange genes with member m. 
		{ unsigned int l = RANDINT(0,ng-1), s = l*len, e = s+l;
		  if (genome.swap(s,e,m.genome)) needsEval = true; }
	unsigned int getGene(const unsigned int i) const				//! returns bitstring section for gene i.
		{ unsigned int b = 0; genome.get(i*len,(i+1)*len,b); return b; }
	void	setGene(const unsigned int i, const unsigned int b)			//! sets bitstring section for gene i from b.
		{ genome.set(i*len,(i+1)*len,b); }
	void	mapValues(vecD<T>& pos) const						//! map genes to parameters, returns pos. 
		{ const T m = T(mask());
		  for (unsigned int i = 0; i < ng; i++) pos[i] = T(getGene(i))/m; }
	void	setFitness(const T f)							//! set fitness to f. 
		{ fitness = f; needsEval = false; }
	T	getFitness() const { return fitness; }					//! returns fitness.
	bool	checkFitness() const { return needsEval; }				//! check if this member needs evaluation of the cost function.
};

//! Implements non-linear minimization in N dimensions using a genetic algorithm

template<class T> class population {
	const unsigned int size;			//!< population size
	const unsigned int maxgen;			//!< maximum generation #
	const T	p_cross;				//!< crossover probability
	const T	p_mut;					//!< mutation probability
	const T	p_exch;					//!< exchange probability
	const T	rankMin;				//!< ranking parameter
	std::vector<member<T>> mem;			//!< actual generation
	std::vector<member<T>> old;			//!< last generation

	void	select();
	virtual T computeFitness(const member<T>& m) = 0;				// interface to cost function.
	virtual void mapValues(const member<T>& m, vecD<T>& pos) = 0;			// maps values from member to position.
public:
	population<T>(const unsigned int _size, const unsigned int ngen, const unsigned int n, 
		const unsigned int l, const T pc, const T pm, const T pe)		//! allocates a genetic optimizer with size members over ngen generations.
		 : size(_size), maxgen(ngen), p_cross(pc), p_mut(pm), p_exch(pe),
		rankMin(T(0.6)), mem(), old()
		{ for (unsigned int i = 0; i < size; i++) mem.push_back(member<T>(n,l)); 
		  for (unsigned int i = 0; i < size; i++) old.push_back(member<T>(n,l)); }
	virtual ~population() { }
	T	generate(vecD<T>& pos, const bool verbose = false);
};

template <typename T>
void population<T>::select()
//! select members for next generation.
{	std::swap(old,mem); uvecD sample(size), rank(size);				// rank = size-1 for best, rank = 0 for worst
	for (unsigned int i = 0; i < size; i++) { unsigned int r = 0;			// find the ith best structure
		for (unsigned int j = 0; j < size; j++)
			if (old[j].fitness < old[i].fitness) r++;
		rank[i] = size-r-1; };							// mark best structure with its rank
	const unsigned int rlim = FTOU(rankMin*size);					// sample best rlim members
	for (unsigned int i = 0; i < size; i++) {
		if (rank[i] > rlim) sample[i] = i;
		else {	unsigned int j = RANDINT(0, size-1);
			for (unsigned int k = 0; k < size; k++) {
				unsigned int l = (k+j)%size; if (rank[l] < rlim) continue;
				j = l; break; };
			sample[i] = j; } };
	for (unsigned int i = 0; i < size-1; i++)	
		std::swap(sample[i],sample[RANDINT(i,size-1)]);				// shuffle
	for (unsigned int i = 0; i < size; i++) mem[i] = old[sample[i]];		// build next generation
}

template <typename T>
T population<T>::generate(vecD<T>& pos, const bool verbose)
//! form next population from the old one via genetic operators.
{	T fmin = std::numeric_limits<T>::max(), fmax = -fmin; member<T> b = mem[0];
	for (unsigned int t = 0; t < maxgen; t++) {					// for maxgen generations...
		select();
		for (unsigned int i = 0; i < FTOU(p_mut*size); i++)			// mutate members...
			mem[RANDINT(1,size-1)].mutate();
		for (unsigned int i = 0; i < FTOU(p_cross*size); i++) {			// cross-over between neighbors...
			const unsigned int j = RANDINT(1,size-1), k = RANDINT(1,size-1);
			if (j != k) mem[j].crossover(mem[k]); };
		for (unsigned int i = 0; i < FTOU(p_exch*size); i++) {			// exchange genes...
			const unsigned int j = RANDINT(1,size-1), k = RANDINT(1,size-1);
			if (j != k) mem[j].exchange(mem[k]); };
		for (unsigned int i = 0; i < size; i++) { T f = mem[i].getFitness();	// evaluate results...
			if (mem[i].checkFitness()) {
				f = computeFitness(mem[i]); mem[i].setFitness(f); };
			if (f < fmin) { fmin = f; b = mem[i]; }; fmax = std::max(f,fmax); };
		if (verbose) { printf("[%3d] %-8.4e %-8.4e\r", t, fmin, fmax); fflush(stdout); } };
       if (verbose) printf("\n");
	mapValues(b,pos); return fmin;							// return result of best member
}
#endif
