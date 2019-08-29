#ifndef mincutOPTIMIZER_H
#define mincutOPTIMIZER_H

/*
 * mincutOptimizer.h: optimization using the mincut algorithm with multiple labels
 * BRIAN Software Package Version 3.0
 *
 * $Id: mincutOptimizer.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (21/06/15): initial version
 * 0.11 (12/07/15): debugged & templated
 * v406 (28/09/16): bumped to version 3.0
 *
 * This code is based on an implementation by Olga Veksler and Andrew Delong
 * but was completely re-written. Find the original copyright below.
 *
 */

/* ##################################################################

2. License & disclaimer.

    Copyright 2007-2010 Olga Veksler  <olga@csd.uwo.ca>
                        Andrew Delong <andrew.delong@gmail.com>

    This software and its modifications can be used and distributed for 
    research purposes only. Publications resulting from use of this code
    must cite publications according to the rules given above. Only
    Olga Veksler has the right to redistribute this code, unless std::expressed
    permission is given otherwise. Commercial use of this code, any of 
    its parts, or its modifications is not permited. The copyright notices 
    must not be removed in case of any modifications. This Licence 
    commences on the date it is electronically or physically delivered 
    to you and continues in effect unless you fail to comply with any of 
    the terms of the License and fail to cure such breach within 30 days 
    of becoming aware of the breach, in which case the Licence automatically 
    terminates. This Licence is governed by the laws of Canada and all 
    disputes arising from or relating to this Licence must be brought 
    in Toronto, Ontario.


    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY std::expRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

##################################################################
*/

/*! \file
    \brief Templated max flow graph and energy optimization classes.

    \details For detailed information, refer to
    (\cite Boykov04, <a href="ref58.pdf" target="_blank">PDF</a>).
*/

using uvec = std::vector<unsigned int>;

//! Provides an optimizer using min cut / max flow problems on general graphs.

template<typename T>
class mincutOptimizer {
	const unsigned int nl;			//! number of labels (classes)
	const unsigned int ns;			//! number of sites (nodes)
	uvecD	lbl;				//! vector of labels
	uvecD	siteVar;			//! holds variable indices during a move
	uvecD	table;    			//! label order for std::expansion and swaps
	T	energy;				//! total energy
	bool	random;				//! true if labels are randomized

	void	addTerm1(mincutEnergy<T>& e, const unsigned int i, const T c0, const T c1) //! adds an unary term to the energy function.
		{ energy += c1; e.addTerm1(i,c0,c1); }
	void	addTerm2(mincutEnergy<T>& e, const unsigned int i, const unsigned int j,
			const T c0, const T c1, const T c2, const T c3)			//! adds a binary term to the energy function.
		{ energy += c3; e.addTerm2(i,j,c0,c1,c2,c3); }
	uvec	getActiveSites(const unsigned int a) const				//! returns a vector of active sites.
		{ uvec act; for (unsigned int s = 0; s < ns; s++)
			if (lbl(s) != a) act.push_back(s);
		  return act; }
	void	permuteTable()								//! performs a random mutation of current labeling.
		{ if (random == false) return; uniformDist<T> ud;
		  for (unsigned int l = 0; l < nl; l++) {
			const unsigned int n = std::floor(ud()*(nl-l));
			std::swap(table[l],table[n+l]); } }
	T	oneExpansion()								//! permutes labels, performs a single std::expansion over all labels, and re-computes energy.
		{ permuteTable();
		  for (unsigned int l = 0; l < nl; l++) alphaExpansion(table[l]);	// each cycle is exactly one pass over the labels
		  return computeEnergy(); }
	T	oneSwap()								//! permutes labels, performs a single swap over all label pairs, and re-computes energy.
		{ permuteTable();
		  for (unsigned int l = 0; l < nl; l++)
			for (unsigned int l1 = nl-1;  l1 != std::numeric_limits<unsigned int>::max(); l1--)
				if (table[l] < table[l1]) alphaBetaSwap(table[l], table[l1]);
		  return computeEnergy(); }
	bool	optimizeExpansion(const unsigned int a, const uvec& act);
	bool	alphaExpansion(unsigned int a);
	void	optimizeSwap(const unsigned int a, const unsigned int b, const uvec& act);
	void	alphaBetaSwap(unsigned int a, unsigned int b);
	virtual uvec getNeighbors(const unsigned int i) = 0;
protected:
	void	setLabelOrder(const bool r)						//! sets a random labeling.
		{ random = r; for (unsigned int l = 0; l < nl; l++) table[l] = l; }
	const uvecD& getLabels() const							//! accesses current labeling.
		{ return lbl; }
	void	setLabelOrder(const uvecD& ord)						//! sets labels from vector ord.
		{ if (ord.N > nl) printf("setLabelOrder has too many labels.\n");
		  random = false; table = std::numeric_limits<unsigned int>::max();
		  for (unsigned int i = 0; i < ord.N; i++) {
			if (ord(i) >= nl) printf("Invalid label %u in setLabelOrder.\n", ord(i));
			table[i] = ord(i); } }
	T	computeEnergy()								//! computes total energy of this labeling.
		{ T e = 0; for (unsigned int i = 0; i < ns; i++) {			// for all sites i
			e += unaryCosts(i,lbl[i]); uvec nb = getNeighbors(i); 		// compute unary cost term; get neighbors
			for (const auto j : nb) {					// for all neighbors j
				if (j < i) e += binaryCosts(i,j,lbl[i],lbl[j]); } };	// and sum up binary cost terms
		  return e; }
	T	expansion(const unsigned int maxit = std::numeric_limits<unsigned int>::max());
	T	swap(const unsigned int maxit = std::numeric_limits<unsigned int>::max());
public:
	mincutOptimizer(const unsigned int _ns, const unsigned int _nl) 
		: nl(_nl), ns(_ns), lbl(ns), siteVar(ns), table(nl),
		energy(0), random(false)						//! provides a context for optimization using the minimum cut principle.
		{ lbl = 0; siteVar = std::numeric_limits<unsigned int>::max(); setLabelOrder(false); }
	virtual ~mincutOptimizer() { }
	virtual T unaryCosts(const unsigned int i, const unsigned int a) const = 0;
	virtual T binaryCosts(const unsigned int i, const unsigned int j,
			const unsigned int a, const unsigned int b) const = 0;
};

template<typename T>
bool mincutOptimizer<T>::optimizeExpansion(const unsigned int a, const uvec& act)
//! optimizes expansion for label a on active sites in act.
{	energy = 0; const unsigned int n = act.size();
	mincutEnergy<T> e; e.addVariable(n);						// create binary variables for each remaining site
	for (unsigned int s = 0; s < n; s++) { const unsigned int i = act[s];		// add the data costs
		addTerm1(e,s,unaryCosts(i,a),unaryCosts(i,lbl[i])); };
	for (unsigned int s = n-1; s < n; s--) {					// and compute the smooth costs between variables.
		const unsigned int i = act[s]; const auto nb = getNeighbors(i);
		for (const auto j : nb) { const unsigned int b = lbl[j];
			if (siteVar[j] == std::numeric_limits<unsigned int>::max()) {
				addTerm1(e,s,binaryCosts(i,j,a,b),binaryCosts(i,j,b,b)); }
			else if (j < i) {
				const T c0 = binaryCosts(i,j,a,a), c1 = binaryCosts(i,j,a,b);
				const T c2 = binaryCosts(i,j,b,a), c3 = binaryCosts(i,j,b,b);
				addTerm2(e,s,siteVar[j],c0,c1,c2,c3); } } };
	const T en = e.minimize(); const bool succ = (en/energy)+1e-6 < 1.0;
	if (succ) { for (unsigned int s = 0; s < act.size(); s++) {			// apply new labeling
		if (e.getVar(s) == 0) lbl[act[s]] = a; } };
	return succ;
}

template<typename T>
bool mincutOptimizer<T>::alphaExpansion(const unsigned int l)
//! performs a single alpha expansion for label l.
{	if (l == std::numeric_limits<unsigned int>::max()) return false;
	const auto act = getActiveSites(l); if (act.size() == 0) return false;		// get list of active sites
	unsigned int i = 0; for (const auto s : act) siteVar[s] = i++;			// initialize reverse-lookup
	const bool succ = optimizeExpansion(l,act);
	for (const auto s : act) siteVar[s] = std::numeric_limits<unsigned int>::max();					// restore siteVar
	return succ;
}

template<typename T>
T mincutOptimizer<T>::expansion(const unsigned int maxit)
//! performs an iterative expansion over all labels.
{	permuteTable();
	if (maxit == std::numeric_limits<unsigned int>::max()) { uvec qs; qs.push_back(nl);				// loop over labels that successfully reduce energy
		for (unsigned int n = 0; qs.empty() == false; ) {
			unsigned int nq = qs.back(), s = n;				// pass over the unchecked labels in the current queue
			while (n < nq) {
				if (alphaExpansion(table[n])) n++;
				else std::swap(table[n],table[--nq]); }			// don't put this label in a new queue
			if (n == s) { n = qs.back(); qs.pop_back(); }			// no success: try more labels from the previous queue
			else if (nq < qs.back()/2) { n = 0; qs.push_back(nq); }		// some success: focus on them in a new queue
			else n = 0; };							// all completed: make another sweep
		return computeEnergy(); }
	else {	T ne = computeEnergy(), oe = ne+1;
		for (unsigned int it = 0; it < maxit; it++) {
			oe = ne; ne = oneExpansion(); 
			if (ne == oe) return ne; } };
	return 0;
}

template<typename T>
void mincutOptimizer<T>::optimizeSwap(const unsigned int a, const unsigned int b, const uvec& act)
//! optimizes swap for labels (a,b) on active sites in act.
{	energy = 0; const unsigned int n = act.size();
	mincutEnergy<T> e; e.addVariable(n);						// create binary variables for each remaining site
	for (unsigned int s = 0; s < act.size(); s++) { const unsigned int i = act[s];	// add the data costs,
		const T c0 = unaryCosts(i,a), c1 = unaryCosts(i,b);
		addTerm1(e,s,c0,c1); }
	for (unsigned int s = n-1; s != std::numeric_limits<unsigned int>::max(); s--) {				// and compute the smooth costs between variables.
		const unsigned int i = act[s]; const auto nb = getNeighbors(i);
		for (const auto j : nb) { const unsigned int l = lbl[j];
			if (siteVar[ns] == std::numeric_limits<unsigned int>::max()) {
				const T c0 = binaryCosts(i,j,a,l);
				const T c1 = binaryCosts(i,j,b,l);
				addTerm1(e,s,c0,c1); }
			else if (ns < s) {
				const T c0 = binaryCosts(i,j,a,a);
				const T c1 = binaryCosts(i,j,a,b);
				const T c2 = binaryCosts(i,j,b,a);
				const T c3 = binaryCosts(i,j,b,b);
				addTerm2(e,s,siteVar[j],c0,c1,c2,c3); } } };
	e.minimize();
	unsigned int i = 0; for (const auto s : act) {					// apply new labeling
		lbl[s] = e.getVar(i) == 0? a: b; siteVar[s] = std::numeric_limits<unsigned int>::max(); i++; }
}

template<typename T>
void mincutOptimizer<T>::alphaBetaSwap(const unsigned int a, const unsigned int b)
//! performs a single alpha-beta swap for labels (a,b).
{	assert(a < nl && b < nl); unsigned int n = 0; uvec act;				// determine the list of active sites for this swap
	for (unsigned int s = 0; s < ns; s++) {
		if (lbl[s] == a || lbl[s] == b) { act.push_back(s); siteVar[s] = n++; } }
	if (n) optimizeSwap(a,b,act);							// if any found, perform swap
}

//! Provides an neighborhood iterator that returns a list of relative indices to valid neighbors.

class indexIterator : public nbIterator {
public: 
	indexIterator(const uvec3& ex, const connectivity cn)				//! allocates a index iterator for extent ex and connectivity cn.
		 : nbIterator(ex, cn) { }
	uvec	neighbors(const unsigned int t)						// returns a vector of indices relative to t
		{ setloc(t); uvec nb; for (i = 0; i < n; i++) { 			// use i from nbIterator
			const unsigned int sx = ITOU(int(loc.x)+dx[i]); if (sx >= ext.x) continue;
			const unsigned int sy = ITOU(int(loc.y)+dy[i]); if (sy >= ext.y) continue;
			const unsigned int sz = ITOU(int(loc.z)+dz[i]); if (sz >= ext.z) continue;
			nb.push_back((sz*ext.y+sy)*ext.x+sx); }; return nb; }
};

//! Provides an optimizer using min cut / max flow problems on images.

template<typename T>
class mincutImageOptimizer : public mincutOptimizer<T> {
	indexIterator ix;
	uvec	getNeighbors(const unsigned int i)					//! provides access to neighbors for the mincutOptimizer.
		{ return ix.neighbors(i); }
protected:
	const uvec3 ex;
public:
	mincutImageOptimizer(const uvec3 _ex, const connectivity cn, const unsigned int nl)
		: mincutOptimizer<T>(_ex.x*_ex.y*_ex.z,nl), ix(_ex,cn), ex(_ex)		//! provides a context for optimizing image energies.
		{ }
	limage	work(const unsigned int nit = 4, const bool verbose = false)		//! optimizes an image labeling using alpha expansion.
		{ const T e0 = this->computeEnergy(); this->expansion(nit);
		  if (verbose) printf("Energy %.3e -> %.3e\n", e0, this->computeEnergy());
		  limage dst(ex); dst = 0; const uvecD& lbl = this->getLabels();
		  for (unsigned int i = 0; i < lbl.N; i++) dst(i) = lbl(i);
		  return dst; }
};

//! Provides an optimizer using min cut / max flow problems on graphs.

template<typename T, typename N>
class mincutGraphOptimizer : public mincutOptimizer<T> {
	const graph<N*>& g;
	uvec	getNeighbors(const unsigned int i)					//! provides access to neighbors for the mincutOptimizer.
		{ uvec nb; const N* n = g(i); if (n) {
			for (const auto& e : n->edges) nb.push_back(e.id); };
		  return nb; }
public:
	mincutGraphOptimizer(const graph<N*>& _g, const unsigned int nl)
		: mincutOptimizer<T>(_g.size(),nl), g(_g)				//! provides a context for optimizing graph energies.
		{ }
	uvecD	work(const unsigned int nit = 4, const bool verbose = false)		//! optimizes an image labeling using alpha expansion.
		{ const T e0 = this->computeEnergy(); this->expansion(nit);
		  if (verbose) printf("Energy %.3e -> %.3e\n", e0, this->computeEnergy());
		  return this->getLabels(); }
};
#endif
