#ifndef pFeProblem_H
#define pFeProblem_H

/*
 *
 * pFeProblem.h: basic definitions & classes for finite element models, partitioned version
 * BRIAN Software Package Version 3.0
 *
 * $Id: pFeProblem.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (18/02/11): first version UH&  FK
 * 0.20 (14/05/11): first working FK
 * 0.30 (23/10/11): first release for BRIAN 2.2
 * 0.40 (31/12/12): released version 2.4
 * 0.50 (22/11/14): partitioned version 2.7
 * v406 (28/09/16): bumped to version 3.0
 *
 */

//! Defines a boundary condition for FE problems

template <typename T> struct condition {
	unsigned int id;			//!< condition type
	vecD<T>	data;				//!< vector of condition data

	condition(unsigned int _id, vecD<T>& _data)					//! stores a boundary condition.
		: id(_id), data(_data) { }
};

/*! \file
    \brief Basic definitions & classes for finite element models, partitioned version.
*/

//! Base class for solving FE problems

template <typename T> class pFeProblem {
protected:
	const unsigned int length = 1024;	//!< maximum line length
	volumeMesh& m;				//!< volumetric mesh
	const ptmap pm;				//!< partition map
 	const unsigned int nv;			//!< number of vertices
 	const unsigned int np;			//!< number of primitives
	const unsigned int dof;			//!< degrees of freedom per node
	const unsigned int dim;			//!< dimension of equation system
	const T eps;				//!< solution tolerance
	const partition pt;			//!< partition
	pmatS<T> K;				//!< stiffness matrix
	pvecD<T> u;				//!< solution
	pvecD<T> f;				//!< right-hand side
	std::vector<condition<T>> fixedNodes;	//!< vector of fixed nodes
	std::vector<condition<T>> loadedNodes;	//!< vector of loaded nodes
	std::vector<condition<T>> timepoints;	//!< vector of timepoints

	bool	readNodeId(unsigned int& id) const					//! read node id.
		{ char* s = strtok(nullptr,"\t\n "); id = s? ATOU(s): 0;
		  return s != nullptr; }
	bool	readRHS(vecD<T>& v) const						//! read RHS.
		{ char* s = nullptr; unsigned int i = 0; v = 0;
		  for allDOF(d) { s = strtok(nullptr, "\t\n ");
			if (s) { v[i++] = T(atof(s)); } };
		  return i == dof; }
	bool	readTimepoint(T& t, T& v) const						//! force factor at time point.
		{ char* s = strtok(nullptr,"\t\n "); t = s? T(atof(s)): T(0);
		  s = strtok(nullptr,"\t\n "); v = s? T(atof(s)): T(0);
		  return s != nullptr; }
	void	readConditions(const char* bc)						//! read boundary conditions from file* bc.
		{ getFixedNodes(bc); getLoadedNodes(bc); getTimepoints(bc); }
	unsigned int eqno(const unsigned int id)					//! returns equation number for node id.
		{ assert(id < nv); return pm.toEq(id)*dof; }
	unsigned int idno(const unsigned int eq)					//! returns node id for equation number eq.
		{ assert(eq < dim); return pm.toId(eq/dof); }
	void	getFixedNodes(const char* bc);
	void	getLoadedNodes(const char* bc);
	void	getTimepoints(const char* bc);
	element<T>* createElement(const cell& c) const;
	unsigned int initRow(vecS<T>& row, const unsigned int id);
	void	allocateLHS(pmatS<T>& A);
	void	updateMatrix(pmatS<T>& A, const cell *c, const dmatD& k);
	void	setFixedNodes(pmatS<T>& A);
	void	setRHS(const T fac);
	unsigned int inPartition(const cell* c);
public:
	pFeProblem(volumeMesh& _m, const unsigned int _dof, const T _eps, const char *bc, const limage& pi)
	//! allocates a context for solving FE problems.
	/*! Parameters are: m refers to the volumetric mesh, dof to the degrees of
	    freedom per node, eps to the solution tolerance, and bc refers to a file
	    with boundary conditions.
	 */
		 : m(_m), pm(pi), nv(m.nVertices()), np(m.nPrimitives()), dof(_dof), dim(dof*nv), eps(_eps),
		pt(pm.offsets()*int(dof)), K(), u(), f(), fixedNodes(), loadedNodes(), timepoints()
		{ m.vertexToPrimitiveMap(); allocateLHS(K); readConditions(bc); }
	~pFeProblem() { }
	void	solve(const bool verbose = false)					//! solve a FE system.
		{ pprecBlockJacobi<T> pre(K,dof); pgmres<T> sv(K,pre,20,eps); sv.solve(u,f,verbose); }
};

template <typename T>
element<T>* pFeProblem<T>::createElement(const cell& c) const
//! create an element for cell c.
{	switch (c.vtx.N) {
	case  4: return new tetElement4<T>(c, m, dof);
	case 10: return new tetElement10<T>(c, m, dof);
	case  8: return new hexElement8<T>(c, m, dof);
	case 20: return new hexElement20<T>(c, m, dof);
	default: throw optException("Invalid element type %u", int(c.vtx.N)); };
}

template <typename T>
void pFeProblem<T>::getFixedNodes(const char* bc)
//! read set of fixed nodes from file *bc.
{	FILE *fp = openFile(bc, "rb"); char line[length];
	unsigned int id = 0; vecD<T> fx(dof); fx = 0;
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {							// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "fix", 3) == 0) { 
			if (readNodeId(id)) fixedNodes.push_back({id,fx});
			else throw optException("Illegal information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void pFeProblem<T>::getLoadedNodes(const char* bc)
//! read set of nodes with loads from file *bc.
{	FILE *fp = openFile(bc, "rb"); char line[length];
	unsigned int id = 0; vecD<T> fx(dof); fx = 0;
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {							// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "mov", 3) == 0) {
			if (readNodeId(id) && readRHS(fx)) loadedNodes.push_back({id,fx});
			else throw optException("Illegal information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void pFeProblem<T>::getTimepoints(const char* bc)
//! read set of time points from file *bc.
{	FILE *fp = openFile(bc, "rb"); char line[length];
	unsigned int id = 0; vecD<T> tp(2); tp = 0;
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {							// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "tp", 2) == 0) {
			if (readTimepoint(tp[0], tp[1])) timepoints.push_back({id++,tp});
			else throw optException("Illegal information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void pFeProblem<T>::setRHS(const T fac)
//! set boundary conditions in RHS.
{	const unsigned int s = K.p.st(), e = K.p.end();
	u.repartition(K.p); u = 0; f.repartition(K.p); f = 0;				// init u & f using K's partitioning scheme
	for (unsigned int i = 0; i < fixedNodes.size(); i++) {				// for all fixed BCs
		const unsigned int eq = eqno(fixedNodes[i].id);
		for allDOF(d) {	if (eq+d < s || eq+d >= e) continue;			// skip if outside of this partition
			f[eq-s+d] = 0.0; } };						// set force to 0
	for (unsigned int i = 0; i < loadedNodes.size(); i++) {				// apply loads to RHS
		condition<T>& c = loadedNodes[i]; const unsigned int eq = eqno(c.id);
		for allDOF(d) {	if (eq+d < s || eq+d >= e) continue;			// skip if outside of this partition
			f[eq-s+d] += c.data[d]*fac; } };				// update force vector
}

template <typename T>
unsigned int pFeProblem<T>::initRow(vecS<T>& row, const unsigned int id)
//! initialize a row of the stiffness matrix.
{	const auto vt = m.verticesAt(id);  						// collect a list of neighbor vertices including this one
	row.resize(vt.size()*dof); unsigned int j = 0;
	for (const auto i : vt)	{ for allDOF(d) row.set(j++, eqno(i)+d, T(0)); };	// set coefficients for rows 
	return row.N;
}

template <typename T>
void pFeProblem<T>::allocateLHS(pmatS<T>& A)
//! allocate entries in partitioned sparse matrix A for mesh m.
{	const unsigned int s = pt.st(), e = pt.end();
	vecS<T> row; unsigned int nz = 0, pid = 0, pel = 0;
	for (unsigned int eq = s; eq < e; eq += dof) {					// for all equations in this partition
		for allDOF(d) {	if (eq+d < s || eq+d >= e) continue;			// skip if outside of this partition
			const unsigned int id = idno(eq+d);
			if (id != pid) { pel = initRow(row,id); pid = id; };		// allocate dof rows per node 
			nz += pel; } };
	A.resize(dim, dim, nz, pt);							// allocate stiffness matrix
	for (unsigned int eq = s; eq < e; eq += dof) {
		for allDOF(d) {	if (eq+d < s || eq+d >= e) continue;			// skip if outside of this partition
			const unsigned int id = idno(eq+d);
			if (id != pid) { pel = initRow(row,id); pid = id; };		// allocate dof rows per node 
			A.initRow(eq+d,row); } };					// allocate sparse row times dof
}

template <typename T>
void pFeProblem<T>::updateMatrix(pmatS<T>& A, const cell* c, const dmatD& k)
//! update global stiffness matrix A using local matrix k and cell c.
{	const unsigned int s = pt.st(), e = pt.end();
	uvecD eq(c->vtx.N*dof); unsigned int n = 0;					// allocate vector of eq #
	for (unsigned int i = 0; i < c->vtx.N; i++)					// for all vertices of this cell
		for allDOF(d) { eq[n++] = eqno(c->vtx(i))+d; };				// compute equation # of each dof
	for (unsigned int i = 0; i < eq.N; i++) {					// assemble k into global matrix A
		if (eq[i] < s || eq[i] >= e) continue;					// row is outside of this partition
		for (unsigned int j = 0; j < eq.N; j++) A(eq[i]-s,eq[j]) += k(i,j); }
}

template <typename T>
void pFeProblem<T>::setFixedNodes(pmatS<T>& A)
//! set rows & columns for fixed nodes in stiffness matrix.
{	const unsigned int s = pt.st(), e = pt.end(); bvecD fix(A.N); fix = false;
	for (unsigned int i = 0; i < fixedNodes.size(); i++) {				// for all fixed BCs
		const unsigned int eq = eqno(fixedNodes[i].id);
		for allDOF(d) {	fix[eq+d] = true;
			if (eq+d < s || eq+d >= e) continue;				// skip if outside of this partition
			for (unsigned int t = A.ri[eq-s+d]; t < A.ri[eq-s+d+1]; t++) {	// set row to identity
				unsigned int c = A.ci[t]; A.x[t] = (c == eq+d)? 1.0: 0.0; } } };
	for (unsigned int r = s; r < e; r++) {						// and zero out columns
		for (unsigned int t = A.ri[r-s]; t < A.ri[r-s+1]; t++) {
			const unsigned int c = A.ci[t];
			if (fix(c) && r != c) A.x[t] = 0; } };
}

template <typename T>
unsigned int pFeProblem<T>::inPartition(const cell* c)
//! returns # of vertices inside this partition.
{	const unsigned int s = pt.st(), e = pt.end(); unsigned int n = 0;		// init row for node id
	for (unsigned int i = 0; i < c->vtx.N; i++) {					// for all vertices of this cell
		const unsigned int eq = eqno(c->vtx(i));				// compute equation #
		for allDOF(d) {	if (eq+d >= s && eq+d < e) n++;	} };			// count eqs inside this partition
	return n;
}

#endif
