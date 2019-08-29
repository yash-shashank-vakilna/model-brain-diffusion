#ifndef FEPROBLEM_H
#define FEPROBLEM_H

/*
 *
 * feProblem.h: basic definitions & classes for finite element models
 * BRIAN Software Package Version 3.0
 *
 * $Id: feProblem.h 506 2017-03-19 01:39:21Z frithjof $
 *
 * 0.10 (18/02/11): first version UH&  FK
 * 0.20 (14/05/11): first working FK
 * 0.30 (23/10/11): first release for BRIAN 2.2
 * 0.40 (31/12/12): released version 2.4
 * 0.50 (21/11/14): templated & joined with feProblem.C
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Basic definitions & classes for finite element models.
*/

//! Defines a boundary condition for FE problems

template <typename T> struct condition {
	unsigned int id;			//!< condition type
	vecD<T>	data;				//!< vector of condition data

	condition(unsigned int _id, vecD<T>& _data)					//! stores a boundary condition.
		: id(_id), data(_data) { }
};

//! Base class for solving FE problems

template <typename T> class feProblem {
public:
	const unsigned int length = 1024;	//!< maximum line length
	volumeMesh& m;				//!< volumetric mesh
 	const unsigned int nv;			//!< number of vertices
 	const unsigned int np;			//!< number of primitives
	const uvecD map;			//!< vertex mapping
	unsigned int dof;			//!< degrees of freedom per node
	unsigned int dim;			//!< dimension of equation system
	const T eps;				//!< solution tolerance
	matS<T> K;				//!< stiffness matrix
	vecD<T> u;				//!< solution
	vecD<T> f;				//!< right-hand side
	std::vector<condition<T>> fixedNodes;	//!< vector of fixed nodes
	std::vector<condition<T>> loadedNodes;	//!< vector of loaded nodes
	std::vector<condition<T>> timepoints;	//!< vector of timepoints

	bool	readNodeId(unsigned int& id) const					//! read node id.
		{ char* s = strtok(nullptr,"\t\n "); id = s? ATOU(s): 0;
		  if (id >= nv) throw rtException("invalid node id %u in boundary condition", id);
		  return s != nullptr; }
	bool	readRHS(vecD<T>& v) const						//! read RHS.
		{ char* s = nullptr; unsigned int i = 0; v = T(0);
		  for allDOF(d) { s = strtok(nullptr, "\t\n ");
			if (s) { v[i++] = T(atof(s)); } };
		  return i <= dof; }
	bool	readTimepoint(T& t, T& v) const						//! force factor at time point.
		{ char* s = strtok(nullptr,"\t\n "); t = s? T(atof(s)): T(0);
		  s = strtok(nullptr,"\t\n "); v = s? T(atof(s)): T(0);
		  return s != nullptr; }
	void	readConditions(const char* bc)						//! read boundary conditions from file* bc.
		{ getFixedNodes(bc); getLoadedNodes(bc); getTimepoints(bc); }
	unsigned int eqno(unsigned int id) const					//! returns equation number for node id.
		{ assert(id < nv); return id*dof; }
	void	getFixedNodes(const char* bc);
	void	getLoadedNodes(const char* bc);
	void	getTimepoints(const char* bc);
	element<T>* createElement(const cell& c) const;
	unsigned int initRow(vecS<T>& row, const unsigned int id);
	void	allocateLHS(matS<T>& A);
	void	updateMatrix(matS<T>& A, const cell* c, const matD<T>& k);
	void	clearRC(matS<T>& A, const unsigned int id);
public:
	feProblem(volumeMesh& _m, const unsigned int _dof, const T _eps, const char *bc, const uvecD& _map)
		//! allocates a context for solving FE problems, mapped version.
		/*! Parameters are: m refers to the volumetric mesh, dof to the degrees of
		    freedom per node, eps to the solution tolerance, and bc refers to a file
		    with boundary conditions.
		 */
		 : m(_m), nv(m.nVertices()), np(m.nPrimitives()), map(_map), dof(_dof), dim(dof*nv),
		eps(_eps), K(), u(), f(), fixedNodes(), loadedNodes(), timepoints()
		{ m.vertexToPrimitiveMap(); allocateLHS(K);
		  f.resize(dim); f = T(0); readConditions(bc); }
	virtual ~feProblem() { }
	void	applyBoundaryConditions(const T fac);
	void	applyBoundaryConditions(const T fac, const uvecD& map);
	void	solveSystem(const bool verbose = false)					//! solve a FE system.
		{ 
//		  precBlockJacobi<T> pre(K,dof); 
//		  precAMG<T> pre(K,T(0.01)); 
//		  precSAS<T> pre(K,T(0.1));
		  precRAS<T> pre(K,m,40,1);
		  u = cg<T>(K,pre,eps).solve(f,verbose); }
};

template <typename T>
element<T>* feProblem<T>::createElement(const cell& c) const
//! create an element for cell c.
{	switch (c.vtx.N) {
	case  4: return new tetElement4<T>(c,m,dof);
	case 10: return new tetElement10<T>(c,m,dof);
	case  8: return new hexElement8<T>(c,m,dof);
	case 20: return new hexElement20<T>(c,m,dof);
	default: throw optException("Invalid element type %u", int(c.vtx.N)); };
}

template <typename T>
void feProblem<T>::getFixedNodes(const char* bc)
//! read set of fixed nodes from file* bc.
{	FILE *fp = openFile(bc, "rb"); char line[length]; unsigned int id; vecD<T> fx;
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {						// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "fix", 3) == 0) { 
			if (readNodeId(id)) fixedNodes.push_back({id,fx});
			else throw optException("Invalid information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void feProblem<T>::getLoadedNodes(const char* bc)
//! read set of nodes with loads from file* bc.
{	FILE *fp = openFile(bc, "rb"); char line[length];
	unsigned int id = 0; vecD<T> fx(dof); fx = T(0);
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {						// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "mov", 3) == 0) {
			if (readNodeId(id) && readRHS(fx)) loadedNodes.push_back({id,fx});
			else throw optException("Invalid information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void feProblem<T>::getTimepoints(const char* bc)
//! read set of time points from file* bc.
{	FILE *fp = openFile(bc, "rb"); char line[length];
	unsigned int id = 0; vecD<T> tp(2); tp = T(0);
	if (fp == 0) throw optException("Cannot open file %s\n", bc);
	while (fgets(line, length, fp)) {						// first run: determine matrix size
		char *s = strtok(line, "\t\n "); if (s == nullptr || *s == '#') continue;
		if (strncmp(s, "tp", 2) == 0) {
			if (readTimepoint(tp[0], tp[1])) timepoints.push_back({id++, tp});
			else throw optException("Invalid information on line (%s)", line); } };
	closeFile(fp);
}

template <typename T>
void feProblem<T>::applyBoundaryConditions(const T fac)
//! apply fixed boundary conditions to stiffness matrix, mapped version.
{	for (const auto& n : fixedNodes) { 						// for all fixed BCs
		const unsigned int id = map.N? map(n.id): n.id; if (id == ~0u) continue;
		vecS<T> row; initRow(row,n.id); unsigned int eq = eqno(id);		// init row for node id
		for allDOF(d) { row = 0; row(eq+d) = T(1);				// for all dofs
			K.setRow(eq+d, row); K.setCol(eq+d, row); };			// set row&  column to identity
                for allDOF(d) f[eqno(id)+d] = T(0); };					// set RHS to 0
        for (const auto& n : loadedNodes) {						// apply loads to RHS
		const unsigned int id = map.N? map(n.id): n.id; if (id == ~0u) continue;
                for allDOF(d) f[eqno(id)+d] += n.data(d)*fac; };			// add values from force vector
}

template <typename T>
unsigned int feProblem<T>::initRow(vecS<T>& row, const unsigned int id)
//! initialize a row of the stiffness matrix.
{       const auto vt = m.verticesAt(id); unsigned int n = 0, j = 0; 			// collect a list of neighbor vertices including this one
	if (map.N == 0) { row.resize(vt.size()*dof);					// un-mapped version
		for (const auto i : vt)							// for this node and its neighbors
			for allDOF(d) row.set(j++, eqno(i)+d, T(0)); }			// set coefficients for rows 
        else {	for (const auto i : vt) { if (map(i) != ~0u) n++; }			// count mapped vertices
		if (n == 0) return 0;
		row.resize(n*dof);
		for (const auto i : vt) { 						// for the node and its neighbors
			const unsigned int id = map(i); if (id == ~0u) continue;	// skip if not mapped here
                	for allDOF(d) row.set(j++, eqno(id)+d, 0); } };			// init this column for all dofs
        return row.N;
}

template <typename T>
void feProblem<T>::allocateLHS(matS<T>& A)
//! allocate entries in sparse matrix A for mesh m.
{       vecS<T> row; unsigned int nz = 0; dim = 0;
	for (unsigned int i = 0; i < nv; i++) {	if (map.N && map(i) == ~0u) continue;	// for all nodes
		const unsigned int n = initRow(row, i);					// get nz elements (0 means no equation)
		if (n) { dim += dof; nz += n*dof; } };					// allocate dof rows per node
        A.resize(dim, dim, nz);								// allocate matrix
        for (unsigned int i = 0; i < nv; i++) {
		const unsigned int id = map.N? map(i): i; if (id == ~0u) continue;
		initRow(row, i);
                for allDOF(d) A.initRow(eqno(id)+d, row); };				// allocate sparse row times dof
printf("LHS alloc: %u %u\n", dim, nz);
}

template <typename T>
void feProblem<T>::updateMatrix(matS<T>& A, const cell* c, const matD<T>& k)
//! update global stiffness matrix A using local matrix k and cell c.
{	uvecD eq(c->vtx.N*dof); unsigned int n = 0;					// allocate vector of eq #
	for (unsigned int i = 0; i < c->vtx.N; i++) {					// for all vertices of this cell
		const unsigned int j = c->vtx(i);
		const unsigned int id = map.N? map(j): j; if (id == ~0u) continue;
		for allDOF(d) { eq[n++] = eqno(id)+d; } };				// compute equation # of each dof
	assert(k.n() == n);								// just to be sure...
	for (unsigned int i = 0; i < eq.N; i++) 					// assemble k into global matrix A
		for (unsigned int j = 0; j < eq.N; j++) A(eq[i],eq[j]) += k(i,j);
}

template <typename T>
void feProblem<T>::clearRC(matS<T>& A, const unsigned int id)
//! set rows & columns for node id to identity.
{	vecS<T> row; initRow(row,id); unsigned int eq = eqno(id);			// init row for node id
	for allDOF(d) { row = 0; row(eq+d) = T(1);					// for all dofs
		A.setRow(eq+d, row); A.setCol(eq+d, row); };				// set row&  column to identity
}

#endif
