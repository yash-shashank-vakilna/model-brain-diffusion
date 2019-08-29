#ifndef VGRAPH_H
#define VGRAPH_H

/*
 *
 * vectorGraph.h: graph structure representing voxel-wise vector properties
 * BRIAN Software Package Version 2.5
 *
 * $Id: vectorGraph.h 351 2016-05-16 23:16:45Z kruggel $
 *
 * 0.10 (26/06/13): initial version FK
 * 0.20 (23/09/13): renamed & integrated into bv
 * 0.40 (16/12/13): documented
 *
 */

/*! \file
    \brief Implements a graph structure representing voxel-wise vector properties.
*/

inline void vec2mat(fmatD& D, const fvecD& t)						// convert vector to tensor
{ D(0,0) = t(0); D(0,1) = t(1); D(0,2) = t(2);
  D(1,0) = t(1); D(1,1) = t(3); D(1,2) = t(4);
  D(2,0) = t(2); D(2,1) = t(4); D(2,2) = t(5); }

//! Implements a fiber descriptor record.

struct fiberDesc {				// 5 floats
	float	pf;				//!< volume fraction
	float	kappa;				//!< concentration
	fvec3	mu;				//!< direction

	fiberDesc(const float _pf = 0.0f, const float _kappa = 0.0f, const fvec3& _mu = 0.0f)
		//! allocates a fiber descriptor for fraction pf, concentration kappa and direction mu.
		: pf(_pf), kappa(_kappa), mu(_mu) { }
	void	print()									//! prints the record.
		{ printf("pf %e kappa %e mu %e %e %e\n", pf, kappa, mu.x, mu.y, mu.z); }
};

//! Implements a tensor record.

struct tensorDesc : public node {		// 6 floats
	float	v[6];				//!< tensor components xx, xy, xz, yy, yz, zz

	tensorDesc()
		{ for (unsigned int i = 0; i < 6; i++) v[i] = 0.0f; }
	void	print()									//! prints the record.
		{ printf("v %e %e %e %e %e %e\n", v[0], v[1], v[2], v[3], v[4], v[5]); }
	fmatD	getTensor() const							//! returns a 3x3 tensor matrix.
		{ fvecD t(6); for (unsigned int i = 0; i < 6; i++) t[i] = v[i];
		  fmatD D(3,3); vec2mat(D,t); return D; }
	fvecD	decompose(fmatD& R) const						//! performs eigenvalue decomposition of tensor.
		{ const fmatD D = getTensor(); const auto S = sev(D);
		  R = S.U; return S.d/S.d.sum(); }
	fvec3	getColor() const							//! converts largest eigenvector to RGB.
		{ fmatD R; fvecD e = decompose(R);					// convert to tensor & decompose
		  fvec3 c(std::abs(R(0,2)), std::abs(R(1,2)), std::abs(R(2,2)));
		  return c.normalize()*e[2]; }
	fvecD	getCoefficients() const							//! converts to spharm (l=2).
		{ fvecD cf(6);								// see Descoteaux, thesis p. 81
		  cf[0] = float((2.0*std::sqrt(M_PI)/3.0)*(v[0]+v[3]+v[5]));
		  cf[1] = float((2.0*std::sqrt(M_PI/15.0))*(v[0]-v[3]));
		  cf[2] = float((4.0*std::sqrt(M_PI/15.0))*v[2]);
		  cf[3] = float((-2.0*std::sqrt(M_PI/45.0))*(v[0]+v[3]-2.0*v[5]));
		  cf[4] = float((4.0*std::sqrt(M_PI/15.0))*v[4]);
		  cf[5] = float((4.0*std::sqrt(M_PI/15.0))*v[1]);
		  return cf; }
	std::vector<fiberDesc> getFibers() const					//! returns a fiber descriptor.
		{ fmatD R; fvecD e = decompose(R);
		  fiberDesc fd1(e[2], 1.0, fvec3(R(0,2),R(1,2),R(2,2))*e[2]);
		  fiberDesc fd2(e[1], 1.0, fvec3(R(0,1),R(1,1),R(2,1))*e[1]);
		  std::vector<fiberDesc> vfd; vfd.push_back(fd1); vfd.push_back(fd2); 
		  return vfd; }
};

//! Implements a spherical harmonics decomposition record of order 4.

struct spharmDesc : public node {
	float	v[];				//!< variable number of coefficients
	void	print()									//! prints the record.
		{ printf("v %e %e %e %e %e %e\n", v[0], v[1], v[2], v[3], v[4], v[5]); }
	fvecD	getCoefficients(const unsigned int n) const				//! returns element vector.
		{ fvecD t(n); for (unsigned int i = 0; i < n; i++) t[i] = v[i];
		  return t; }
	spharmDesc();									// forbids calling default constructor
	static spharmDesc* create(const unsigned int n)					// this is ugly...but works.
		{ return reinterpret_cast<spharmDesc*>(new(n*repnTable[static_cast<int>(repnType::flt)].precision/8) node()); }
	float	generalizedAnisotropy(const unsigned int n)
		{ float c0 = SQR(v[0]), cs = 0;
		  for (unsigned int i = 0; i < n; i++) cs += SQR(v[i]);
		  return std::sqrt(1.0f-c0/cs); }
};

//! Implements a fiber mixture record.

struct fiberMixtureDesc : public node {		// 19 floats
	float	piso;				//!< isotropic fraction
	float	s0;				//!< signal intensity
	float	si2;				//!< precision
	float	nf;				//!< number of fibers
	fiberDesc fd[3];			//!< up to three fiber descriptors

	fiberMixtureDesc()								//! allocates an empty fiber mixture record.
		 : piso(0.0), s0(0.0), si2(0.0), nf(0.0) { }
	void	print()									//! prints the record.
		{ printf("s0 %e si2 %e piso %e\n", s0, si2, piso);
		  for (unsigned int i = 0; i < nf; i++) fd[i].print(); }
	std::vector<fiberDesc> getFibers() const					//! returns a vector of fiber descriptors.
		{ std::vector<fiberDesc> vfd;
		  for (unsigned int i = 0; i < nf; i++) vfd.push_back(fd[i]);
		  return vfd; }
};

//! Implements a record for the Zeppelin-Cylinder-Dot model.

struct fiberZCDDesc : public node {		// 11 floats
	float	s0;				//!< gradient-free signal
	float	si2;				//!< gradient-free signal
	float	dpar;				//!< parallel diffusivity
	float	dper;				//!< perpendicular diffusivity
	float	rad;				//!< cylinder radius
	float	mu[3];				//!< fiber direction
	float	p[3];				//!< volume fractions

	fiberZCDDesc()									//! allocates an empty ZCD record.
		 : s0(0.0), si2(0.0), dpar(0.0), dper(0.0), rad(0.0) { }
	void	print()									//! prints the record.
		{ printf("s0 %e si2 %e piso %e\n", s0, si2, p[2]);
		  printf("ic %e rad %e dpar %e mu %e %e %e\n", p[0], rad, dpar, mu[0], mu[1], mu[2]);
		  printf("ec %e dper %e\n", p[1], dper); }
	std::vector<fiberDesc> getFibers() const					//! returns a vector of fiber descriptors.
		{ std::vector<fiberDesc> vfd;
		  for (unsigned int i = 0; i < 3; i++) { if (p[i] == 0.0f) continue;
			fiberDesc fd(p[i], 20.0f, mu[i]); vfd.push_back(fd); };
		  return vfd; }
};

//! Implements a graph variant with vectors of specified length.

class vectorGraph : public graph {
protected:
	unsigned int nx;			//!< grid dimension x
	unsigned int ny;			//!< grid dimension y
	unsigned int nz;			//!< grid dimension z
public:
   	vectorGraph()									//! allocates an empty vector graph.
		 : graph(), nx(0), ny(0), nz(0) { }
   	vectorGraph(const fimage& src, unsigned int sz, const char* ip)			//! allocates a vector graph for image src, field size sz, and type *ip.
		 : graph(src.nel(), sz, repnType::flt, false), nx(src.nx), ny(src.ny), nz(src.nz)
		{ char buf[40]; const char *s;
		  sprintf(buf, "%u %u %u", nx, ny, nz); at.add(attribute("extent", buf));
		  fvec3 v = 1; src.getVoxelSize(v); sprintf(buf, "%.4f %.4f %.4f", v.x, v.y, v.z);
		  at.add(attribute("voxel", buf));
		  s = src.at.lookup("patient"); if (s) at.add(attribute("patient", s));
		  s = src.at.lookup("date"); if (s) at.add(attribute("date", s));
		  s = src.at.lookup("convention"); if (s) at.add(attribute("convention", s));
		  at.add(attribute("interpretation", ip)); }
   	vectorGraph(const graph& g) : graph(g), nx(0), ny(0), nz(0) { setExtent(); }	//! copies from graph g.
   	vectorGraph(const vectorGraph& g) : graph(g), nx(g.nx), ny(g.ny), nz(g.nz) { }	//! copies from vector graph g.
   	vectorGraph(vectorGraph&& g)
		 : graph(std::move(g)), nx(g.nx), ny(g.ny), nz(g.nz) { }		//! moves from vector graph g.
  	virtual ~vectorGraph();
	vectorGraph& operator=(const vectorGraph& g)					//! copy assigns from vector graph g.
		{ if (this != &g) { graph::operator=(static_cast<graph const&>(g));
		  	nx = g.nx; ny = g.ny; nz = g.nz; }; return *this; }
	vectorGraph& operator=(vectorGraph&& g)						//! move assigns from vector graph g.
		{ assert(this != &g); graph::operator=(static_cast<graph const&>(g));
		  nx = g.nx; ny = g.ny; nz = g.nz; return *this; }
	virtual void read(FILE *) = 0;							// provided by subclass
	void	read(const char *fname)							//! reads vector graph from file named fname.
		{ FILE *fp = openFile(fname, "r"); read(fp); fclose(fp); setExtent(); }
	void 	save(FILE *fp)								//! saves vector graph in file* fp. 
		{ graph::save(fp); }
	void	save(const char *fname)							//! saves vector graph in file named fname.
		{ FILE *fp = openFile(fname, "w"); graph::save(fp); fclose(fp); }
	void	addNode(const node& n, const uvec3& s)					//! adds node n at site s.
		{ unsigned int i = index(s); delete table[i]; table[i] = copyNode(n); }
	node*	operator()(const uvec3& s) const					//! returns node at site s.
		{ unsigned int i = index(s); return i <= size? table[i]: 0; }
	node*	operator()(const unsigned int i) const					//! returns node at index i.
		{ return i <= size? table[i]: 0; }
	uvec3	extent() const								//! returns extent of underlying image.
		{ return uvec3(nx,ny,nz); }
	void	setExtent()								//! sets extent attribute.
		{ const char *s = at.lookup("extent");
		  if (s) sscanf(s, "%u %u %u", &nx, &ny, &nz);
		  else throw rtException("Extent in vectorGraph missing"); }
	unsigned int nel() const							//! returns number of sites in vector graph.
		{ return nx*ny*nz; }
	unsigned int index(const uvec3& s) const					//! returns index of site s.
		{ return (s.z*ny+s.y)*nx+s.x; }
	uvec3	site(const unsigned int i) const					//! returns site of index i.
		{ unsigned int z = i/(nx*ny), r = i%(nx*ny); return uvec3(r%nx,r/nx,z); }
	bool	getVoxelSize(fvec3 &v) const						//! returns voxel size in read world dimensions.
		{ const char *s = at.lookup("voxel");
		  if (s) sscanf(s, "%f%f%f", &v.x, &v.y, &v.z); return s != 0; }
	unsigned int nVertices() const							//! returns index of last node in use.
		{ return last; }
	virtual float intensityAt(const uvec3&) const					//! returns any sensible intensity at s.
		{ return 0.0; }
	virtual fvec3 colorAt(const uvec3&) const					//! returns any sensible RGB value at s.
		{ return fvec3(0.0); }
	virtual fmatD tensorAt(const uvec3&) const					//! returns any sensible 3x3 tensor at s.
		{ return fmatD(); }
	virtual fvecD spharmCoeffsAt(const uvec3&) const				//! returns spharm coefficients at s.
		{ return fvecD(); }
	virtual renderMethod defaultMode()						//! display hint.
		{ return renderMethod::off; }
	virtual std::vector<fiberDesc> fibersAt(const uvec3&) const			//! returns any sensible fiber descriptors at s.
		{ return std::vector<fiberDesc>(); }
};

//! Implements a graph variant where nodes contain fiber segments.

class streamlineGraph : public vectorGraph {
public:
   	streamlineGraph()								//! allocates an empty streamline graph.
		 : vectorGraph() { useWeights = true; }
   	streamlineGraph(const fimage& src)						//! allocates a streamline graph for image src.
		 : vectorGraph(src, 0, "streamline") { useWeights = true; }
   	streamlineGraph(const graph& g)							//! copies from graph g.
		 : vectorGraph(g) { useWeights = true; }
   	streamlineGraph(FILE* fp)							//! reads from FILE* fp.
		 : vectorGraph() { read(fp); setExtent(); }
  	virtual ~streamlineGraph();
	void	read(FILE *fp)								//! reads streamline graph from file* fp.
		{ bundleList list; readFile(fp, list); graph::read(*list.begin());
		  const char *s = at.lookup("interpretation");		// check for graph type
		  if (strncmp(s, "streamline", 10)) throw rtException("illegal input for streamlineGraph"); 
		  setExtent(); }
	bool	addSegment(const uvec3& s, const unsigned int j, const unsigned int lc)	//! adds segment from lc to site s.
		{ node n; unsigned int i = index(s); n.weight = j? 0.0f: float(lc);	// indicate line start
		  if (addNodeAt(n,i) == false) return false;
		  return j? linkNodesUni(j,i,float(lc)): true; }
	void	resize(unsigned int n)							//! resizes graph to n fields.
		{ graph::resize(n, 0, repnType::flt, true); }
	renderMethod defaultMode()							//! display hint.
		{ return renderMethod::spline; }
};

//! Implements a graph variant where nodes represent tensors.

class tensorGraph : public vectorGraph {
public:
   	tensorGraph()									//! allocates an empty tensor graph.
		 : vectorGraph() { }
   	tensorGraph(const fimage& src)							//! allocates a tensor graph for image src.
		 : vectorGraph(src, 6, "tensor") { }
   	tensorGraph(const graph& g)							//! copies from graph g.
		 : vectorGraph(g) { }
   	tensorGraph(FILE* fp)								//! reads from FILE* fp.
		 : vectorGraph() { read(fp); setExtent(); }
  	virtual ~tensorGraph();
	void	read(FILE *fp)								//! reads tensor graph from file* fp.
		{ bundleList list; readFile(fp, list); graph::read(*list.begin());
		  const char *s = at.lookup("interpretation");				// check for graph type
		  if (strncmp(s, "tensor", 5)) throw rtException("illegal input for tensorGraph"); 
		  setExtent(); }
	renderMethod defaultMode()							//! display hint.
		{ return renderMethod::tensor; }
	fvec3	colorAt(const uvec3& s) const						//! returns RGB code from largest eigenvector.
		{ const tensorDesc* td = reinterpret_cast<const tensorDesc*>((*this)(s));
		  return td? td->getColor(): fvec3(0.0); }
	fmatD	tensorAt(const uvec3& s) const						//! returns 3x3 tensor at s.
		{ const tensorDesc* td = reinterpret_cast<const tensorDesc*>((*this)(s));
		  return td? td->getTensor(): fmatD(); }
	fvecD	spharmCoeffsAt(const uvec3& s) const					//! returns spharm coefficients at s.
		{ const tensorDesc* td = reinterpret_cast<const tensorDesc*>((*this)(s));
		  return td? td->getCoefficients(): fvecD(); }
	std::vector<fiberDesc> fibersAt(const uvec3& s) const				//! returns fiber descriptors at s.
		{ const tensorDesc* td = reinterpret_cast<const tensorDesc*>((*this)(s));
		  return td? td->getFibers(): std::vector<fiberDesc>(); }
};

//! Implements a graph variant where nodes represent a spherical harmonics parametrization.

class spharmGraph : public vectorGraph {
public:
   	spharmGraph()									//! allocates an empty spharm graph.
		 : vectorGraph() { }
   	spharmGraph(const fimage& src, const unsigned int nfields)			//! allocates a spharm graph for image src.
		 : vectorGraph(src, nfields, "spharm") { }
   	spharmGraph(const graph& g)							//! copies from graph g.
		 : vectorGraph(g) { }
   	spharmGraph(FILE* fp)								//! reads from FILE* fp.
		 : vectorGraph() { read(fp); setExtent(); }
  	virtual ~spharmGraph();
	void	read(FILE *fp)								//! reads spharm graph from file* fp.
		{ bundleList list; readFile(fp, list); graph::read(*list.begin());
		  const char *s = at.lookup("interpretation");				// check for graph type
		  if (strncmp(s, "spharm", 5)) throw rtException("illegal input for spharmGraph"); }
	renderMethod defaultMode()							//! display hint.
		{ return renderMethod::spharm; }
	fvecD	spharmCoeffsAt(const uvec3& s) const					//! returns spharm coefficients at s.
		{ const spharmDesc* sd = reinterpret_cast<const spharmDesc*>((*this)(s));
		  return sd? sd->getCoefficients(nfields): fvecD(); }
	unsigned int degree() const							//! returns order of SH decomposition with nc coefficients.
		{ return FTOU(std::sqrt(8*nfields+1)-3)/2; }
};

//! Implements a graph variant where nodes contain a fiber mixture.

class fiberMixtureGraph : public vectorGraph {
public:
   	fiberMixtureGraph()								//! allocates an empty fiber mixture graph.
		 : vectorGraph() { }
   	fiberMixtureGraph(const fimage& src)						//! allocates a fiber mixture graph for image src.
		 : vectorGraph(src, 19, "fiberMixture") { }
   	fiberMixtureGraph(const graph& g)						//! copies from graph g.
		 : vectorGraph(g) { }
   	fiberMixtureGraph(FILE* fp)							//! reads from FILE* fp.
		 : vectorGraph() { read(fp); setExtent(); }
  	virtual ~fiberMixtureGraph();
	void	read(FILE *fp)								//! reads fiber mixture graph from file* fp.
		{ bundleList list; readFile(fp, list); graph::read(*list.begin());
		  const char *s = at.lookup("interpretation");				// check for graph type
		  if (strncmp(s, "fiberMixture", 12)) throw rtException("illegal input for fiberMixtureGraph"); }
	renderMethod defaultMode()							//! display hint.
		{ return renderMethod::cylinder; }
	std::vector<fiberDesc> fibersAt(const uvec3& s) const				//! returns fiber descriptors at s.
		{ const fiberMixtureDesc* fm = reinterpret_cast<const fiberMixtureDesc*>((*this)(s));
		  return fm? fm->getFibers(): std::vector<fiberDesc>(); }
};

//! Implements a graph variant where nodes contain Zeppelin-Cylinder-Dot model parameters.

class fiberZCDGraph : public vectorGraph {
public:
   	fiberZCDGraph()									//! allocates an empty fiber ZCD graph.
		 : vectorGraph() { }
   	fiberZCDGraph(const fimage& src)						//! allocates a fiber ZCD graph for image src.
		 : vectorGraph(src, 12, "fiberZCD") { }
   	fiberZCDGraph(const graph& g)							//! copies from graph g.
		 : vectorGraph(g) { }
   	fiberZCDGraph(FILE* fp)								//! reads from FILE* fp.
		 : vectorGraph() { read(fp); setExtent(); }
  	virtual ~fiberZCDGraph();
	void	read(FILE *fp)								//! reads fiber ZCDe graph from file* fp.
		{ bundleList list; readFile(fp, list); graph::read(*list.begin());
		  const char *s = at.lookup("interpretation");				// check for graph type
		  if (strncmp(s, "fiberZCD", 8)) throw rtException("illegal input for fiberZCDGraph"); }
	renderMethod defaultMode()							//! display hint.
		{ return renderMethod::cylinder; }
	std::vector<fiberDesc> fibersAt(const uvec3& s) const				//! returns fiber descriptors at s.
		{ const fiberZCDDesc* fm = reinterpret_cast<fiberZCDDesc*>((*this)(s));
		  return fm? fm->getFibers(): std::vector<fiberDesc>(); }
};

#endif
