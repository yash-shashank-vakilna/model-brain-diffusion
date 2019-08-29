#ifndef SPLINE_H
#define SPLINE_H

/*
 *
 * spline.h: splines in 1d and 3d
 * BRIAN Software Package Version 3.0
 *
 * $Id: spline.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (16/05/10): initial version
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.50 (28/08/15): reimplemented based on fsl_splines.h (Copyright (C) 2007 University of Oxford) 
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for splines in 1d and 3d.
*/

//! Implements functions for third order spline interpolation in one dimension.

template<typename T> class spline1d {
	unsigned int sp;			//!< knot spacing.
	unsigned int dv;			//!< degree of derivative.
	vecD<T> cf;  				//!< vector of coefficients.

	T	valueAt(const T x) const						//! returns interpolated value at x in [-2,2] for spline derivative dv.
		{ const T ax = std::abs(x);
		  switch (dv) {
		  case 0: if (ax <= T(1.0)) return T(2.0/3.0)+SQR(ax)*(T(0.5)*ax-T(1.0));
			  else if (ax < T(2.0)) return T(1.0/6.0)*SQR(T(2.0)-ax)*(T(2.0)-ax);
			  else return T(0);
		  case 1: if (ax < T(1e-6)) return 0;
			  else if (ax <= T(1.0)) return (x/ax)*(T(1.5)*SQR(ax)-T(2.0)*ax)/sp;
			  else if (ax < T(2.0)) return (x/ax)*(T(-0.5)*(T(2.0)-ax)*(T(2.0)-ax))/sp;
			  else return T(0);
		  case 2: if (ax <= T(1.0)) return (T(3.0)*ax-T(2.0))/SQR(sp);
			  else if (ax < T(2.0)) return (T(2.0)-ax)/SQR(sp);
			  else return T(0);
		  default: throw rtException("spline1d: only differentiable three times"); } }
	vecD<T>	coefficients() const							//! returns vector of spline values.
		{ vecD<T> p(4*sp-1); const T c = T(0.5)*(p.N-1);			// center of array in sites
		  for (unsigned int i = 0; i < p.N; i++) p[i] = valueAt((T(i)-c)/sp);
		  return p; }
 	int	index(const unsigned int k) const					//! returns the voxel index of knot k.
 		{ return int(k*sp)-int(3*sp/2); }
public:
	spline1d(const unsigned int _sp = 1, const unsigned int _dv = 0)
		 : sp(_sp), dv(_dv), cf(coefficients()) 
		{ if (sp < 1) throw rtException("spline1d: knot-spacing has to be positive"); }
	unsigned int ksp() const							//! returns knot spacing in sites.
		{ return sp; }
	unsigned int nc() const								//! returns number of spline coefficients.
		{ return cf.N; }
	unsigned int nk(const unsigned int nx) const					//! returns number of knots that represent data of length nx.
		{ unsigned int k = nx/sp+1;
		  while (index(k)-cf.N/2 < nx-1) k++;
		  return k; }
	void	range(const unsigned int k, const unsigned int nx,
			unsigned int& f, unsigned int& l) const				//! returns site range (f,l) that knot k supports.
		{ f = std::max(0,index(k)-int(cf.N/2));
		  l = std::min(int(nx),index(k)+int(cf.N/2)+1); }
	unsigned int coefBegin(const unsigned int k) const				//! returns the array index of the first coefficient of knot k.
		{ unsigned int i = std::max(0,int(cf.N/2)-index(k));
		  return std::min(i,cf.N); }
	void	knotRange(const T x, const unsigned int cx,
			unsigned int& f, unsigned int& l) const				//! returns knot range (f,l) that contributes to position x.
		{ const T dx = x-index(0), c2 = T(0.5)*(cf.N+1);
		  f = std::max(0,int(std::ceil((dx-c2)/sp)));
		  l = std::min(int(cx),std::max(0,int(std::ceil((dx+c2)/sp)))); }
	void	coefRange(const unsigned int k, const unsigned int nx,
			unsigned int& f, unsigned int& l) const				//! returns coefficient range (f,l) that contributes to knot k.
		{ f = k < 4? 0: k-3; l = std::min(nk(nx),k+4); }
	void	coefRange(const unsigned int k1, const unsigned int nx,
			const spline1d& s, unsigned int& f, unsigned int& l) const	//! returns coefficient range (f,l) that overlaps with spline s.	
		{ unsigned int f1, l1; range(k1,nx,f1,l1); int k0 = -1; 
		  for (unsigned int k2 = 0; k2 < s.nk(nx); k2++) {
			unsigned int f2, l2; s.range(k2,nx,f2,l2);
			const int ov = std::min(int(l1)-1,int(l2)-1)-std::max(int(f1),int(f2));
			if (ov >= 0) { if (k0 == -1) k0 = k2; else l = k2+1; }
			else { if (k0 != -1) break; } };
		  if (k0 == -1) { f = 0; l = 0; } else f = k0; }
	T	operator()(const unsigned int i) const					//! returns spline coefficient i.
		{ return cf(i); }
	T	operator()(const T x, const unsigned int k) const			//! interpolates spline at position x, relative to knot k.
		{ return valueAt((x-index(k))/sp); }
	void	overlapsIn(const unsigned int k1, const unsigned int nx,
			const spline1d& s, const unsigned int k2,
			unsigned int& f1, unsigned int& f2, unsigned int& l) const	//! returns coefficient range (f1,f2,l) that overlap between this spline and s.
		{ unsigned int o1 =   coefBegin(k1), if1, il1;   range(k1,nx,if1,il1);	// coefficient offset and range of sites (image pixels) in this spline
		  unsigned int o2 = s.coefBegin(k2), if2, il2; s.range(k2,nx,if2,il2);	// coefficient offset and range of sites (image pixels) in spline s
		  const unsigned int d1 = std::max(0,int(if2)-int(if1)); f1 = std::min(o1+d1,cf.N);
		  const unsigned int d2 = std::max(0,int(if1)-int(if2)); f2 = std::min(o2+d2,s.cf.N);
		  const unsigned int len = std::max(0,int(std::min(il1,il2))-int(if1));
		  l = std::min(o1+len,cf.N); }
	matD<T>	reconMatrix(const unsigned int nx, const unsigned int nk) const		//! returns reconstruction matrix: A*coef = data.
		{ matD<T> A(nx,nk); A = T(0);
		  for (unsigned int k = 0; k < nk; k++) {
			unsigned int s = coefBegin(k), f, l; range(k,nx,f,l);
				for (unsigned int i = f; i < l; i++) A(i,k) = cf(s++); };
		  return A; }
};

//! Implements functions for third order spline interpolation in three dimensions.

template<typename T> class spline3d {
	spline1d<T> sx, sy, sz;			//!< 1d splines in x,y,z direction
	image<T> cf;				//!< 3d spline coefficients
	image<T> cfm;				//!< 3d spline coefficients, possibly pre-multiplied
public:
	spline3d()
		: sx(), sy(), sz(), cf(), cfm()
		{ };
	spline3d(const uvec3& sp, const uvec3& dv = 0)					//! allocates a 3d spline with knot spacing sp and derivative dv.
		 : sx(spline1d<T>(sp.x,dv.x)), sy(spline1d<T>(sp.y,dv.y)), sz(spline1d<T>(sp.z,dv.z)),
		cf(), cfm()
		{ uvec3 ex(sx.nc(),sy.nc(),sz.nc()); cf = image<T>(ex); 
		  allSites(s,0) cf(s) = sx(s.x)*sy(s.y)*sz(s.z);
		  cfm = cf; }
	uvec3	ksp() const								//! returns knot spacing in voxels.
		{ return uvec3(sx.ksp(),sy.ksp(),sz.ksp()); }
	uvec3	nk(const uvec3& ex) const						//! returns the extent of knots that specify a 3d spline in (0,ex).
		{ return uvec3(sx.nk(ex.x),sy.nk(ex.y),sz.nk(ex.z)); }
	unsigned int nnz(const uvec3& ex) const						//! returns number of non-zero coefficients.
		{ unsigned int n = 0, f = 0, l = 0; const uvec3 cx = nk(ex);
		  for (unsigned int z = 0; z < cx.z; z++) {
			sz.coefRange(z,ex.z,f,l); const unsigned int zr = l-f;
			for (unsigned int y = 0; y < cx.y; y++) {
				sy.coefRange(y,ex.y,f,l); const unsigned int yr = l-f;
				for (unsigned int x = 0; x < cx.x; x++) {
					sx.coefRange(x,ex.x,f,l); n += zr*yr*(l-f); } } };
		  return n; }
	unsigned int nnz(const uvec3& ex, const spline3d& s) const			//! returns number of non-zero coefficients overlapping with 3d spline s.
		{ unsigned int n = 0, f = 0, l = 0; const uvec3 cx = nk(ex);
		  for (unsigned int z = 0; z < cx.z; z++) {
			sz.coefRange(z,ex.z,s.sz,f,l); const unsigned int zr = l-f;
			for (unsigned int y = 0; y < cx.y; y++) {
				sy.coefRange(y,ex.y,s.sy,f,l); const unsigned int yr = l-f;
				for (unsigned int x = 0; x < cx.x; x++) {
					sx.coefRange(x,ex.x,s.sx,f,l); n += zr*yr*(l-f); } } };
		  return n; }
	T	operator()(const unsigned int x, const unsigned int y,
			const unsigned int z) const					//! returns (pre-muliplied) spline coefficient (x,y,z).
		{ return cfm(x,y,z); }
	T	operator()(const vec3<T>& p, const unsigned int x, const unsigned int y,
			const unsigned int z) const					//! interpolates spline at p, relative to knot t.
		{ return sx(p.x,x)*sy(p.y,y)*sz(p.z,z); }
	spline3d& operator*=(const T b)							//! multiplies this spline by scalar b.
		{ cfm *= b; return *this; }
	spline3d& operator/=(const T b)							//! divides this spline by scalar b.
		{ cfm /= b; return *this; }
	void	range(const uvec3& t, const uvec3& ex, uvec3& f, uvec3& l) const	//! returns the range of data for which site t has support.
		{ sx.range(t.x,ex.x,f.x,l.x); sy.range(t.y,ex.y,f.y,l.y);
		  sz.range(t.z,ex.z,f.z,l.z); }         
	void	knotRange(const vec3<T>& t, const uvec3& cx, uvec3& f, uvec3& l) const	//! returns the range of data for which site t has support.
		{ sx.knotRange(t.x,cx.x,f.x,l.x); sy.knotRange(t.y,cx.y,f.y,l.y);
		  sz.knotRange(t.z,cx.z,f.z,l.z); }         
	void	coefRange(const uvec3& t, const uvec3& ex, uvec3& f, uvec3& l) const 	//! returns the range of coefficients that contribute to site t. 
		{ sx.coefRange(t.x,ex.x,f.x,l.x); sy.coefRange(t.y,ex.y,f.y,l.y);
		  sz.coefRange(t.z,ex.z,f.z,l.z); }
	T	multiplyBy(const uvec3& t1, const uvec3& t2, const uvec3& ex, const spline3d& s3) const
		{ uvec3 f1, f2, l; T v = T(0);
		  sx.overlapsIn(t1.x,ex.x,s3.sx,t2.x,f1.x,f2.x,l.x);
		  sy.overlapsIn(t1.y,ex.y,s3.sy,t2.y,f1.y,f2.y,l.y);
		  sz.overlapsIn(t1.z,ex.z,s3.sz,t2.z,f1.z,f2.z,l.z);
		  for (unsigned int z1 = f1.z, z2 = f2.z; z1 < l.z; z1++, z2++)
			for (unsigned int y1 = f1.y, y2 = f2.y; y1 < l.y; y1++, y2++)
				for (unsigned int x1 = f1.x, x2 = f2.x; x1 < l.x; x1++, x2++)
					v += cfm(x1,y1,z1)*s3.cfm(x2,y2,z2);
		  return v; }
	uvec3	coefBegin(const uvec3& s) const						//! returns offset in (x,y,z) direction.
		{ return uvec3(sx.coefBegin(s.x),sy.coefBegin(s.y),sz.coefBegin(s.z)); } 
	void	premul(const uvec3& t, const uvec3& ex, const fimage& src)		//! multiply coefficient array by float image src.
		{ uvec3 c = coefBegin(t), f, l; range(t,ex,f,l); cfm = T(0);
		  for (unsigned int z = f.z, cz = c.z; z < l.z; z++, cz++)
		  	for (unsigned int y = f.y, cy = c.y; y < l.y; y++, cy++)
		  		for (unsigned int x = f.x, cx = c.x; x < l.x; x++, cx++) {
					const unsigned int i = cf.index(cx,cy,cz);
					cfm(i) = cf(i)*src(x,y,z); } }
};

//! Implements functions for spline interpolation of three-dimensional images.
template<typename T> class splineField {

	//! Implements a helper class for energy computation.
	class enHelper {
		image<T> im;			//!< energy image.
		uvec3	c;			//!< center of coefficient space.
	public:
		enHelper(const spline3d<T>& sp)						//! allocates a structure for energy computation on spline sp.
			: im(), c()
			{ im = T(0); const uvec3 ex = sp.ksp()*1000u, s = sp.nk(ex)/2;	// fake a "really large" FOV
			  uvec3 f,l; sp.coefRange(s,ex,f,l);				// get indices of overlapping splines
			  uvec3 dx = l-f; c = s-f; im = image<T>(dx);			// initialize this structure
			  for (unsigned int z = 0; z < dx.z; z++)			// get values for all overlaps
				for (unsigned int y = 0; y < dx.y; y++)
					for (unsigned int x = 0; x < dx.x; x++)
						im(x,y,z) = sp.multiplyBy(s,f+uvec3(x,y,z),ex,sp); }
		T	operator()(const int x, const int y, const int z) const		//! returns coefficient (x,y,z).
			{ return im(x+c.x,y+c.y,z+c.z); }
		enHelper& operator*=(const T b)						//! multiplies coefficients by scalar b.
			{ im *= b; return *this; }
		enHelper& operator+=(const enHelper& h)					//! adds coefficients of another helper to this one.
			{ im += h.im; return *this; }
	};

	//! Implements a helper class for extrapolating spline coefficients.
	class volHelper {								//! Implements a helper class for extending spline coefficients.
		image<T> im;			//!< coefficient image.

		unsigned int len(const unsigned int dim) const				//! returns length of coefficient image in dimension dim.
			{ switch (dim) {
			  case 0: return im.nx;
			  case 1: return im.ny;
			  case 2: return im.nz;
			  default: throw optException("volHelper::step: Unknown dimension"); } }
		unsigned int start(const unsigned int x, const unsigned int y,
			const unsigned int dim) const					//! returns start index into coefficient image (x,y) in dimension dim.
			{ switch (dim) {
			  case 0: return (y*im.ny+x)*im.nx;
			  case 1: return y*im.ny*im.nx+x;
			  case 2: return y*im.nx+x;
			  default: throw optException("volHelper::step: Unknown dimension"); } }
		unsigned int end(const unsigned int x, const unsigned int y,
			const unsigned int dim) const					//! returns end index into coefficient image (x,y) in dimension dim.
			{ const unsigned int l = len(dim);
			  return start(x,y,dim)+l*step(dim); }
		unsigned int step(const unsigned dim) const				//! returns index increment for coefficient image in dimension dim.
			{ switch (dim) {
			  case 0: return 1;
			  case 1: return im.nx;
			  case 2: return im.nx*im.ny;
			  default: throw optException("volHelper::step: Unknown dimension"); } }
		vecD<T>	getColumn(const unsigned int x, const unsigned int y,
				const unsigned int dim) const				//! returns coefficient column at (x,y) in dimension dim.
			{ vecD<T> c(len(dim)); c = T(0); unsigned int t = 0;
			  const unsigned int e = end(x,y,dim), inc = step(dim);
			  for (unsigned int i = start(x,y,dim); i < e; i += inc) c[t++] = im(i);
			  return c; }
		void	setColumn(const unsigned int x, const unsigned int y,
				const unsigned int dim, const vecD<T>& c)		//! sets coefficient column c at (x,y) in dimension dim.
			{ const unsigned int e = end(x,y,dim), inc = step(dim); unsigned int j = 0;
			  for (unsigned int i = start(x,y,dim); i < e; i += inc) im(i) = c(j++); }
		matD<T>	getSmoothing(const unsigned int l, const unsigned int n) const	//! regularization for extending coefficients
			{ matD<T> S(l,n); S = T(0);
			  S(0,0) = T(2.0); S(0,1) = T(-1.0); S(0,n-1) = T(-1.0);
			  S(l-1,0) = T(-1.0); S(l-1,n-2) = T(-1.0); S(l-1,n-1) = T(2.0);
			  return S; }
	public:
		volHelper(const image<T>& src)						//! allocates a structure for extending spline coefficients from coefficient image src.
			: im(src) { }
		volHelper(const uvec3& ex)						//! allocates a structure for extending spline coefficients for extent ex.
			: im(image<T>(ex)) { im = T(0); }
		volHelper getCoef(const unsigned int dim, const unsigned int n,
				const unsigned int ksp)	const				//! extends coefficients in dimension dim to n elements, given knot spacing ksp.
			{ uvec3 ex = im.extent(); unsigned int dx, dy, l = len(dim);
			  switch (dim) {
			  case 0: ex.x = n; dx = im.ny; dy = im.nz; break;
			  case 1: ex.y = n; dx = im.nx; dy = im.nz; break;
			  case 2: ex.z = n; dx = im.nx; dy = im.ny; break;
			  default: throw optException("getCoef: Unknown direction"); }
			  volHelper dst(ex); spline1d<T> sp(ksp);
			  const matD<T> A = sp.reconMatrix(l,n);
			  matD<T> AS = A; AS.joinBelow(T(0.005)*getSmoothing(l,n));
			  const matD<T> AStAS = trp(AS)*AS, M = inv(AStAS)*trp(A);
			  for (unsigned int x = 0; x < dx; x++)
				for (unsigned int y = 0; y < dy; y++) {
					const vecD<T> c = getColumn(x,y,dim);
					dst.setColumn(x,y,dim,M*c); };
			  return dst;  }
		image<T> cf() const							//! return coefficient image.
			{ return im; }
	};

	uvec3	ex;			//!< image extent.
	vec3<T>	vs;			//!< image voxel size.
	spline3d<T> s3;			//!< spline interpolant in three dimensions.
	image<T> cf;			//!< volume of spline coefficients.
	std::vector<image<T>> fd;	//!< spline-interpolated image and derivatives.

	T	getVal(const matS<T>& M, const unsigned int r, const unsigned int c) const
		{ const unsigned int t = M.ri[c], n = M.ri[c+1]-t, *a = M.ci+t; 
		  assert(n > 0); if (r < a[0] || r > a[n-1]) return T(0);
		  int i, lo = -1, up = n;  while (up-lo > 1) {
			i = (lo+up) >> 1; if (r >= a[i]) lo = i; else up = i; };
		  return a[lo] == r? M.x[t+lo]: T(0); }
	void	updateFields()								//! updates coefficient fields from coefficient vector cf.
		{ fd[0] = makeField(spline3d<T>(s3.ksp()));
		  fd[1] = makeField(spline3d<T>(s3.ksp(),uvec3(1,0,0)));
		  fd[2] = makeField(spline3d<T>(s3.ksp(),uvec3(0,1,0)));
		  fd[3] = makeField(spline3d<T>(s3.ksp(),uvec3(0,0,1))); }
	image<T> makeField(const spline3d<T>& sp) const;
	vecD<T>	getJte(const spline3d<T>& sp, const fimage& src) const;
	matS<T>	getJtJ(const spline3d<T>& ds1, const spline3d<T>& ds2, const fimage& src) const;
	matS<T> symJtJ(const spline3d<T>& sp, const fimage& src) const;
public:
	splineField()
		: ex(), vs(), s3(), cf(), fd(4)
		{ }
	splineField(const uvec3& _ex, const vec3<T>& _vs, const uvec3& _sp)		//! allocates a context for spline interpolation of an image of extent _ex, voxel size _vs and knot spacing _sp.
		: ex(_ex), vs(_vs), s3(_sp), cf(s3.nk(ex)), fd(4)
		{ cf = T(0); }
	uvec3	extent() const								//! returns image extent.
		{ return ex; }
	vec3<T>	getVoxelSize() const							//! returns voxel size.
		{ return vs; }
	void	setCoefficients(const image<T>& _cf)					//! updates coefficients from vector _cf and updates fields.
		{ cf = _cf; updateFields(); }
	image<T> getCoefficients() const						//! retrieves coefficient image cf.
		{ return cf; }
	T	interpolateAt(const vec3<T>& p) const					//! returns image interpolated at position p.
		{ T v = T(0); if (cf.nx == 0) return v;
		  uvec3 f, l; s3.knotRange(p,cf.extent(),f,l);
		  for (unsigned int z = f.z; z < l.z; z++)
			for (unsigned int y = f.y; y < l.y; y++)
            			for (unsigned int x = f.x; x < l.x; x++)
					v += cf(x,y,z)*s3(p,x,y,z);
		  return v; }
	void	updateCoef(const image<T>& src)						//! extends coefficient images based on coefficient image src.
		{ if (ex != src.extent()) throw optException("splineField: spf mismatch");
		  const volHelper h(src); const uvec3 sp = s3.ksp();
		  volHelper c = h.getCoef(0,cf.nx,sp.x);
		  c = c.getCoef(1,cf.ny,sp.y); c = c.getCoef(2,cf.nz,sp.z);
		  cf = c.cf(); updateFields(); }
	T	operator()(const unsigned int f, const unsigned int i) const		//! returns value of field f at index i.
		{ return fd[f](i); }
	void	scale(const T b)							//! scales coefficients by scalar b and updates fields.
		{ cf *= b; updateFields(); }
	vecD<T>	Jte(const fimage& src, const uvec3& dv = 0) const			//! returns gradient of bending energy for image src, masked by mask.
		{ return dv == uvec3(0)? getJte(s3,src): getJte(spline3d<T>(s3.ksp(),dv),src); }
	matS<T>	JtJ(const fimage& src, const uvec3& dv = 0) const			//! returns symmetric Hessian of bending energy for image src.
		{ return dv == uvec3(0)? symJtJ(s3,src): symJtJ(spline3d<T>(s3.ksp(),dv),src); }
	matS<T>	JtJ(const fimage& src, const uvec3& dv1, const uvec3& dv2) const	//! returns asymmetric Hessian of bending energy for image src.
		{ return getJtJ(spline3d<T>(s3.ksp(),dv1),spline3d<T>(s3.ksp(),dv2),src); }
	T	energy() const								//! returns bending energy of this spline field.
		{ const vecD<T> g = energyGrad(); return dot(cf.x,g.x,g.N)*T(0.5); }
	vecD<T>	energyGrad() const;							//! returns gradient of bending energy of this spline field.
	matS<T>	energyHess() const;							//! returns Hessian of bending energy of this spline field.
	fimage	save() const;
	void	read(const fimage& src);
	image<T> getField(const unsigned int f) const					//! retrieves field f.
		{ return fd[f]; }
};

template <typename T>
image<T> splineField<T>::makeField(const spline3d<T>& sp) const
//! allocates spline coefficient field.
{	image<T> fld(ex); fld = T(0);
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) { const T cs = cf(s); 
		uvec3 c = sp.coefBegin(s), f, l; sp.range(s,ex,f,l);
		for (unsigned int z = f.z, cz = c.z; z < l.z; z++, cz++)
			for (unsigned int y = f.y, cy = c.y; y < l.y; y++, cy++)
            			for (unsigned int x = f.x, cx = c.x; x < l.x; x++, cx++)
					fld(x,y,z) += sp(cx,cy,cz)*cs; };
	return fld;
}

template <typename T>
vecD<T> splineField<T>::getJte(const spline3d<T>& sp, const fimage& src) const
//! returns gradient of bending energy of spline field.
{	vecD<T> g(cf.nel()); g = T(0); unsigned int t = 0;
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) { T v = T(0);
		uvec3 c = sp.coefBegin(s), f, l; sp.range(s,ex,f,l);
		for (unsigned int z = f.z, cz = c.z; z < l.z; z++, cz++)
			for (unsigned int y = f.y, cy = c.y; y < l.y; y++, cy++)
            			for (unsigned int x = f.x, cx = c.x; x < l.x; x++, cx++)
					v += sp(cx,cy,cz)*src(x,y,z);
		g[t++] = v; };
	return g;
}

template <typename T>
matS<T> splineField<T>::getJtJ(const spline3d<T>& rsp, const spline3d<T>& sp, const fimage& src) const
//! returns Hessian of bending energy of spline field.
{	spline3d<T> csp(sp); matS<T> M(cf.nel(),cf.nel(),rsp.nnz(ex,csp));
	unsigned int t = 0; M.ri[0] = t;	
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) { 
		csp.premul(s,ex,src); uvec3 f,l; csp.coefRange(s,ex,f,l);
		for (unsigned int z = f.z; z < l.z; z++)
			for (unsigned int y = f.y; y < l.y; y++)
				for (unsigned int x = f.x; x < l.x; x++) {
					M.ci[t] = cf.index(x,y,z);
					M.x[t] = csp.multiplyBy(s,uvec3(x,y,z),ex,rsp); t++; };
		M.ri[cf.index(s)+1] = t; };
	return M;									// N.B. is actually transpose
}

template <typename T>
matS<T> splineField<T>::symJtJ(const spline3d<T>& sp, const fimage& src) const
//! returns Hessian of bending energy of spline field.
{	spline3d<T> csp(sp); matS<T> M(cf.nel(),cf.nel(),sp.nnz(ex,csp)); unsigned int t = 0;;	
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) {
		const unsigned int c = cf.index(s); M.ri[c] = t;
		csp.premul(s,ex,src); uvec3 f,l; csp.coefRange(s,ex,f,l);
		for (unsigned int z = f.z; z < s.z; z++) {				// use first-level (diagonal) symmetry
			for (unsigned int y = f.y; y < l.y; y++) {
				for (unsigned int x = f.x; x < l.x; x++) {
					const unsigned int r = cf.index(x,y,z);
					M.ci[t] = r; M.x[t++] = getVal(M,c,r); } } };
		for (unsigned int z = s.z; z < l.z; z++) {
			for (unsigned int y = f.y; y < s.y; y++)			// use second-level symmetry
				for (unsigned int x = f.x; x < l.x; x++) {
					M.ci[t] = cf.index(x,y,z);
					M.x[t++] = getVal(M,cf.index(s.x,s.y,z),cf.index(x,y,s.z)); }
			for (unsigned int y = s.y; y < l.y; y++) {
				for (unsigned int x = f.x; x < s.x; x++) {		// use second-level symmetry
					M.ci[t] = cf.index(x,y,z);
					M.x[t++] = getVal(M,cf.index(s.x,y,z),cf.index(x,s.y,s.z)); }
				for (unsigned int x = s.x; x < l.x; x++) { 		// these must be computed
					M.ci[t] = cf.index(x,y,z);
 					M.x[t++] = csp.multiplyBy(s,uvec3(x,y,z),ex,sp); } } } };
	M.ri[M.M] = t; return M;
}

template <typename T>
vecD<T> splineField<T>::energyGrad() const
//! returns gradient of bending energy of this spline field.
{	const spline3d<T> sp(s3.ksp()); enHelper sum(sp); sum *= T(0);
	vecD<T> vxs(3); vxs(0) = vs.x; vxs(1) = vs.y; vxs(2) = vs.z;
	for (unsigned int d1 = 0; d1 < 3; d1++) {
		for (unsigned int d2 = d1; d2 < 3; d2++) {
			uvecD dv(3); dv = 0; dv[d1]++; dv[d2]++;
			spline3d<T> ds(s3.ksp(),uvec3(dv(0),dv(1),dv(2)));
			ds /= (vxs(d1)*vxs(d2));					// derivative in mm^{-1}
			enHelper h(ds); if (d1 != d2) h *= T(2.0); sum += h; } };
	vecD<T> g(cf.nel()); g = T(0); unsigned int t = 0;           
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) {
		uvec3 f, l; sp.coefRange(s,ex,f,l); T v = T(0);
		int fx = int(f.x)-int(s.x), fy = int(f.y)-int(s.y), fz = int(f.z)-int(s.z);
		int lx = int(l.x)-int(s.x), ly = int(l.y)-int(s.y), lz = int(l.z)-int(s.z);
		for (int z = fz; z < lz; z++)
			for (int y = fy; y < ly; y++)
				for (int x = fx; x < lx; x++) {
					const unsigned int i = t+(z*static_cast<int>(cf.ny)+y)*static_cast<int>(cf.nx)+x;
					assert(i < cf.nel()); v += sum(x,y,z)*cf(i); };
		g[t++] = v; };
	return g*T(2.0);
}

template <typename T>
matS<T> splineField<T>::energyHess() const
//! returns Hessian of bending energy of this spline field.
{	std::vector<enHelper> hp;							// get helpers with values for all possible overlaps
	vecD<T> vxs(3); vxs(0) = vs.x; vxs(1) = vs.y; vxs(2) = vs.z;
	for (unsigned int d1 = 0; d1 < 3; d1++) {
		for (unsigned int d2 = d1; d2 < 3; d2++) {
			uvecD dv(3); dv = 0; dv[d1]++; dv[d2]++;
			spline3d<T> ds(s3.ksp(),uvec3(dv(0),dv(1),dv(2)));
			ds /= (vxs(d1)*vxs(d2)); enHelper h(ds); 			// derivative in mm^{-1}
			h *= d1 == d2? T(2.0): T(4.0); hp.push_back(h); } };
	matS<T> M(cf.nel(),cf.nel(),s3.nnz(ex));
	unsigned int t = 0; M.ri[0] = t; spline3d<T> sp(s3.ksp());
	fwIterator fi(cf.extent(),0); uvec3 s; while (fi(s)) {
		uvec3 f,l; sp.coefRange(s,ex,f,l);
		for (unsigned int z = f.z; z < l.z; z++)
			for (unsigned int y = f.y; y < l.y; y++)
				for (unsigned int x = f.x; x < l.x; x++) { T v = T(0); 
					for (const auto& h : hp) v += h(x-s.x,y-s.y,z-s.z);
					M.ci[t] = cf.index(x,y,z); M.x[t] = v; t++; };
		M.ri[cf.index(s)+1] = t; };
	return M;
}

template <typename T>
fimage splineField<T>::save() const
//! saves coefficients of spline field into BRIAN file named fname.
{	fimage dst(cf.extent());
	for (unsigned int i = 0; i < cf.nel(); i++) dst(i) = float(cf(i));
	char buf[128]; uvec3 sp = s3.ksp();
	sprintf(buf, "%d %d %d", ex.x, ex.y, ex.z); dst.at.add({"field", buf});		// save extent of spline field
	sprintf(buf, "%d %d %d", sp.x, sp.y, sp.z); dst.at.add({"spacing", buf});	// save knot spacing
	dst.setVoxelSize(vs); return dst;						// save voxel size and write to file
}

template <typename T>
void splineField<T>::read(const fimage& src)
//! reads coefficients of spline field from BRIAN file named fname.
{	if (src.nel() == 0) return; cf = image<T>(src.extent());
	for (unsigned int i = 0; i < src.nel(); i++) cf(i) = T(src(i));
	uvec3 sp; fvec3 _vs;								// read spline field from file
	const char* s1 = src.at.lookup("field"); 					// read extent of spline field
	if (s1 == nullptr) throw rtException("Image is not a spline field");
	sscanf(s1, "%d%d%d", &ex.x, &ex.y, &ex.z);
	const char* s2 = src.at.lookup("spacing"); 					// read knot spacing
	if (s2 == nullptr) throw rtException("Image is not a spline field");
	sscanf(s2, "%d%d%d", &sp.x, &sp.y, &sp.z);
	vs = _vs; s3 = spline3d<T>(sp); updateFields();
}
#endif

