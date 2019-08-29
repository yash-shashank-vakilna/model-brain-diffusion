#ifndef IMAGE_H
#define IMAGE_H

/*
 *
 * image.h: templated 3D image
 * BRIAN Software Package Version 3.0
 *
 * $Id: image.h 488 2017-02-16 00:40:01Z kruggel $
 *
 * 0.10 (17/04/09): initial version
 * 0.20 (10/10/09): sync'ed with vector.h & matrix.h
 * 0.30 (12/05/10): FFTs added
 * 0.40 (14/05/10): filterGaussian fixed
 * 0.50 (26/11/10): transpose added
 * 0.51 (19/11/11): copyAttributes added to copy constructor
 * 0.60 (10/12/12): const added & interface changes for BRIAN 2.4
 * 0.61 (22/04/13): transpose dimensions corrected
 * 0.70 (16/12/13): documented
 * 0.80 (12/09/14): renamed field to image
 * 0.90 (07/11/15): masked Gaussian and NLM filter added
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Templated 3D image class.
*/

#define allElem(i) (unsigned int i = 0; i < nel(); i++) 

#define unaryOpScalarF(sym) \
	image<T>& operator sym (const T b) \
	{ for allElem(i) x[i] sym b; return *this; }
#define binaryOpScalarF(sym) \
	image<T> operator sym (const T b) const \
	{ image<T> a(nx,ny,nz); for allElem(i) a.x[i] = x[i] sym b; return a; }
#define unaryOpImageF(sym) \
	image<T>& operator sym (const image<T>& b) \
	{ checkdims(b); for allElem(i) x[i] sym b.x[i]; return *this; }
#define binaryOpImageF(sym) \
	image<T> operator sym (const image<T>& b) const \
	{ checkdims(b); image<T> a(nx,ny,nz); for allElem(i) a.x[i] = x[i] sym b.x[i]; return a; }

#define prepare() \
	if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny || p.z < 0 || p.z >= nz) return false; \
	float q; \
	const float dx1 = std::modf(p.x,&q), dx2 = 1.0f-dx1; const unsigned int xi = ITOU(q); \
	const float dy1 = std::modf(p.y,&q), dy2 = 1.0f-dy1; const unsigned int yi = ITOU(q); \
	const float dz1 = std::modf(p.z,&q), dz2 = 1.0f-dz1; const unsigned int zi = ITOU(q);
#define getNeighbors() \
	const unsigned int dx = p.x >= nx-1? 0: 1, dy = p.y >= ny-1? 0: 1, dz = p.z >= nz-1? 0: 1; \
	T k222 = x[index(xi, yi, zi)], k122 = x[index(xi+dx, yi, zi)]; \
	T k212 = x[index(xi, yi+dy, zi)], k112 = x[index(xi+dx, yi+dy, zi)]; \
	T k221 = x[index(xi, yi, zi+dz)], k121 = x[index(xi+dx, yi, zi+dz)]; \
	T k211 = x[index(xi, yi+dy, zi+dz)], k111 = x[index(xi+dx, yi+dy, zi+dz)];

//! Base class of all 3D images.

template <class T> class image  {
public:
	unsigned int	nx;			//!< dimension in x direction
	unsigned int	ny;			//!< dimension in y direction
	unsigned int	nz;			//!< dimension in z direction
	T		*x;			//!< image data array
	attrList	at;			//!< attribute list
private:
	void		checksize(const unsigned int _nx, const unsigned int _ny,
				const unsigned int _nz) const 				//! checks that image dimensions are less than 4GB.
			{ assert(static_cast<unsigned long>(_nx)*static_cast<unsigned long>(_ny)*static_cast<unsigned long>(_nz) < std::numeric_limits<unsigned int>::max()); }
	void		checkdims(const image<T>& b) const 				//! checks that image dimensions are compatible with image b.
			{ assert(matchesExtent(b)); }
	vec3<T>		centralGradient(const unsigned int i) const 			//! computes central gradient at index i.
			{ return T(0.5)*vec3<T>(x[i+1]-x[i-1], x[i+nx]-x[i-nx], x[i+nx*ny]-x[i-nx*ny]); }
	T		centralGradientMagnitude(const unsigned int i) const  		//! computes central gradient magnitude at index i.
			{ return norm(centralGradient(i)); }
	vec3<T>		centralGradient(const unsigned int _x, const unsigned int _y,
				const unsigned int _z) const  				//! computes central gradient at site (x,y,z).
			{ return centralGradient(index(_x,_y,_z)); }
	T		centralGradientMagnitude(const unsigned int _x,
				const unsigned int _y, const unsigned int _z) const 
			{ return centralGradientMagnitude(index(_x,_y,_z)); }
	T		centralGradientMagnitude(const uvec3& s) const   		//! computes central gradient at site s.
			{ return centralGradientMagnitude(index(s)); }
	void		convolveX(T* d, const T* s, const fvecD& wt) const;
	void		convolveY(T* d, const T* s, const fvecD& wt) const;
	void		convolveZ(T* d, const T* s, const fvecD& wt) const;
	void		maskedConvolveX(T* d, const T* s, const bool* m, const fvecD& wt) const;
	void		maskedConvolveY(T* d, const T* s, const bool* m, const fvecD& wt) const;
	void		maskedConvolveZ(T* d, const T* s, const bool* m, const fvecD& wt) const;
	void		flipX();
	void		flipY();
	void		flipZ();
	image<T>	transposeXZY(fvec3& d) const;
	image<T>	transposeYXZ(fvec3& d) const;
	image<T>	transposeZXY(fvec3& d) const;
	image<T>	transposeYZX(fvec3& d) const;
	image<T>	transposeZYX(fvec3& d) const;
protected:
	void		alloc(const unsigned int _nx, const unsigned int _ny,
				const unsigned int _nz) 				//! allocates space for image of dimensions (nx,ny,nz).
			{ nx = _nx; ny = _ny; nz = _nz; checksize(nx,ny,nz); x = new T [nel()]; }
	void		assign(const image<T>& b)					//! (re)allocates & copy from image b.
			{ if (matchesExtent(b) == false) { clear(); alloc(b.nx,b.ny,b.nz); };
		 	  for allElem(i) x[i] = b.x[i];
			  copyAttributes(b.at); }
public:
	void		resize(const unsigned int _nx, const unsigned int _ny,
				const unsigned int _nz) 				//! resize for image of dimensions (nx,ny,nz).
			{ clear(); alloc(_nx,_ny,_nz); }
	image<T>()									//! allocate an empty image.
			: nx(0), ny(0), nz(0), x(nullptr), at() { }
	image<T>(const unsigned int _nx, const unsigned int _ny, const unsigned int _nz) //! allocate an image of dimensions (nx,ny,nz).
			: nx(0), ny(0), nz(0), x(nullptr), at() { alloc(_nx,_ny,_nz); }
	image<T>(const uvec3& s) 							//! allocate an image of extent s.
			: nx(0), ny(0), nz(0), x(nullptr), at()
			{ alloc(s.x,s.y,s.z); }
	image<T>(const image<T>& b) 							//! allocate an image from image b.
			: nx(0), ny(0), nz(0), x(nullptr), at() { assign(b); }		// N.B. does not copy attributes (for speed)
	image<T>(image<T>&& b) 								//! moves from image b.
			: nx(b.nx), ny(b.ny), nz(b.nz), x(b.x), at(b.at) { b.x = nullptr; }
	image<T>(FILE* fp) 								//! allocate an image from FILE* fp.
			: nx(0), ny(0), nz(0), x(nullptr), at()
			{ read(fp); }
	~image<T>()	{ clear(); }
	image<T>&	operator=(const image<T>& b)					//! assigns from image b.
			{ if (this != &b) assign(b); return *this; }
	image<T>&	operator=(image<T>&& b)						//! move assigns from image b.
			{ assert(this != &b); at = std::move(b.at);
			  nx = b.nx; ny = b.ny; nz = b.nz;
			  delete [] x; x = b.x; b.x = nullptr; return *this; }
	unsigned int	nel() const							//! returns number of elements in this image.
			{ return nx*ny*nz; }
	unsigned int	index(const unsigned int x, const unsigned int y,
				const unsigned int z) const				//! returns index of element at site (x,y,z).
			{ return (z*ny+y)*nx+x; }
	unsigned int	index(const uvec3& s) const					//! returns index of element at site s.
			{ return index(s.x,s.y,s.z); }
	unsigned int	index(const fvec3& p) const					//! returns index of element at site s.
			{ return index(FTOU(p.x),FTOU(p.y),FTOU(p.z)); }
	T		operator()(const unsigned int _x, const unsigned int _y,
				const unsigned int _z) const				//! returns element at site (x,y,z).
			{ return x[index(_x,_y,_z)]; }
	T&		operator()(const unsigned int _x, const unsigned int _y,
				const unsigned int _z)					//! returns ref to element at site (x,y,z).
			{ return x[index(_x,_y,_z)]; }
	T		operator()(const unsigned int i) const				//! returns element at index i.
			{ return x[i]; }
	T&		operator()(const unsigned int i)				//! returns ref to element at index i.
			{ return x[i]; }
	T		operator()(const uvec3& s) const				//! returns element at site s.
			{ return x[index(s)]; }
	T&		operator()(const uvec3& s)					//! returns ref to element at site s.
			{ return x[index(s)]; }
	bool		matchesExtent(const image<T>& b) const				//! checks if this image has same dimensions as b.
			{ return nx == b.nx && ny == b.ny && nz == b.nz; }
	bool		interpolateAt(const fvec3& p, T& v) const			//! interpolate image at point p, returns result in v.
			{ v = 0; prepare(); getNeighbors(); 
			  k111 = k111*dx1+k211*dx2; k121 = k121*dx1+k221*dx2; k112 = k112*dx1+k212*dx2; k122 = k122*dx1+k222*dx2;
			  k111 = k111*dy1+k121*dy2; k112 = k112*dy1+k122*dy2; v = k111*dz1 + k112*dz2; return true; }
	bool		interpolateGradientAt(const fvec3& p, T& v, vec3<T>& g) const
			//! interpolate image and gradient at point p, returns intensity in v and gradient in g.
			{ v = 0; g = 0; prepare(); getNeighbors(); 
			  g.x = (((k111-k211)*dy1+(k121-k221)*dy2))*dz1+(((k112-k212)*dy1+(k122-k222)*dy2))*dz2;
			  k111 = k111*dx1+k211*dx2; k121 = k121*dx1+k221*dx2; k112 = k112*dx1+k212*dx2; k122 = k122*dx1+k222*dx2;
			  g.y = (k111-k121)*dz1 + (k112-k122)*dz2; k111 = k111*dy1+k121*dy2; k112 = k112*dy1+k122*dy2;
			  g.z = k111-k112; v = k111*dz1+k112*dz2; return true; }
	T		gradX(const uvec3& s) const					//! computes central gradient in x direction at site s.
			{ const unsigned int i = index(s); return x[i+1]-x[i-1]; }
	T		gradY(const uvec3& s) const					//! computes central gradient in y direction at site s.
			{ const unsigned int i = index(s); return x[i+nx]-x[i-nx]; }
	T		gradZ(const uvec3& s) const 					//! computes central gradient in z direction at site s.
			{ const unsigned int i = index(s); return x[i+nx*ny]-x[i-nx*ny]; }
	T		laplaceX(const uvec3& s) const					//! computes Laplacian in x direction at site s.
			{ const unsigned int i = index(s); return x[i+1]-x[i]*2.0f+x[i-1]; }
	T		laplaceY(const uvec3& s) const					//! computes Laplacian in y direction at site s.
			{ const unsigned int i = index(s); return x[i+nx]-x[i]*2.0f+x[i-nx]; }
	T		laplaceZ(const uvec3& s) const 					//! computes Laplacian in z direction at site s.
			{ const unsigned int i = index(s); return x[i+nx*ny]-x[i]*2.0f+x[i-nx*ny]; }
	T		gradXY(const uvec3& s) const					//! computes central gradient in xy direction at site s.
			{ const unsigned int i = index(s); return x[i-nx-1]-x[i-nx+1]+x[i+nx+1]-x[i+nx-1]; }
	T		gradXZ(const uvec3& s) const 					//! computes central gradient in xz direction at site s.
			{ const unsigned int i = index(s); return x[i-nx*ny-1]-x[i-nx*ny+1]+x[i+nx*ny+1]-x[i+nx*ny-1]; }
	T		gradYZ(const uvec3& s) const 					//! computes central gradient in yz direction at site s.
			{ const unsigned int i = index(s); return x[i-nx*ny-nx]-x[i-nx*ny+nx]+x[i+nx*ny+nx]-x[i+nx*ny-nx]; }
	vec3<T>		centralGradient(const uvec3& s) const				//! computes central gradient at site s.
			{ return centralGradient(index(s)); }
	image<T>	scale(uvec3 ex)							//! returns a copy scaled by extent ex.
			{ float fx = float(nx)/ex.x, fy = float(ny)/ex.y, fz = float(nz)/ex.z;
			  image<T> a(ex); unsigned int i = 0;
			  for (unsigned int zi = 0; zi < ex.z; zi++)
			  	  for (unsigned int yi = 0; yi < ex.y; yi++)
			  		  for (unsigned int xi = 0; xi < ex.x; xi++, i++) {
							fvec3 p(xi*fx, yi*fy, zi*fz); interpolateAt(p, a.x[i]); };
			  *this = a; return *this; }
			unaryOpScalarF(=)
			unaryOpScalarF(+=)
			unaryOpScalarF(-=)
			unaryOpScalarF(*=)
			unaryOpScalarF(/=)
			binaryOpScalarF(+)
			binaryOpScalarF(-)
			binaryOpScalarF(*)
			binaryOpScalarF(/)
			unaryOpImageF(+=)
			unaryOpImageF(-=)
			binaryOpImageF(+)
			binaryOpImageF(-)
			binaryOpImageF(*)						// voxelwise product
	uvec3		extent() const							//! returns image extent.
			{ return uvec3(nx,ny,nz); }
	uvec3		site(const unsigned int i) const				//! returns site at index i.
			{ unsigned int z = i/(nx*ny), r = i%(nx*ny); return uvec3(r%nx,r/nx,z); }
	void		copyAttributes(const attrList& list)				//! copy attributes from list.
			{ at = list; }
	bool		getVoxelSize(fvec3& v) const					//! returns voxel size in real world dimensions.
			{ const char *s = at.lookup("voxel");
			  if (s) sscanf(s, "%f%f%f", &v.x, &v.y, &v.z);
			  return s != nullptr; }
	void		setVoxelSize(const fvec3& v)					//! sets voxel size in real world dimensions.
			{ char buf[40]; sprintf(buf, "%.2f %.2f %.2f", v.x, v.y, v.z);
			  at.update("voxel", buf); }
	bool		getRepnType(repnType& t) const					//! returns image value representation in t.
			{ const char *s = at.lookup("repn");
			  if (s) t = lookupRepn(s);
			  return s != nullptr; }
	float		sumValues() const						//! returns scalar sum of all elements.
			{ T s = 0; for allElem(i) { s += x[i]; }; return s.scalar(); }
	float		sum() const 							//! returns float sum of all elements.
			{ float s = 0; for allElem(i) s += x[i]; return s; }
	void		clear() { delete [] x; x = nullptr; nx = 0; ny = 0; nz = 0; at.clear(); }
	image<vec3<T> > centralGradient() const;
	image<T>	centralGradientMagnitude() const;
	image<T> 	filterGaussian(const float sigma, const unsigned int deg = 0) const;
	image<T> 	maskedFilterGaussian(const image<bool>& mask, const float sigma, const unsigned int deg = 0) const;
	image<T>	cut(const unsigned int b, uvec3& org) const;
	image<T>	paste(const uvec3& org, const uvec3& ex) const;
	image<unsigned int> label(const connectivity cn, unsigned int& lmax) const;
	image<float>	chamferDTF() const;
	image<float>	euclidianDTF() const;
	image<bool>	selectBig(const unsigned int nc) const;
	image<bool>	selectAt(const connectivity cn, const uvec3& org);
	image<bool> 	inv() const;
	image<bool>	binarize(const T vmin, const T vmax = std::numeric_limits<T>::max()) const;
	image<bool> 	rankFilter(const connectivity cn) const;
	image<bool>	dilate(const float dist) const;
	image<bool>	erode(const float dist) const;
	image<bool>	open(const float dist) const;
	image<bool>	close(const float dist) const;
	void		range(T& vmin, T& vmax) const;
	fvecD		histogram(const unsigned int nbins) const;
	fvecD		histogram(const unsigned int nbins, const float vmin, const float vmax) const;
	float		estimateNoise(const unsigned int nsamp, const unsigned int nbins) const;
	image<float>	anisotropicDiffusion(const float lambda, const unsigned int nit) const;
	image<fvec3>	gradient() const;						// Zucker-Hummel gradient filter
	image<fvec3>	filteredGradient(const float sigma) const;			// Gaussian-filtered Zucker-Hummel gradients
	image<float>	gradientMagnitude() const;
	image<float>	filteredGradientMagnitude(const float sigma) const;
	float		ssd(const image<float>& b) const;
	float		cc(const image<float>& ref, const fmat3& m) const;
	float		histogramScale(const unsigned int nbins, const float fr) const;
	float		nmi(const image<float>& ref, const fmat3& m) const;
	image<T>	deform(const image<fvec3>& f) const;
	image<T>	deformLabel(const image<fvec3>& f) const;
	image<T>	inverseDeform(const image<fvec3>& f) const;
	image<T>	inverseLabelDeform(const image<fvec3>& f) const;
	bool		interpolateAt(const fvec3& p) const;
	int		genus(const connectivity cn) const;
	bool		isSimple(const uvec3& s, const connectivity cn) const;
	bool		isSimple(const unsigned int id, const connectivity cn) const;
	bool		isSimpleInBG(const uvec3& s, const connectivity cn) const;
	unsigned int	removeLayer(const connectivity cn);
	unsigned int	growLayer(const image<bool>& dom, const image<float>& dist, const connectivity cn);
	image<bool>	growZero(const connectivity cn) const;
	image<float>	bsplineGradientMagnitude(const uvec3& ex, float& gmax) const;
	image<float>	transform(const uvec3& ex, const fmat3& m) const;
	image<bool>	thinning() const;
	void		read(is& in, const repnType r);
	void		fromBundle(const bundle& b);
	void		read(FILE* fp)							//! reads the first image in file* fp
			{ bundleList bl; readFile(fp,bl); fromBundle(*bl.begin()); }
	void		read(const char *fname)
			{ FILE* fp = openFile(fname,"r"); read(fp); closeFile(fp); }
	bundle		toBundle(const char* tag, repnType r = repnType::unknown) const;
	void		save(FILE* fp, const repnType r = repnType::unknown) const	//! saves an image in file* fp
			{ bundleList bl; bl.push_back(toBundle("image",r)); writeFile(fp,bl); }
	void		save(const char* fname, const repnType r = repnType::unknown) const //! saves an image in file* fname
			{ FILE* fp = openFile(fname,"w"); save(fp,r); closeFile(fp); }
	image<bool>	zeroCrossings() const;
	image<float>	constrainedDTF(const unsigned int in, const unsigned int out) const;
	image<bool>	boundary(const connectivity cn) const;
	image<T>	sumOver(const unsigned int w) const;
	image<unsigned int> classifyImage(const unsigned int nbins, const unsigned int nc,
				const image<bool>& mask = image<bool>()) const;
	void		computeCSFparams(float& mu, float& std) const;
	void		filterGaussianU(const float sigma, const unsigned int deg = 0);
	void		maskedFilterGaussianU(const image<bool>& mask, const float sigma, const unsigned int deg = 0);
	image<bool>	selectBig(const connectivity cn);
	image<float>	mask(const image<bool>& b, const bool bg) const;
	void		adaptIntensity(const image<T>& b);
	image<bool>	fillHoles();
	void		moments(FILE *outf, const connectivity cn, const bool detail) const;
	void		transpose(char *code);
	bool		interpolateLabelAt(const fvec3& p, unsigned int& lmax) const;
	bool		interpolateLabelAt(const fvec3& p, bool& b) const;
	image<T>	transformLabels(const uvec3& ex, const fmat3& m) const;
	void		displayRange(float& qmin, float& qmax) const;
	image<float>	medianFilter(const connectivity cn) const; 
	image<float>	standardize() const; 
	image<T>	getSubimage(const uvec3& s) const;
	bool		finite() const							//! checks that all elements are finite.
			{ for allElem(i) if (std::isfinite(x[i]) == false) return false; return true; }
	image<float>	affineTransform(const uvec3& ex, const fmat3& m) const;
	float		estimateNoise() const;
	image<float>	filterMean(const unsigned int d) const;
	image<float>	filterVar(const image<float>& mn, const unsigned int d) const;
	image<T>	filterLee(const unsigned int w, const float d) const;
	float		intensityDist(const uvec3& s, const uvec3& t, const unsigned int d) const;
	image<T>	filterNLM(const unsigned int m = 5, const unsigned int d = 2, const bool verbose = false) const;
	image<T>	filterNLMRician(const unsigned int m = 5, const unsigned int d = 2, const bool verbose = false) const;
	float		ssd(const image<T>& ref, const image<bool>& mask) const;
	image<T>	insertMatrix(const matD<T>& M, const unsigned int z);
};

template<typename T> image<std::complex<T>> fft(const image<T>&);			// compute a forward FFT of an image, see fft.C
template<typename T> image<T> ifft(const image<std::complex<T>>&);			// compute an inverse FFT of a complex image, see fft.C
template<typename T> image<vec3<std::complex<T>>> fft(const image<vec3<T>>&);		// compute a forward FFT of a vector image, see fft.C
template<typename T> image<vec3<T>> ifft(const image<vec3<std::complex<T>>>&);		// compute an inverse FFT of a complex vector image, see fft.C

#undef unaryOpScalarF
#undef binaryOpScalarF
#undef unaryOpImageF
#undef binaryOpImageF
#undef prepare
#undef getNeighbors
#undef allElem

using dimage = image<double>;
using fimage = image<float>;
using limage = image<unsigned int>;
using bimage = image<bool>;
using cimage = image<fcomplex>;
using fvec3image = image<fvec3>;
using cvec3image = image<cvec3>;
using vbimage = std::vector<bimage>;
using vlimage = std::vector<limage>;
using vfimage = std::vector<fimage>;
#endif
