#ifndef IMAGEFUNCS_H
#define IMAGEFUNCS_H

/*
 *
 * imageFuncs.h: templated 3D image functions
 * BRIAN Software Package Version 3.0
 *
 * $Id: imageFuncs.h 504 2017-03-16 23:28:00Z kruggel $
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

template <typename T>
void image<T>::fromBundle(const bundle& b)
//! reads an image from bundle b.
{	const char *s = b.at.lookup("repn");
	if (s == 0) throw optException("Unknown representation for fimage");
	repnType r = lookupRepn(s);
	const unsigned int _nx = std::max(b.at.lookupNumber("ncolumns"), 1u);
	const unsigned int _ny = std::max(b.at.lookupNumber("nrows"), 1u);
	const unsigned int _nz = std::max(b.at.lookupNumber("nbands"), 1u);
	resize(_nx,_ny,_nz); is in(b.data,b.length);
	in.readArray(x,nel(),r); copyAttributes(b.at);
}

template <typename T>
bundle image<T>::toBundle(const char* tag, repnType r) const
//! saves this image in bundle b.
{	os out; out.saveArray(x,nel(),r); char buf[LINELEN];
	bundle b; b.copyStream(out.data()); b.copyAttributes(at);
	b.repn = repnType::image; b.updateValue(tag);					// add tag for image type
	sprintf(buf, "%ld", 0L); b.at.update("data", buf);				// updated when bundle is written
	sprintf(buf, "%zd", b.length); b.at.update("length", buf);			// updated when bundle is written
	sprintf(buf, "%u", nx); b.at.update("ncolumns", buf);				// update automatic attributes
	sprintf(buf, "%u", ny); b.at.update("nrows", buf);
	sprintf(buf, "%u", nz); b.at.update("nbands", buf);
	b.at.update("nframes", buf);
	if (r == repnType::unknown) {
		if (typeid(T) == typeid(unsigned char)) r = repnType::byte;
		else if (typeid(T) == typeid(short int)) r = repnType::shrt;
		else if (typeid(T) == typeid(unsigned int)) r = repnType::sint;
		else if (typeid(T) == typeid(float)) r = repnType::flt; };
	b.at.update("repn", repnTable[static_cast<int>(r)].name);
	return b;
}

template <typename T>
void image<T>::convolveX(T *d, const T *s, const fvecD &wt) const 
//! 1d convolution in x direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int zi = 0; zi < nz; zi++) {
		for (unsigned int yi = 0; yi < ny; yi++) {
			for (unsigned int xi = 0; xi < nx; xi++) {
				const unsigned int x0 = xi < wn2? 0: xi-wn2;
				const unsigned int x1 = xi+wn2 > nx-1? nx-1: xi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(xi-x0)];
				for (unsigned int ti = x0; ti <= x1; ti++, wi++) {
					v += s[index(ti,yi,zi)]*(*wi); w += *wi; }
				const unsigned int i = index(xi,yi,zi);
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

template <typename T>
void image<T>::convolveY(T *d, const T *s, const fvecD &wt) const 
//! 1d convolution in y direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int zi = 0; zi < nz; zi++) {
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int yi = 0; yi < ny; yi++) {
				const unsigned int y0 = yi < wn2? 0: yi-wn2;
				const unsigned int y1 = yi+wn2 > ny-1? ny-1: yi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(yi-y0)];
				for (unsigned int ti = y0; ti <= y1; ti++, wi++) {
					v += s[index(xi,ti,zi)]*(*wi); w += *wi; }
				const unsigned int i = index(xi,yi,zi);
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

template <typename T>
void image<T>::convolveZ(T *d, const T *s, const fvecD &wt) const 
//! 1d convolution in z direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int yi = 0; yi < ny; yi++) {
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int zi = 0; zi < nz; zi++) {
				const unsigned int z0 = zi < wn2? 0: zi-wn2;
				const unsigned int z1 = zi+wn2 > nz-1? nz-1: zi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(zi-z0)];
				for (unsigned int ti = z0; ti <= z1; ti++, wi++) {
					v += s[index(xi,yi,ti)]*(*wi); w += *wi; }
				const unsigned int i = index(xi,yi,zi);
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

template <typename T>
void image<T>::maskedConvolveX(T *d, const T *s, const bool *m, const fvecD &wt) const 
//! masked 1d convolution in x direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int zi = 0; zi < nz; zi++) {
		for (unsigned int yi = 0; yi < ny; yi++) {
			for (unsigned int xi = 0; xi < nx; xi++) {
				const unsigned int i = index(xi,yi,zi); if (m[i] == false) continue;
				const unsigned int x0 = xi < wn2? 0: xi-wn2;
				const unsigned int x1 = xi+wn2 > nx-1? nx-1: xi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(xi-x0)];
				for (unsigned int ti = x0; ti <= x1; ti++, wi++) {
					const unsigned int j = index(ti,yi,zi);
					if (m[j] == false) continue;
					v += s[j]*(*wi); w += *wi; };
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

template <typename T>
void image<T>::maskedConvolveY(T *d, const T *s, const bool *m, const fvecD &wt) const 
//! masked 1d convolution in y direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int zi = 0; zi < nz; zi++) {
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int yi = 0; yi < ny; yi++) {
				const unsigned int i = index(xi,yi,zi); if (m[i] == false) continue;
				const unsigned int y0 = yi < wn2? 0: yi-wn2;
				const unsigned int y1 = yi+wn2 > ny-1? ny-1: yi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(yi-y0)];
				for (unsigned int ti = y0; ti <= y1; ti++, wi++) {
					const unsigned int j = index(xi,ti,zi);
					if (m[j] == false) continue;
					v += s[j]*(*wi); w += *wi; };
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

template <typename T>
void image<T>::maskedConvolveZ(T *d, const T *s, const bool *m, const fvecD &wt) const 
//! masked 1d convolution in z direction.
{	const unsigned int wn2 = wt.N/2;
	for (unsigned int yi = 0; yi < ny; yi++) {
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int zi = 0; zi < nz; zi++) {
				const unsigned int i = index(xi,yi,zi); if (m[i] == false) continue;
				const unsigned int z0 = zi < wn2? 0: zi-wn2;
				const unsigned int z1 = zi+wn2 > nz-1? nz-1: zi+wn2;
				T v{0}; float w = 0.0f; const float *wi = &wt.x[wn2-(zi-z0)];
				for (unsigned int ti = z0; ti <= z1; ti++, wi++) {
					const unsigned int j = index(xi,yi,ti);
					if (m[j] == false) continue;
					v += s[j]*(*wi); w += *wi; };
				d[i] = w > 0.0f? v/w: s[i]; } } };
}

inline fvecD gaussianParams(const float sigma, const unsigned int deg)
//! compute Gaussian derivative filter kernel, derivative degree [0,2].
{	unsigned int n = std::max(FTOU(6.0f*std::abs(sigma)),3u); if (n % 2 == 0) n++;	// ensure unary filter width
	fvecD m(n); float s = std::sqrt(float(2.0*M_PI));				// allocate coefficient array
	switch (deg) {									// normalization factor for degree deg
	case 0: s *= sigma; break;
	case 1: s *= std::pow(sigma,3.0f); break;
	case 2: s *= std::pow(sigma,5.0f); break;
	default: throw optException("gaussianParams: Unimplemented derivative degree %u", int(deg)); };
	for (unsigned int i = 0; i < n; i++) {
		float v = float(int(i)-int(n/2)), t = std::exp(-0.5f*SQR(v/sigma)); 	// compute coefficients
		switch (deg) {								// modify for degree deg
		case 0: break;
		case 1: t *= -v; break;
		case 2: t *= (v-sigma)*(v+sigma); break; };
		m[i] = t/s; };								// and normalize
	return m/m.sum();								// return coefficient vector
}

template <typename T>
image<T> image<T>::filterGaussian(const float sigma, const unsigned int deg) const 
//! returns image filtered by a Gaussian derivative filter kernel of size sigma and degree [0,2].
{	image<T> a(*this); a.filterGaussianU(sigma, deg); return a;
}

template <typename T>
void image<T>::filterGaussianU(const float sigma, const unsigned int deg)
//! filters this image by a Gaussian derivative filter kernel of size sigma and degree [0,2].
{	if (sigma == 0) return;
	const fvecD wt = gaussianParams(sigma, deg); image<T> t(*this);
	convolveX(x, t.x, wt); convolveY(t.x, x, wt); convolveZ(x, t.x, wt);
}

template <typename T>
image<T> image<T>::maskedFilterGaussian(const image<bool>& mask, const float sigma, const unsigned int deg) const 
//! returns image filtered by a Gaussian derivative filter kernel of size sigma and degree [0,2], masked by boolean image mask.
{	image<T> a(*this); a.maskedFilterGaussianU(mask, sigma, deg); return a;
}

template <typename T>
void image<T>::maskedFilterGaussianU(const image<bool>& mask, const float sigma, const unsigned int deg)
//! filters this image by a Gaussian derivative filter kernel of size sigma and degree [0,2], masked by boolean image mask.
{	if (sigma == 0) return;
	const fvecD wt = gaussianParams(sigma, deg); image<T> t(*this);
	maskedConvolveX(x, t.x, mask.x, wt);
	maskedConvolveY(t.x, x, mask.x, wt);
	maskedConvolveZ(x, t.x, mask.x, wt);
}

//! \relates image
template <typename T>
image<vec3<T> > image<T>::centralGradient() const
//! returns central gradient image.
{	uvec3 ex = extent(); image<vec3<T> > dst(ex); dst = T(0.0);
	allSites(s,1) dst(s) = centralGradient(s);
	return dst;
}

//! \relates image
template <typename T>
image<T> image<T>::centralGradientMagnitude() const
//! returns central gradient magnitude image.
{	uvec3 ex = extent(); image<T> dst(ex); dst = T(0.0);
	allSites(s,1) dst(s) = centralGradientMagnitude(s);
	return dst;
}

template <typename T>
image<T> image<T>::cut(const unsigned int b, uvec3& org) const
//! determines bounding box of binary objects in this image, adds margin b; returns subimage and origin of box.
{	unsigned int xmin = nx, xmax = 0, ymin = ny, ymax = 0, zmin = nz, zmax = 0;
	for (unsigned int zi = 0; zi < nz; zi++) {
		for (unsigned int yi = 0; yi < ny; yi++) {
			for (unsigned int xi = 0; xi < nx; xi++) {
				if ((*this)(xi,yi,zi) == 0) continue;
				if (xi < xmin) xmin = xi;
				if (xi > xmax) xmax = xi;
				if (yi < ymin) ymin = yi;
				if (yi > ymax) ymax = yi;
				if (zi < zmin) zmin = zi;
				if (zi > zmax) zmax = zi; } } };
	xmin = (xmin < b)? 0: xmin-b; xmax = (xmax+b > nx)? nx: xmax+b;
	ymin = (ymin < b)? 0: ymin-b; ymax = (ymax+b > ny)? ny: ymax+b;
	zmin = (zmin < b)? 0: zmin-b; zmax = (zmax+b > nz)? nz: zmax+b;
	unsigned int ex = xmax-xmin, ey = ymax-ymin, ez = zmax-zmin;
	image<T> dst(ex, ey, ez);
	for (unsigned int zi = zmin; zi < zmax; zi++)
		for (unsigned int yi = ymin; yi < ymax; yi++)
			for (unsigned int xi = xmin; xi < xmax; xi++)
				dst(xi-xmin,yi-ymin,zi-zmin) = (*this)(xi,yi,zi);
	org = uvec3(xmin,ymin,zmin); return dst;
}

template <typename T>
image<T> image<T>::paste(const uvec3& org, const uvec3& ex) const
//! pastes this image at site org into a new image of extent ex.
{	assert(org.x < ex.x && org.y < ex.y && org.z < ex.z);
	unsigned int xmax = std::min(nx,ex.x-org.x), ymax = std::min(ny,ex.y-org.y), zmax = std::min(nz,ex.z-org.z);
	image<T> dst(ex); dst = T(0.0);
	for (unsigned int zi = 0; zi < zmax; zi++)
		for (unsigned int yi = 0; yi < ymax; yi++)
			for (unsigned int xi = 0; xi < xmax; xi++)
				dst(xi+org.x,yi+org.y,zi+org.z) = (*this)(xi,yi,zi);
	return dst;
}

template <typename T>
void image<T>::flipX()
//! flip image along x in-place.
{	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx/2; xi++) { 
				T u = (*this)(xi,yi,zi);
				(*this)(xi,yi,zi) = (*this)(nx-xi-1,yi,zi);
				(*this)(nx-xi-1,yi,zi) = u; };
}

template <typename T>
void image<T>::flipY()
//! flip image along y in-place.
{	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny/2; yi++)
			for (unsigned int xi = 0; xi < nx; xi++) {
				T u = (*this)(xi,yi,zi);
				(*this)(xi,yi,zi) = (*this)(xi,ny-yi-1,zi);
				(*this)(xi,ny-yi-1,zi) = u; };
}

template <typename T>
void image<T>::flipZ()
//! flip image along z in-place.
{	for (unsigned int zi = 0; zi < nz/2; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++) {
				T u = (*this)(xi,yi,zi);
				(*this)(xi,yi,zi) = (*this)(xi,yi,nz-zi-1);
				(*this)(xi,yi,nz-zi-1) = u; };
}

template <typename T>
image<T> image<T>::transposeXZY(fvec3& d) const
//! transpose xzy out-of-place.
{	float t = d.z; d.z = d.y; d.y = t; image<T> dst(nx, nz, ny); unsigned int i = 0;
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++, i++) dst(xi,zi,yi) = x[i];
	return dst;
}

template <typename T>
image<T> image<T>::transposeYXZ(fvec3& d) const
//! transpose image yxz out-of-place.
{	float t = d.x; d.x = d.y; d.y = t; image<T> dst(ny, nx, nz); unsigned int i = 0;
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++, i++) dst(yi,xi,zi) = x[i];
	return dst;
}

template <typename T>
image<T> image<T>::transposeZXY(fvec3& d) const
//! transpose image zxy out-of-place.
{	float t = d.z; d.z = d.y; d.y = d.x; d.x = t; image<T> dst(nz, nx, ny); unsigned int i = 0;
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++, i++) dst(zi, xi, yi) = x[i];
	return dst;
}

template <typename T>
image<T> image<T>::transposeYZX(fvec3& d) const
//! transpose image yzx out-of-place.
{	float t = d.z; d.z = d.x; d.x = d.y; d.y = t; image<T> dst(ny, nz, nx); unsigned int i = 0;
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++, i++) dst(yi,zi,xi) = x[i];
	return dst;
}

template <typename T>
image<T> image<T>::transposeZYX(fvec3& d) const
//! transpose image zyx out-of-place.
{	float t = d.z; d.z = d.x; d.x = t; image<T> dst(nz, ny, nx); unsigned int i = 0;
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++)
			for (unsigned int xi = 0; xi < nx; xi++, i++) dst(zi,yi,xi) = x[i];
	return dst;
}

inline bool checkCode(const unsigned char c)
{	return c == 'R' || c == 'P' || c == 'I'; }

inline void changeCode(unsigned char&& c)
{	if (c == 'R') c = 'L'; else if (c == 'P') c = 'A'; else if (c == 'I') c = 'S'; }

template <typename T>
void image<T>::transpose(char* code)
//! transpose image according to code into orientation LAS.
{	fvec3 d = 1; const char *s = at.lookup("voxel"); if (s) sscanf(s, "%f%f%f", &d.x, &d.y, &d.z);
	image<T> dst; bool cp = false;
	if (checkCode(code[0])) { flipX(); changeCode(code[0]); };
	if (checkCode(code[1])) { flipY(); changeCode(code[1]); };
	if (checkCode(code[2])) { flipZ(); changeCode(code[2]); };
	if (code[0] == 'L') { if (code[1] != 'A') { dst = transposeXZY(d); cp = true; } }
	else if (code[0] == 'A') { dst = (code[1] == 'L')? transposeYXZ(d): transposeZXY(d); cp = true; }
	else if (code[0] == 'S') { dst = (code[1] == 'L')? transposeYZX(d): transposeZYX(d); cp = true; };
	if (cp) { dst.copyAttributes(at); *this = dst; };
	setVoxelSize(d); at.update("orientation", "LAS"); strcpy(code, "LAS");
}

template <typename T>
T sumWindow(const vecD<T>& line, const unsigned int i0, const unsigned int i1)
{	T s = 0; for (unsigned int i = i0; i <= i1; i++) s += line(i); return s; }

template <typename T>
image<T> image<T>::sumOver(const unsigned int w) const
//! returns sum over 3D window with half-width w, fast separable version.
{	image<T> dst(nx, ny, nz); dst = *this; vecD<T> lx(nx), ly(ny), lz(nz); 
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int yi = 0; yi < ny; yi++) {
			for (unsigned int xi = 0; xi < nx; xi++) lx[xi] = dst(xi,yi,zi);
			for (unsigned int xi = 0; xi < nx; xi++) {
				unsigned int i0 = xi < w? 0: xi-w;
				unsigned int i1 = xi+w > nx-1? nx-1: xi+w;
				dst(xi,yi,zi) = sumWindow(lx,i0,i1); } };
	for (unsigned int zi = 0; zi < nz; zi++)
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int yi = 0; yi < ny; yi++) ly[yi] = dst(xi,yi,zi);
			for (unsigned int yi = 0; yi < ny; yi++) {
				unsigned int i0 = yi < w? 0: yi-w;
				unsigned int i1 = yi+w > ny-1? ny-1: yi+w;
				dst(xi,yi,zi) = sumWindow(ly,i0,i1); } };
	for (unsigned int yi = 0; yi < ny; yi++) 
		for (unsigned int xi = 0; xi < nx; xi++) {
			for (unsigned int zi = 0; zi < nz; zi++) lz[zi] = dst(xi,yi,zi);
			for (unsigned int zi = 0; zi < nz; zi++) {
				unsigned int i0 = zi < w? 0: zi-w;
				unsigned int i1 = zi+w > nz-1? nz-1: zi+w;
				dst(xi,yi,zi) = sumWindow(lz,i0,i1); } };
	return dst;
}

template <typename T>
image<T> image<T>::getSubimage(const uvec3& s) const
//! returns 3x3x3 subimage at site s.
{	image<T> dst(3,3,3); dst = 0; unsigned int i = 0;
	if (s.x < 1 || s.x >= nx-1 || s.y < 1 || s.y >= ny-1 || s.z < 1 || s.z >= nz-1) return dst;
	for (unsigned int z = s.z-1; z <= s.z+1; z++)
		for (unsigned int y = s.y-1; y <= s.y+1; y++)
			for (unsigned int x = s.x-1; x <= s.x+1; x++, i++)
				dst(i) = (*this)(x,y,z);
	return dst;
}

template <typename T>
float image<T>::ssd(const image<T>& ref, const image<bool>& mask) const
//! computes similarity of src and ref.
{	float s = 0.0f; unsigned int n = 0;
	for (unsigned int i = 0; i < nel(); i++) { if (mask(i) == false) continue;
		const float d = static_cast<float>(x[i]-ref(i)); s += SQR(d); n++; };
	return n? s/static_cast<float>(n): 0.0f;
}

template <typename T>
float image<T>::estimateNoise() const
//! returns variance of Gaussian noise in image src.
// Gasse, T., Sroka L., Steinmetz C. (1986)
// Residual variance and residual pattern in non linear regression.
// Biometrika 73, 625-633.
{	const uvec3 ex = extent(); float s2 = 0.0f; unsigned int n = 0;
	allSites(s,1) { nbIterator ni(ex,s,connectivity::c6); uvec3 t;
		float v = 6.0f*static_cast<float>((*this)(s));
		while (ni(t)) v -= static_cast<float>((*this)(t));
		s2 += SQR(v); n++; };
	return n? s2/static_cast<float>(n*42): 0.0f;
}

template <typename T>
image<float> image<T>::filterMean(const unsigned int d) const
//! computes local means using box width d; returns float image.
{	const uvec3 ex = extent(); image<float> dst(ex); dst = 0.0f;
	allSites(s,0) { boxIterator bx(ex,0,s,d); uvec3 t; 
		float m = 0.0f; unsigned int n = 0;
		while (bx(t)) { m += static_cast<float>((*this)(t)); n++; };
		dst(s) = n? m/static_cast<float>(n): static_cast<float>((*this)(s)); };
	return dst;
}

template <typename T>
image<float> image<T>::filterVar(const image<float>& mn, const unsigned int d) const
//! computes local variances using box width d.
{	const uvec3 ex = extent(); image<float> dst(ex); dst = 0.0f;
	allSites(s,0) { boxIterator bx(ex,0,s,d); uvec3 t; 
		float m = mn(s), s2 = 0.0f; unsigned int n = 0;
		while (bx(t)) { s2 += SQR(static_cast<float>((*this)(t))-m); n++; };
		dst(s) = n > 1? s2/static_cast<float>(n-1): 0; };
	return dst;
}

template <typename T>
image<T> image<T>::filterLee(const unsigned int w, const float d) const
//! performs adaptive filtering in box of extent w with max intensity difference d.
// Lee J.S. (1983)
// Digital image smoothing and the sigma filter.
// Comput. Vis., Graph. Image Process. 24, 255-269.
{	const uvec3 ex = extent(); image<T> dst(ex); dst = T(0);
	allSites(s,0) { boxIterator bx(ex,0,s,w/2); uvec3 t;
		const float v = static_cast<float>((*this)(s));
		float vs = 0; unsigned int n = 0;
		while (bx(t)) { if (std::abs(static_cast<float>((*this)(t))-v) < d) {
			vs += static_cast<float>((*this)(t)); n++; } };
		dst(s) = n? T(vs/n): T(v); };
	return dst;
}

template <typename T>
float image<T>::intensityDist(const uvec3& s, const uvec3& t, const unsigned int d) const
//! returns squared intensity difference between boxes of half-width w around site s and t. 
{	const uvec3 ex = extent(); boxIterator bx(ex,0,s,d); uvec3 u; 
	float s2 = 0; unsigned int n = 0; while (bx(u)) {
		int x = int(t.x+u.x)-int(s.x); if (x < 0 || x >= int(ex.x)) continue;
		int y = int(t.y+u.y)-int(s.y); if (y < 0 || y >= int(ex.y)) continue;
		int z = int(t.z+u.z)-int(s.z); if (z < 0 || z >= int(ex.z)) continue;
		float v = static_cast<float>((*this)(u)-(*this)(x,y,z)); s2 += SQR(v); n++; };
	return n? s2/n: 0;
}

template <typename T>
image<T> image<T>::filterNLM(const unsigned int m, const unsigned int d, const bool verbose) const
//! performs a non-local means filter using neighborhood m and filter box width d.
// Coupe P., Yger P., Prima S., Hellier P., Kervrann C., Barillot C. (2008)
// An Optimized Blockwise Nonlocal Means Denoising Filter for 3D Magnetic Resonance Images
// IEEE TMI 27, 425-441.
{	const float mu1 = 0.95f, si1 = 0.5f, beta = 1.0f;
	const uvec3 ex = extent(); const float s2 = estimateNoise();
	const image<float> mn = filterMean(d), var = filterVar(mn, d);			// get local mean and variance image
	image<T> dst(ex); dst = 0.0f; allSites(s,0) { float vs = 0.0f, ws = 0.0f;		// for all sites s
		if (verbose && s.x == 0 && s.y == 0) { printf("%u\r", s.z); fflush(stdout); };
		if (mn(s) != 0.0f && var(s) != 0.0f) {					// do not work on background or constant regions
			boxIterator bx(ex,0,s,m); uvec3 t; while (bx(t)) {		// for all sites t in a box around s
				const float mr = mn(s)/mn(t);				// select voxels (Eq 7)...
				if (mu1 > mr || mr > 1.0f/mu1) continue;		// ...by means
				const float vr = var(s)/var(t);
				if (si1 > vr || vr > 1.0f/si1) continue;		// ...and by variances
				const float v = intensityDist(s, t, d);			// compute local intensity difference (Eq 5)
				const float w = std::exp(-v/(2.0f*beta*s2));		// compute weight (Eq 6)
				vs += static_cast<float>((*this)(t))*w; ws += w; } };
		dst(s) = ws > 0.0f? T(vs/ws): (*this)(s); };
	return dst;
}

template <typename T>
image<T> image<T>::filterNLMRician(const unsigned int m, const unsigned int d, const bool verbose) const
//! performs a non-local means filter using neighborhood m and filter box width d.
// Wiest-Daessle N., Prima S., Coupe P., Morrissey S.P., Barillot C. (2008)
// Rician Noise Removal by Non-Local Means Filtering for Low Signal-to-Noise Ratio MRI: Applications to DT-MRI.
// Proc MICCAI LCNS 5242, 171-179.
{	const float mu1 = 0.95f, si1 = 0.5f, beta = 1.0f;
	const uvec3 ex = extent(); const float s2 = estimateNoise();
	const image<float> mn = filterMean(d), var = filterVar(mn, d);			// get local mean and variance image
	const unsigned int sz = (2*m+1)*(2*m+1)*(2*m+1); fvecD wi(sz), xi(sz);
	image<T> dst(ex); dst = 0.0f; allSites(s,0) { unsigned int n = 0;			// for all sites s
		if (verbose && s.x == 0 && s.y == 0) { printf("%u\r", s.z); fflush(stdout); };
		if (mn(s) != 0.0f && var(s) != 0.0f) {					// do not work on background or constant regions
			boxIterator bx(ex,0,s,m); uvec3 t; while (bx(t)) {		// for all sites t in a box around s
				const float mr = mn(s)/mn(t);				// select voxels (Eq 7)...
				if (mu1 > mr || mr > 1.0f/mu1) continue;		// ...by means
				const float vr = var(s)/var(t);
				if (si1 > vr || vr > 1.0f/si1) continue;		// ...and by variances
				const float v = intensityDist(s, t, d);			// compute local intensity difference (Eq 5)
				const float w = std::exp(-v/(2.0f*beta*s2));		// compute weight (Eq 6)
				xi[n] = static_cast<float>(SQR((*this)(t))); wi[n] = w; n++; };
			float ws = 0.0f; for (unsigned int i = 0; i < n; i++) ws += wi[i];
			for (unsigned int i = 0; i < n; i++) wi[i] /= ws; };
		float xs = 0.0f; for (unsigned int i = 0; i < n; i++) xs += wi[i]*xi[i];
		xs -= 2*s2; dst(s) = n > 0? T(std::sqrt(xs)): (*this)(s); };
	return dst;
}

template <typename T>
image<T> image<T>::insertMatrix(const matD<T>& M, const unsigned int z)
//! inserts matrix M at plane z0.
{	if (nx < M.M || ny < M.N || z >= nz)
		throw rtException("image dimensions do not match with matrix");
	for (unsigned int y = 0; y < M.N; y++)
		for (unsigned int x = 0; x < M.M; x++) (*this)(x,y,z) = M(x,y);
	return *this;
}

template <typename T>
void readImages(std::vector<image<T>>& vim, const char* fname)
{	FILE* fp = openFile(fname,"r"); readImages(vim,fp); closeFile(fp);
}

template <typename T>
void readImages(std::vector<image<T>>& vim, FILE* fp)
{	bundleList bl; readFile(fp,bl);
	unsigned int n = 0; vim.resize(bl.size());
 	for (const auto& b : bl) { if (b.repn != repnType::image) continue;
		 vim[n].fromBundle(b); n++; };
	assert(n == bl.size());
}

template <typename T>
void saveImages(const std::vector<image<T>>& vim, FILE* fp, const repnType repn = repnType::unknown)
{	bundleList bl; for (const auto& im : vim) bl.push_back(im.toBundle("image",repn));	
	writeFile(fp,bl);
}

template <typename T>
void saveImages(const std::vector<image<T>>& vim, const char* fname, const repnType repn = repnType::unknown)
{	FILE* fp = openFile(fname,"w"); saveImages(vim,fp,repn); closeFile(fp);
}

template <typename T>
image<bool> image<T>::binarize(const T vmin, const T vmax) const
//! threshold a label image between vmin and vmax
{	image<bool> dst(nx,ny,nz);
	for (unsigned int i = 0; i < nel(); i++) dst(i) = (x[i] >= vmin && x[i] < vmax);
	return dst;
}

#endif

