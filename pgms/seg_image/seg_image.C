#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>
#include <spharmDWI.h>
#include <lregImage.h>
#include <pointset.h>

#define allDirs(i)	(unsigned int i = 0; i < nd; i++)
#define allVectors(i)	(unsigned int i = 0; i < nv; i++)

template <typename T>
static vec3<T> toVec3(const std::string& s)
 {   std::stringstream ss(s); vec3<T> a;
     ss >> a.x >> a.y >> a.z; return a;
 }

template<typename T> class new_mvNormalDist {
     static_assert(std::is_floating_point<T>::value, "mvNormalDist requires real float types");
 
     const vecD<T> mu;           
     const matD<T> Sigma;            
     const unsigned int N;           
     const matD<T> S;            
 
     matD<T> initS() const
         { matD<T> D(N,N); const auto V = sev(Sigma); D = T(0);
           for (unsigned int i = 0; i < N; i++) D(i,i) = T(std::sqrt(V.d(i)));
           return V.U*D; }
 public:
     new_mvNormalDist<T>(const vecD<T>& _mu, const matD<T>& _Sigma)          
         : mu(_mu), Sigma(_Sigma), N(mu.N), S(initS())
         { if (Sigma.M != N || Sigma.N != N)
             throw optException("mvNormalDist: Invalid argument"); }
     vecD<T> operator()() const                          
         { vecD<T> rn(N); for (unsigned int i = 0; i < N; i++)
             rn[i] = normalDist<T>::sample();
           return mu+S*rn; }
     T   lpdf(const vecD<T>& x) const                        
         { if (x.N != N) throw optException("mvNormalDist::lpdf: Invalid argument");
           const vecD<T> d = x-mu; const T t = T(0.5)*dot(d,inv(Sigma)*d);
           const T c = T(0.5)*(N*std::log(2.0*M_PI)+std::log(det(Sigma)));
           return -t-c; }
     T   pdf(const vecD<T>& x) const                     
         { return std::exp(lpdf(x)); }
     vecD<T> mean() const                                
         { return mu; }
     matD<T> var() const                             
         { return Sigma; }
     T   entropy() const                             
         { const T t = std::log(det(Sigma));
           return T(0.5)*(N+N*std::log(T(2)*M_PI)+t); }
 };
 


/********Pre-built gmmVB from Brian 3.1 starts**************************/
template<typename T> class new_gmmVB {
     static_assert(std::is_floating_point<T>::value, "new_gmmVB requires real float types");
 
     const T eps = T(1e-8);
     const std::vector<vecD<T>>& data;   
     const unsigned int N;           
     const unsigned int D;           
     unsigned int nc;            
     T   alpha0, beta0, nu0;
     vecD<T> alpha, beta, nu;
     vecD<T> mu0;
     matD<T> Wi0;
     std::vector<vecD<T>> mu, xm;
     std::vector<matD<T>> S, W;
     vecD<T> pi, lambda, nk;
     matD<T> R;
     const T removalCriterion;       
 
     T   lnDirichletC(const vecD<T>& alpha, const unsigned int nc)       
         { T as = T(0), cn = T(0);
           for allClasses(c) { as += alpha(c); cn -= std::lgamma(alpha(c)); }
           return cn+std::lgamma(as); }
     T   lnWishartB(const matD<T>& W, const T n)                 
         { unsigned int p = W.M; T b = T(-0.25)*T(std::log(M_PI))*p*(p-T(1));
           for (unsigned int i = 0; i < p; i++) b -= std::lgamma(T(0.5)*(n-i));  // \Gamma_p(n/2)
           b -= T(0.5)*n*std::log(det(W));
           b += T(0.5)*n*p*std::log(T(2)); return b; }               // |W|^(n/2) * 2^((n*p)/2)
     void    computeNk()
         { for allClasses(c) { T s = T(0); for allX(n) { s += R(n,c); };
             nk[c] = s; } }                          // compute Nk from r
     void    init();
     void    computeLambda();
     bool    removeClasses();
     void    eStep();
     void    mStep();
     T   costFunction();
 public:
     new_gmmVB(const std::vector<vecD<T>>& _data, const unsigned int _nc, const T rc)
          : data(_data), N(ITOU(data.size())), D(data[0].N), nc(_nc), alpha0(T(1)),
         beta0(T(1)), nu0(D), alpha(nc), beta(nc), nu(nc), mu0(D), Wi0(), mu(nc),
         xm(nc), S(nc), W(nc), pi(nc), lambda(nc), nk(nc), R(N,nc), removalCriterion(rc)
         { init();
           computeLambda(); }
     std::vector<std::pair<new_mvNormalDist<T>,T>> work(const unsigned int nit = 100, const bool verbose = false);
     uvecD   classify() const;
 };
 
 template <typename T>
 void new_gmmVB<T>::init()
 {   R = T(1)/nc; mu0 = T(0); nu = D;                        // init responsibilities
     uniformDist<float> ud; for allClasses(c) mu[c] = data[FTOU(ud()*N)];        // select random points as mean
     Wi0.resize(D,D); Wi0 = T(0); Wi0.addToD(T(4)/D); for allClasses(c) W[c] = Wi0;  // init precision matrices
     Wi0 = T(0); Wi0.addToD(D/T(4)); alpha = T(1); beta = T(10);             // init helper variables
 }
 
 template <typename T>
 void new_gmmVB<T>::computeLambda()
 {   T as = T(0); for allClasses(c) as += alpha(c); 
     for allClasses(c) pi[c] = std::exp(digamma(alpha[c])-digamma(as));      // compute pi
     for allClasses(c) { T t = T(0);
         for allDim(d) t += digamma(T(0.5)*(nu(c)-d));
         lambda[c] = std::exp(t+D*std::log(T(2.0))+std::log(det(W[c]))); };  // compute lambda
 }
 
 template <typename T>
 bool new_gmmVB<T>::removeClasses()
 {   unsigned int l = 0, r = 0;
     for allClasses(c) { if (nk(c) < removalCriterion*N) { r++; continue; };     // remove this column
         if (l != c) { for allX(n) R(n,l) = R(n,c); }; l++; };           // else move column
     assert(nc == l+r); nc = l; return r != 0;
 }
 
 template <typename T>
 void new_gmmVB<T>::mStep()
 {   computeNk(); if (removeClasses()) computeNk();                  // compute Nk from r
     for allClasses(c) { vecD<T> s(D); s = T(0);                 // compute xm from Nk & r
         for allX(n) s += R(n,c)*data[n];
         xm[c] = s/nk(c); };
     for allClasses(c) { matD<T> s(D,D); s = T(0);                   // compute S from xm & Nk
         for allX(n) { vecD<T> d = data[n]-xm[c]; s += R(n,c)*outer(d,d); };
         S[c] = s/nk[c]; };
     for allClasses(c) alpha[c] = alpha0+nk(c);                  // compute alpha from Nk
     for allClasses(c) beta[c] = beta0+nk(c);                    // compute beta from Nk
     for allClasses(c) mu[c] = (beta0*mu0+nk(c)*xm[c])/beta(c);          // compute mu from xm & beta
     for allClasses(c) nu[c] = nu0+nk(c);                        // compute nu from Nk
     for allClasses(c) { vecD<T> d = xm[c]-mu0; matD<T> v = outer(d,d);      // compute W from S, xm & Nk
         v = Wi0+nk[c]*S[c]+v*((beta0*nk(c))/(beta0+nk(c))); W[c] = inv(v); };
 }
 
 template <typename T>
 void new_gmmVB<T>::eStep()
 {   for allX(n) { T s = T(0); bool nm = false;                  // compute responsibilites
         for allClasses(c) { const vecD<T> d = data[n]-mu[c];
             const T t = D/beta(c)+nu(c)*dot(d,W[c]*d);
             T r = pi(c)*std::sqrt(lambda(c))*std::exp(T(-0.5)*t);
             if (r == T(0)) r = T(1);
             s += r; R(n,c) = r; };
         for allClasses(c) { T r = R(n,c)/s;                 // normalization
             if (r < eps) { r = eps; nm = true; }; R(n,c) = r; };        // check for small values
         if (nm) { T t = T(0); for allClasses(c) t += R(n,c);            // re-normalize if a small value was found
             for allClasses(c) R(n,c) /= t; } };
 }
 
 template <typename T>
 uvecD new_gmmVB<T>::classify() const
 {   uvecD cl(N);
     for allX(n) { T dmin = std::numeric_limits<T>::max(); unsigned int cmin = ~0u;  // compute responsibilites
         for allClasses(c) { const vecD<T> d = data[n]-mu[c];
             const T t = D/beta(c)+nu(c)*dot(d,W[c]*d);
             if (t < dmin) { dmin = t; cmin = c; } };
         cl[n] = cmin; };
     return cl;          
 }
 
 template <typename T>
 T new_gmmVB<T>::costFunction()
 {   T cn = T(0); for allX(n) { for allClasses(c) cn += R(n,c)*std::log(R(n,c)); };  // Eq. 61
     for allClasses(c) cn += (alpha(c)-alpha0)*std::log(pi(c));          // Eq. 62
     cn += lnDirichletC(alpha,nc); vecD<T> al(nc); al = alpha0;
     cn -= lnDirichletC(al,nc);                          // Eq. 66
     for allClasses(c) cn += lnWishartB(W[c],nu(c));
     cn -= nc*lnWishartB(inv(Wi0),nu0);
     T t = T(0); for allClasses(c) { const vecD<T> d = mu[c]-mu0;            // Eq. 67
         t += D*(beta0/beta(c)-std::log(beta0/beta(c))-nu(c)-T(1));
         t += beta0*nu(c)*dot(d,W[c]*d);
         t += nu(c)*(Wi0*W[c]).trace()+(nu(c)-nu0)*std::log(lambda(c)); };
     cn += T(0.5)*t; t = T(0); for allClasses(c) { vecD<T> d = xm[c]-mu[c];      // Eq. 64
         T t1 = std::log(lambda(c))+T(2)*std::log(pi(c))-D/beta(c)-D*std::log(T(2)*M_PI);
         T t2 = (S[c]*W[c]).trace()+dot(d,W[c]*d);
             t += nk(c)*(t1-nu(c)*t2); };
     cn -= T(0.5)*t; return cn;  
 }
 
 template <typename T>
 std::vector<std::pair<new_mvNormalDist<T>,T>> new_gmmVB<T>::work(const unsigned int nit, const bool verbose)
 {   T prevCost = std::numeric_limits<T>::max();
     for (unsigned int it = 0; it < nit; it++) {
         eStep();
         mStep();
         computeLambda();
         const T cost = costFunction();
         if (std::abs(prevCost-cost) < eps) break;
         prevCost = cost;
						};
     if (verbose) {  };
     std::vector<std::pair<new_mvNormalDist<T>,T>> theta;
     for allClasses(c) theta.push_back({{mu[c],inv(W[c]*N)},nk[c]/N});
     return theta;
 }
 /********Pre-built gmmVB from Brian 3.1 ends**************************/
 /********Pre-built dwimage from Brian 3.1 starts**************************/
struct new_dwimage {
	vfimage& src;				//!< vector of (nd+1) source images
	const unsigned int nd;			//!< number of gradient directions
	const unsigned int nv;			//!< number of gradient directions
	const uvec3 ex;				//!< image extent
	const vfvec3 g;				//!< vector of gradient directions
	const fvecD b;				//!< vector of gradient strengths

	new_dwimage(vfimage& _src, const unsigned int flipCode)
		: src(_src), nd(ITOU(src.size()-1)), nv(src[nd].nel()), ex(src[nd].extent()),
		g(initGradients(flipCode)), b(initStrengths())
		{  for allDir(d) { for allVoxels(i) { 
			if (src[d](i) <= 0) src[d](i) = 0.01f; } } }			// fix negative values inside mask
	fvec3	gradientDirection(const unsigned int d) const				//! returns gradient direction from attribute.
		{ const std::string& s = src[d].at.lookup("direction");
		  if (s.empty()) 
		  	throw optException("dwImage::gradientDirection: Invalid or missing attribute");
		  return toVec3<float>(s); }
	float	gradientStrength(const unsigned int d) const				//! returns gradient strength from attribute.
		{ const std::string& s = src[d].at.lookup("gradient_strength");
		  if (s.empty()) 
		  	throw optException("dwImage::gradientStrength: Invalid or missing attribute");
		  return std::stof(s); }
	float	getDiffusionTime(const unsigned int d) const				//! returns diffusion time from attribute.
		{ const std::string& s = src[d].at.lookup("diffusion_time");
		  if (s.empty()) 
		  	throw optException("dwImage::getDiffusionTime: Invalid or missing attribute");
		  return std::stof(s); }
	float	getPulseWidth(const unsigned int d) const				//! returns gradient pulse width from attribute.
		{ const std::string& s = src[d].at.lookup("pulse_width");
		  if (s.empty())
		  	throw optException("dwImage::getPulseWidth: Invalid or missing attribute");
		  return std::stof(s);}
	vfvec3	initGradients(const unsigned int flipCode) const			//! returns a vector of all gradient vectors and optionally flips them.
		{ vfvec3 v(nd);
		  for allDir(d) { fvec3 gi = gradientDirection(d);
		  	if (flipCode == 1) gi.x = -gi.x;
		  	else if (flipCode == 2) gi.y = -gi.y;
		  	else if (flipCode == 3) gi.z = -gi.z;
		  	else if (flipCode == 4) { gi.x = -gi.x; gi.y = -gi.y; }
		  	else if (flipCode == 5) { gi.x = -gi.x; gi.z = -gi.z; }
		  	else if (flipCode == 6) { gi.y = -gi.y; gi.z = -gi.z; };
			v[d] = gi; };
		  return v; }
	fvecD	initStrengths() const							//! returns a vector of all gradient strengths.
		{ fvecD v(nd);
		  for allDir(d) v[d] = gradientStrength(d)/1000.0f;
		  return v; }
	fmatD	initBm() const								//! initializes the linear tensor estimation matrix.
		{ fmatD B(nd, 6);
		  for allDir(d) { 
			B(d,0) = b(d)*g[d].x*g[d].x; B(d,1) = 2.0f*b(d)*g[d].x*g[d].y;
			B(d,2) = 2.0f*b(d)*g[d].x*g[d].z; B(d,3) = b(d)*g[d].y*g[d].y;
			B(d,4) = 2.0f*b(d)*g[d].y*g[d].z; B(d,5) = b(d)*g[d].z*g[d].z; }
		  return B; }
	fvecD	adc(const fvec3& s) const						//! returns apparent diffusion constant at site s.
		{ fvecD sm(nd); sm = 0.0f; if (src[nd](s) <= 0.0f) return sm; 
		  const float ls0 = std::log(src[nd](s));
		  for allDir(d) sm[d] = (ls0-std::log(src[d](s)))/b(d);
		  return sm; }
	fvecD	data(const uvec3& s) const						//! returns normalized data at site s.
		{ fvecD sm(nd); float sc = src[nd](s); sc = sc > 0.0f? 1.0f/sc: 0.0f;
		  for allDir(d) sm[d] = src[d](s)*sc;
		  return sm; }
	fvecD	interpolateAt(const uvec3& p) const					//! returns normalized data at position p.
		{ fvecD sm(nd); float sc; src[nd].interpolateAt(p,sc);
		  sc = sc > 0.0f? 1.0f/sc: 0.0f;
		  for allDir(d) { float s; src[d].interpolateAt(p,s); sm[d] = s*sc; };
		  return sm; }
	std::vector<std::pair<float,unsigned int>> getShells(const float tol = 0.1f) const //! returns number of gradient shells and gradients.
		{ std::vector<std::pair<float,unsigned int>> sh; sh.push_back({b(0),1});
		  for (unsigned int i = 1; i < b.N; i++) { unsigned int j;
			for (j = 0; j < sh.size(); j++)
				if (std::abs(sh[j].first-b(i)) < tol) {			// do not distinguish small differences in b
					sh[j].second++; break; };			// found this shell, so count gradient in for this shell
			if (j == sh.size()) sh.push_back({b(i),1}); };			// not found, so add this shell
		  return sh; }
	float	fractionalAnisotropy(const fvecD& s) const				//! returns fractional anisotropy from eigenvalues.
		{ const float md = (s(0)+s(1)+s(2))/3.0f; 
		  const float nom = SQR(s(0)-md)+SQR(s(1)-md)+SQR(s(2)-md);
		  const float den = SQR(s(0))+SQR(s(1))+SQR(s(2));
		  const float fa = den > 0.1f? std::sqrt(1.5f*nom/den): 0.0f;
		  return CLAMP(fa,0.0f,1.0f); };
	float	faAt(const fmatD& Bi, const uvec3& s) const				//! returns fractional anisotropy at site s.
		{ const fvecD t = Bi*adc(s); fmatD D(3,3);				// compute tensor, b(d) is included in Bi.
		  D(0,0) = t(0); D(0,1) = t(1); D(0,2) = t(2);
		  D(1,0) = t(1); D(1,1) = t(3); D(1,2) = t(4);
		  D(2,0) = t(2); D(2,1) = t(4); D(2,2) = t(5);
		  const auto S = sev(D); return fractionalAnisotropy(S.d); }		// decompose and return FA	
	fimage	estimateFA(const bimage& mask) const					//! estimates FA in mask; returns FA image.
		{ assert(mask.extent() == ex); fimage fa(ex); fa = 0.0f;
		  const fmatD Bi = pinv(initBm());
		  for (unsigned int i = 0; i < nv; i++) {
			if (mask(i) == false) continue;					// for all voxels inside mask
			const float f = faAt(Bi,mask.site(i));				// compute fractional anisotropy
			fa(i) = f; };							// store in image fa
		  fa.copyAttributes(src[nd].at); return fa; }
	fimage	b0Image() const								//! returns the b0 image (last one by definition)
		{ return src[nd]; }
	bimage	getMask() const;
	bimage	getForeground() const;							//! returns foreground mask
	limage	segmentCompartments(const std::string& rname = std::string()) const;
};

bimage new_dwimage::getForeground() const
//! returns foreground mask; if segmentation fails, returns "all in" mask.
{	bimage tmp = b0Image().classifyImage(256,2).binarize(1,2);			// select supra-threshold voxels
	for (const auto& si : src) {
		for (size_t i = 0; i < tmp.nel(); i++)
			if (tmp(i) && si(i) <= 0.0f) tmp(i) = false; };
	bimage mask = tmp.open(1.0f).selectBig(connectivity::c6).fillHoles();		// select biggest components and fill noise voxels
	if (mask.sum() == 0) { mask = true; };								// if mask is empty - flag all (for synthetic images)
	mask.copyAttributes(src[0].at); return mask;	
}

static uvecD identifyCompartments(std::vector<std::pair<new_mvNormalDist<double>,double>> theta)
//! maps multivariate Gaussian distributions theta onto compartment classes CSF, GM, WM
// NB: includes several heuristics about intensities & FA - specific to function classifyData() below.
{	const size_t nc = theta.size(); size_t ni = 0;
	uvecD cls(nc); cls = 3; dvecD b0(nc), fa(nc), t1(nc);
	for (size_t i = 0; i < nc; i++) {
		const dvecD& mu = theta[i].first.mean(); ni = mu.N;
		b0[i] = mu(0); fa[i] = mu(1); if (ni > 2) t1[i] = mu(2); };
	const auto b0r = rank(b0), far = rank(fa);
	cls[far(0)] = 1;													// CSF
	cls[far(1)] = 2;													// GM
	cls[far(2)] = 3; cls[far(3)] = 3;									// WM
	return cls;
}	
		  

limage new_dwimage::segmentCompartments(const std::string& rname) const
{	limage dst(ex); dst = 0; dst.copyAttributes(src[0].at);				// allocate destination image
	bimage mask = getForeground();										// mask out foreground
	if (mask.sum() == mask.nel()) { dst = 3; return dst; };				// mask is empty - flag all as WM
	vfimage src(2); src[0] = b0Image(); src[1] = estimateFA(mask);		// set b0 and fa image
	dst.at.remove("dt"); dst.at.remove("direction"); dst.at.remove("readout_time");
	dst.at.remove("gradient_strenDWIgth"); dst.at.remove("phase_encoding");
	if (!rname.empty()) { fimage ref; ref.read(rname.data());			// add T2 image if specified 
		const fvec3 tr = 0.0f, rt = 0.0f, sc = 1.0f;
		lregImage lt(ref,src[0]); lt.work(0, tr, rt, sc);				// register T2 to b0 using correlation
		fimage tref = lt.transform(ref); tref.copyAttributes(src[0].at);	// transform T2-weighted image into DWI space
		for (size_t i = 0; i < tref.nel(); i++) {						// clip DWI data by reference
			if (tref(i) < 1.0f) mask(i) = false; };
		src.push_back(tref);} 
	const size_t ns = std::min(static_cast<size_t>(mask.sum()/10.0f),50000ul);	// number of samples
	const size_t ni = src.size(); assert(ni <= 4);						// number of source channels
	const size_t nv = src[0].nel(); uniformDist<float> ud;				// collect ns samples from foreground
	vdvecD sample; dvecD vox(ni); size_t l = 0;
	for (size_t i = 0; i < ns; i++) {
		while (true) { l = FTOU(nv*ud()); if (mask(l)) break; };
		for (size_t i = 0; i < ni; i++) vox[i] = src[i](l);
		sample.push_back(vox); };
	const auto theta = new_gmmVB<double>(sample,4,0.01f).work(1000,false);		// classify as Gaussian distributions
	const size_t nc = theta.size(); assert(nc == 4);					// number of classes
	const auto cls = identifyCompartments(theta); dvecD p(nc);			// map classes onto compartments
	for (size_t i = 0; i < mask.nel(); i++) { if (mask(i) == false) continue;	// classify all foreground voxels
		for (size_t j = 0; j < ni; j++) vox[j] = src[j](i);				// sample at location i
		for (size_t j = 0; j < nc; j++) p[j] = theta[j].first.pdf(vox);	// compute probabilities
		dst(i) = cls(p.imax()); }										// map max. probability through classification
	return dst;
}
 /********Pre-built dwimage from Brian 3.1 starts**************************/
int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		vfimage vin; 
		option opts[] = { };			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);		//Parsing cmdline input
		readImages(vin,inf); 						//load dwimage
		printf("Image loaded\n");
		new_dwimage dwiI(vin,0);
		limage segmented_image = dwiI.segmentCompartments("102513_t1_al.v");
		segmented_image.save(outf);
		
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
