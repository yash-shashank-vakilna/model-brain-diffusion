

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>
#include <pointset.h>
#define max_sample 881651
#define log_transform 1
#define clus_type 1 //0: GMM, 1: k-means

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
     std::vector<std::pair<mvNormalDist<T>,T>> work(const unsigned int nit = 100, const bool verbose = false);
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
 std::vector<std::pair<mvNormalDist<T>,T>> new_gmmVB<T>::work(const unsigned int nit, const bool verbose)
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
     std::vector<std::pair<mvNormalDist<T>,T>> theta;
     for allClasses(c) theta.push_back({{mu[c],inv(W[c]*N)},nk[c]/N});
     return theta;
 }
 /********Pre-built gmmVB from Brian 3.1 ends**************************/
 static uvecD identifyCompartments(std::vector<fvecD> sample, uvecD clus_vec ,unsigned int nc)
{	
	fvecD FA_mean(nc);uvecD no_ele(nc), cls(nc); size_t sample_size = sample.size();
	FA_mean=0; no_ele=0;
	for (int i = 0; i < sample_size; i++)
	{
		FA_mean(clus_vec(i)) += sample[i](2);
		no_ele(clus_vec(i))++;
	}
	for(unsigned int i = 0; i < nc; i++)
	{
		FA_mean(i) = FA_mean(i)/no_ele(i);
	}
	const auto far = rank(FA_mean);
	cls[far(0)] = 1;													// CSF
	cls[far(1)] = 2;													// GM
	cls[far(2)] = 3; 
	cls[far(3)] = 4;									// WM
	cls[far(4)] = 4;
	return cls;
}

 static float compute_FA(float lp, float lt)
 {
	 float lm = (lp+2*lt)/3;
	 float FA = (lm==0)?0:std::sqrt(1.5*(std::pow(lp-lm,2)+2*std::pow(lt-lm,2))
				/(std::pow(lp,2)+2*std::pow(lt,2)));
	 return FA;
}

static void imageToVec(const vfimage& vim, std::vector< fvecD>& value_vec,
						 std::vector< unsigned int>& ind_vec)
//!Converts image to in input-format of new_gmm
{	
	unsigned int size_im = vim[0].nel();
	for (unsigned int nos = 0; nos < size_im; nos += 1)
	{	
		
		if (vim[0](nos)+vim[1](nos) == 0) 
		{
			
			continue;
		}
		fvecD sample(3);
		sample[0] = vim[0](nos); 
		if(vim[1](nos) == 0)
			{sample[1]=0;}
		else
		{
		sample[1] = log_transform?log(vim[1](nos)):vim[1](nos);
		}
		sample[2]=compute_FA(vim[0](nos),vim[1](nos));
		value_vec.push_back(sample);
		ind_vec.push_back(nos);
		
	}
}

uvecD gmm_cluster( const std::vector< fvecD>& vec, unsigned int nc)
//!Perform kmeans clustering 
{
	new_gmmVB<float> gmfunc(vec,nc,0.01);
	uvecD clus;clus = gmfunc.classify();
	return clus;
}

static void writeToTxt(FILE *file, const std::vector< fvecD>& value_vec, 
						uvecD clus_vec)
//!write image to csv file
{
	for (unsigned int i=0;i<value_vec.size();i++)
	{
		fprintf(file,"%f, %f, %d\n",
		value_vec[i](0), value_vec[i](1), clus_vec(i));
	}
}

static void populateImagespace(fimage& clusimage, const uvecD& clusvec, 
				const std::vector< unsigned int>& ind_vec, uvecD cls)
//!populate clusimage(voxels=clusters) using indices stored in ind_vec and cluster-classification in clus_vec
{
	for (unsigned int i = 0; i < clusvec.N; i += 1)
	{
		clusimage(ind_vec[i]) = cls(clusvec(i));
	}
}

static uvecD km_cluster( const std::vector< fvecD>& vec, unsigned int nc)
//!Perform kmeans clustering 
{
	kmeans<float> kmfunc(vec,nc);
	uvecD clus;
	clus = kmfunc.work();
	return clus;
}

int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr, *csvf = nullptr;
		vfimage vin; 
		option opts[] = { };	
		std::string clus_fname(argv[4],0,6), append_name = "_clus.csv";
		clus_fname.append(append_name);
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);		//Parsing cmdline input
		if(log_transform) printf("lt is log-transformed, to switch off change in #DEFINE\n");
		(clus_type)?printf("performing k-means\n"):printf("performing GMM\n");
		readImages(vin,inf); 
		printf("Image loaded\n");
		std::vector< fvecD> value_vec; 
		std::vector< unsigned int> ind_vec;
		
		imageToVec(vin, value_vec, ind_vec);				//convert to input format for kmeans
		unsigned int nc = 5;
		uvecD clus_vec; 
		clus_vec = clus_type?
		gmm_cluster(value_vec,nc):km_cluster(value_vec,nc); 		//perform clustering 0:gmm & 1:k-means
		uvecD cls = identifyCompartments( value_vec, clus_vec ,nc);
		csvf = fopen(clus_fname.data(),"w");
		//writeToTxt(csvf, value_vec, clus_vec);
		fimage clusimage(vin[0]);					//intialize clusimage
		populateImagespace(clusimage, clus_vec, ind_vec, cls);		//populating clusters in image_space
		clusimage.save(outf);
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
