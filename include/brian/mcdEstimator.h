#ifndef MCDESTIMATOR_H
#define MCDESTIMATOR_H

/*
 *
 * mcdEstimator.h: robust classification
 * BRIAN Software Package Version 3.0
 *
 * $Id: mcdEstimator.h 422 2016-10-19 00:27:19Z frithjof $
 *
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * 0.50 (14/10/14): revised
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements robust classification of multi-variate data using the MCD method.
*/

//! Helper class for an MCD estimator

struct item {
	unsigned int id;			//!< class id
	float	d; 				//!< sum of normalized distances

	item() : id(0), d(0) { }							//! allocates an empty item record.
	item(const unsigned int _id, const fvecD& mn, const fmatD& cv,
		const fvecD& v, const unsigned int nd)					//! allocates an item record i, given multi-variate mean mn and covariance matrix cv.
		 : id(_id), d(0) 
		{ fvecD u = v-mn; for (unsigned int i = 0; i < nd; i++) {
			float s = 0; for (unsigned int j = 0; j < nd; j++) s += u(j)*cv(j,i);
			d += s*u(i); } }
	bool	operator< (const item& rhs) const					//! comparison by increasing d.
		{ return (this->d < rhs.d); }
};

using vitem = std::vector<item>;

//! Collects information for a MCD estimator

class MCDEstimator
{	static const unsigned int maxit = 100; 	//!< maximum number of estimation iterations
	static const unsigned int nsol = 10;	//!< size of solution set
	const vfvecD& set;			//!< vector of observation tuples
	const unsigned int ns;			//!< number of observations
	const unsigned int nd;			//!< size of a tuple
	const float ratio;			//!< fraction of total samples
	const unsigned int setNum;		//!< total number of sets
	unsigned int setSize;			//!< size of each set
	fvecD	**mean;				//!< multi-variate means
	fmatD	**cov;				//!< multi-variate covariance matrices

	void	sample(fvecD& mn, fmatD& cv, const fvecD& v) const			//! sample from v, returns mean and covariance.
		{ for (unsigned int i = 0; i < nd; i++) { float vi = v(i); mn(i) += vi;
			for (unsigned int j = 0; j < nd; j++) cv(i,j) += vi*v(j); } }
	bool	normalize(fvecD& mn, fmatD& cv, const float n) const			//! normalizes mean and covariance.
		{ mn /= n; cv /= n; for (unsigned int i = 0; i < nd; i++) {
			for (unsigned int j = 0; j < nd; j++) cv(i,j) -= mn(i)*mn(j); };
		  return isPSD(cv); }
	bool	initialize(fvecD& mn, fmatD& cv, const vfvecD& s) const			//! determines mean and covariance of sample s.
		{ mn = 0; cv = 0; const unsigned int n = ITOU(s.size());
		  for (unsigned int i = 0; i < n; i++) sample(mn, cv, s[i]);
		  return normalize(mn, cv, n); }
	bool	estimate(fvecD& mn, fmatD& cv, const vfvecD& s)	const			//! estimates mean and covariance by sampling s.
		{ const unsigned int n = ITOU(s.size()), nr = FTOU(n*ratio);
		  const fmatD cvi = inv(cv); vitem ar(n); 
		  for (unsigned int i = 0; i < n; i++) ar[i] = item(i, mn, cvi, s[i], nd); 
		  std::sort(ar.begin(), ar.end()); mn = 0; cv = 0;
		  for (unsigned int i = 0; i < nr; i++) sample(mn, cv, s[ar[i].id]);
		  return normalize(mn, cv, nr); }
	void	randomSample(vfvecD& s, const unsigned int n)				//! draws a random sample s.
		{ for (unsigned int i = 0; i < n; i++) {
			const unsigned int r = FTOU(ns*drand48()); s[i] = set[r]; }; }
public:
	MCDEstimator(const vfvecD& _set, const unsigned int _nd,
		const unsigned int _sn, const unsigned int _ss, const float _r)		//! allocates an MCD estimator for sample set.
		 :  set(_set), ns(ITOU(set.size())), nd(_nd), ratio(_r),
		setNum(_sn), setSize(_ss), mean(0), cov(0)
		{ srand48(time(0)); mean = new fvecD* [setNum]; cov = new fmatD* [setNum];
		  for (unsigned int i = 0; i < setNum; i++) {
			mean[i] = new fvecD [nsol]; cov[i] = new fmatD [nsol]; }; }
	MCDEstimator(const MCDEstimator& m) = delete; 
	MCDEstimator(MCDEstimator&& m) = delete; 
	~MCDEstimator()
		{ for (unsigned int i = 0; i < setNum; i++) { delete [] mean[i]; delete [] cov[i]; }
		  delete [] mean; delete [] cov; }
	MCDEstimator& operator=(const MCDEstimator& m) = delete;
	MCDEstimator& operator=(MCDEstimator&& m) = delete; 
	void	work(fvecD& mn, fmatD& cv);

};

void MCDEstimator::work(fvecD& mn, fmatD& cv)
//! performs classification on sample, returns multi-variate mean and covariance matrix.
{	float d, od, md = 0; fvecD dt(nsol); mn = 0; cv = 0;
	for (unsigned int i = 0; i < setNum; i++) {     						
		vfvecD sample(setSize);	randomSample(sample, setSize);			// randomly draw a new sample from the population
        	dt = std::numeric_limits<float>::max();
        	for (unsigned int j = 0; j < maxit; j++) {				// repeat initalization&  estimation 100 times
			initialize(mn, cv, sample);
			if (estimate(mn, cv, sample) == false) continue;
			const float dj = det(cv);
			if (j < nsol) { dt[j] = dj; mean[i][j] = mn; cov[i][j] = cv; }	// record best NSOL solutions
			else { for (unsigned int k = 0; k < dt.N; k++)
					if (dj < dt[k]) { dt[k] = dj;
						mean[i][k] = mn; cov[i][k] = cv; break; } } } }; 
        dt = std::numeric_limits<float>::max();
	for (unsigned int i = 0; i < setNum; i++)  {     				// evolution in the large set
		vfvecD sample(setNum*setSize); randomSample(sample, setNum*setSize);
		for (unsigned int j = 0; j < nsol; j++)  {				// improve on NSOL solutions
			mn = mean[i][j]; cv = cov[i][j];
			if (estimate(mn, cv, sample) == false) continue;
			const float dj = det(cv);
			if (dj < dt[j]) { dt[j] = dj; mean[0][j] = mn; cov[0][j] = cv; } } };
	for (unsigned int i = 0; i < nsol; i++) {					// evolution in the full set
		vfvecD sample(ns); randomSample(sample, ns);
		d = std::numeric_limits<float>::max();
		for (unsigned int it = 0; it < maxit; it++) {
			if (estimate(mean[0][i], cov[0][i], sample) == false) continue;
			od = d; d = det(cov[0][i]); if (std::abs(d-od)/od < 0.0001f) break; }
        	if (isSmall<float>(md) || d < md) { mn = mean[0][i]; cv = cov[0][i]; md = d; } }; 
}

#undef NSOL
#undef MAXIT
#endif
