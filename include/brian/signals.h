#ifndef SIGNALS_H
#define SIGNALS_H

/*
 *
 * signals.h: signal processing and filtering  
 * BRIAN Software Package Version 3.0
 *
 * $Id: signals.h 504 2017-03-16 23:28:00Z kruggel $
 *
 * 0.10 (18/11/11): for BRIAN2 by FK
 * 0.30 (31/12/12): released version 2.4
 * 0.31 (07/04/13): acf revised
 * 0.40 (16/12/13): documented
 * v396 (15/09/16): mop-up and const correctness, trace removed
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions and classes for signal processing and filtering.
*/

/*! \enum spectrumType
    \brief Symbolic constants for spectrum types in signal processing.
*/
enum class spectrumType {
	relPower,				///< power relative to maximum
	absPower				///< absolute power
};

/*! \enum filterType
    \brief Symbolic constants for filter characteristics in signal processing.
*/
enum class filterType {
	lowpass,				///< lowpass filter
	highpass,				///< highpass filter
	bandpass,				///< bandpass filter
	bandstop,				///< bandstop (notch) filter
	allpass					///< allpass filter
};

using cfl = std::complex<float>;
using cdb = std::complex<double>;

extern const cdb bessel_poles[];

//! Implements a moving average filter

class maFilter { 
	const unsigned int len;			//!< length of moving average
public:
	maFilter(const unsigned int _len)						//! allocates a moving average filter of length len.
		 : len(_len+(_len%2? 0:1)) { }
	fvecD	operator()(const fvecD& b) const					//! applies filter to timeseries b.
		{ const unsigned int l2 = len/2; fvecD u(b.N+len), a(b.N);		// compose extended timeseries
		  if (l2 > b.N) throw optException("Filter length %u too wide", int(len));
		  for (unsigned int i = 0; i < l2; i++) u[i] = b(l2-i-1);		// mirror half length of filter at start
		  for (unsigned int i = 0; i < b.N; i++) u[i+l2] = b(i);		// copy timeseries
		  for (unsigned int i = 0; i < l2+1; i++) u[i+l2+b.N] = b(b.N-i-1);	// mirror half length of filter at end
		  for (unsigned int t = 0; t < b.N; t++) { float v = 0;			// convolve timeseries x[0..n+len] with filter
			for (unsigned int i = 0; i < len; i++) v += u[t+i];
			a[t] = v/len; };
		  return a; }
};

//! Implements a finite impulse response filter

class firFilter { 
	fvecD	cf;				//!< filter coefficients

	void	normalize()								//! normalize filter coefficients.
		{ float n = 0;
		  for (unsigned int i = 0; i < cf.N; i++) n += std::abs(cf[i]);
		  cf /= n; }
	double	hammingWindow(const double m)						//! implements a Hamming window of length m.
		{ return 0.54+0.46*std::cos((M_PI*2.0*m)/double(cf.N-1)); }
public:
	firFilter(const filterType ft, const unsigned int len,
			const float _f1, const float _f2)				//! allocates a FIR filter of length len and cutoff frequencies (f1,f2).
		 : cf(len+(len%2? 0:1)) { double f1 = 2*M_PI*_f1, f2 = 2*M_PI*_f2;
		  for (unsigned int i = 0; i < cf.N; i++) { double m = i-0.5*(cf.N-1), c;
		  	switch (ft) {
			case filterType::lowpass: 
				c = m == 0.0? f1/M_PI: std::sin(m*f1)/(m*M_PI);
				break;
			case filterType::highpass:
				c = m == 0.0? 1.0-(f2/M_PI): -std::sin(m*f2)/(m*M_PI);
				break;
			case filterType::bandpass:
				c = m == 0.0? (f2-f1)/M_PI: (std::sin(m*f2)-std::sin(m*f1))/(m*M_PI);
				break;
			case filterType::bandstop:
				c = m == 0.0? 1.0+(f1-f2)/M_PI: (std::sin(m*f1)-std::sin(m*f2))/(m*M_PI);
				break;
			default: throw optException("Invalid filter type %u", static_cast<int>(ft)); };
			cf[i] = float(c * hammingWindow(m)); } }
	fvecD	operator()(const fvecD& b) const					//! applies filter to timeseries b.
		{ const unsigned int l2 = cf.N/2; fvecD u(b.N+cf.N), a(b.N);		// compose extended timeseries
		  if (l2 > b.N) throw optException("Filter length %u too wide", int(l2));
		  for (unsigned int i = 0; i < l2; i++) u[i] = b(l2-i-1);		// mirror half length of filter at start
		  for (unsigned int i = 0; i < b.N; i++) u[i+l2] = b(i);		// copy timeseries
		  for (unsigned int i = 0; i < l2+1; i++) u[i+l2+b.N] = b(b.N-i-1);	// mirror half length of filter at end
		  for (unsigned int t = 0; t < b.N; t++) { float v = 0;			// convolve timeseries x[0..n+len] with filter
			for (unsigned int i = 0; i < cf.N; i++) v += cf(i)*u[t+i];
			a[t] = v; };
		  return a; }
};

//! Implements a infinite impulse response filter

class iirFilter {
	static const unsigned int maxpz = 512;	//!< maximum number of coefficients
protected:
	const double eps = 1e-10;		//!< lower limit for filter coefficients
	zvecD	spoles, zpoles;			//!< poles
	zvecD	szeros, zzeros;			//!< roots
	dvecD	xc, yc;				//!< filter coefficients
	double	alpha1;				//!< lower frequency
	double	alpha2;				//!< upper frequency
	double	gain;				//!< normalization factor
	unsigned int snp, znp;			//!< pole counts
	unsigned int snz, znz;			//!< root counts

	void	multin(const cdb w, const unsigned int npz, zvecD& c) const		//! multiply coefficients c by w.
		{ const cdb nw = -w;
		  for (unsigned int i = npz; i >= 1; i--) c[i] = (nw*c(i))+c(i-1);
		  c[0] = nw*c(0); }
	cdb	expj(const double th) const						//! returns complex std::exponential.
		{ return cdb(std::cos(th), std::sin(th)); }
	cdb	sqrt(const cdb z) const							//! returns complex square root.
		{ const double r = std::hypot(z.imag(), z.real());
		  cdb v = cdb(std::sqrt(0.5*(r+z.real())), std::sqrt(0.5*(r-z.real())));
		  return (z.imag() >= 0.0)? v: conj(v); }
	cdb	eval(const zvecD& c, const cdb z) const					//! evaluate polynomial in z, substituting for c.
		{ cdb s = 0;
		  for (unsigned int i = c.N-1; i < c.N; i--) s = (s*z)+c(i);
		  return s; }
	cdb	evaluate(const zvecD& top, const zvecD& bot, const cdb z) const		//! evaluate top and bottom polynomial.
		{ return eval(top,z)/eval(bot,z); }
	void	setPole(const cdb z)							//! save pole z.
		{ if (z.real() < 0.0) spoles[snp++] = z; }
	void	expand(const zvecD& pz, zvecD& c) const					//! expand polynomial.
		{ c = cdb(0.0); c[0] = 1.0;
		  for (unsigned int i = 0; i < c.N-1; i++) multin(pz(i), c.N-1, c);
		  for (unsigned int i = 0; i < c.N; i++) assert(std::abs(c[i].imag()) < eps); }
	void	expandpoly(const filterType t);
	void	normalize(const filterType t);
	iirFilter(const double low, const double high, const unsigned n = maxpz)	//! allocates a generic IIR filter.
		 : spoles(n), zpoles(n), szeros(n), zzeros(n),
		xc(n+1), yc(n+1), alpha1(low), alpha2(high), 
		gain(0), snp(0), znp(0), snz(0), znz(0) { }
	virtual ~iirFilter() { }
public:
	fvecD	operator()(const fvecD& b) const					//! applies filter to timeseries b.
		{ dvecD xv(znp+1), yv(znp+1); xv = 0.0; yv = 0.0; fvecD a(b.N);
		  for (unsigned int i = 0; i < b.N; i++) {
			for (unsigned int j = 0; j < znz; j++) {
				xv[j] = xv(j+1); yv[j] = yv(j+1); };
			xv[znz] = b(i)/gain; yv[znz] = 0.0;
			for (unsigned int j = 0; j <= znz; j++) yv[znz] += xc(j)*xv(j);
			for (unsigned int j = 0; j < znz; j++) yv[znz] += yc(j)*yv(j);
			a[i] = float(yv(znz)); }; return a; }
};

//! Implements a bandpass filter

class iirResonator : public iirFilter {
	static const unsigned int nit = 50;		//!< iteration limit
	cdb	reflect(const cdb z) const						//! reflects coefficient z.
		{ const double r = std::hypot(z.imag(),z.real()); return z/SQR(r); }
public:
	iirResonator(const filterType ft, const double q, const double low, const double high);
};

//! Implements a Chebyshev filter

class iirChebyshev : public iirFilter {
public:
	iirChebyshev(const filterType ft, const unsigned int ord, const double r, 
			const double low, const double high)				//! allocates a Chebyshev filter of order ord with frequencies (low,high).
		 : iirFilter(low, high)
		{ for (unsigned int i = 0; i < 2*ord; i++) {
			const double th = ord&1? (i*M_PI)/ord: (i+0.5)*M_PI/ord;
			setPole(expj(th)); }
		  const double rip = std::pow(10.0, r/10.0), eps = std::sqrt(rip-1.0);
		  const double y = std::asinh(1.0/eps)/ord; assert(y > 0.0);
		  for (unsigned int i = 0; i < snp; i++)
			spoles[i] = cdb(spoles[i].real()*std::sinh(y), spoles[i].imag()*std::cosh(y));
		  normalize(ft); expandpoly(ft); }
};

//! Implements a Bessel filter

class iirBessel : public iirFilter {
	static const cdb bpoles[];							// see table in signal.C
public:
	iirBessel(const filterType ft, const unsigned int ord,
			const double low, const double high)				//! allocates a Bessel filter of type ft, order ord with frequencies (low,high).
		 : iirFilter(low, high)
		{ unsigned int p = SQR(ord)/4; if (ord & 1) setPole(bpoles[p++]);
		  for (unsigned int i = 0; i < ord/2; i++, p++) {
			setPole(bpoles[p]); setPole(conj(bpoles[p])); }
		  normalize(ft); expandpoly(ft); }
};

//! Implements a Butterworth filter

class iirButterworth : public iirFilter {
public:
	iirButterworth(const filterType ft, const unsigned int ord,
			const double low, const double high)				//! allocates a Butterworth filter of type ft, order ord with frequencies (low,high).
		 : iirFilter(low, high)
		{ for (unsigned int i = 0; i < 2*ord; i++) {
			const double th = ord&1? i*M_PI/ord: (i+0.5)*M_PI/ord;
			setPole(expj(th)); }
		  normalize(ft); expandpoly(ft); }
};

#define allChannels(c)	(unsigned int c = 0; c < nc; c++)

//! Implements a set of timeseries vectors

class signal {
public:
	vfvecD	ch;				//!< vector of channels
	unsigned int nc;			//!< number of channels
	float	sampleInterval;			//!< timestep per sample
	attrList at;				//!< attributes
private:
	void	checkUnits()								//! checks if this is a signal.
		{ const char *s = at.lookup("sampleInterval");
		  sampleInterval = s? float(atof(s)): 1.0f; 
		  s = at.lookup("xAxisLabel"); if (s == 0) return;			// assume sample interval "as is"
		  if (strncmp(s, "ms", 2) == 0) sampleInterval *= 0.001f;		// convert to seconds
		  else if (strncmp(s, "mHz", 2) == 0) sampleInterval *= 0.001f; }	// convert to Hz
	float	cphase(const cfl c) const						//! computes complex phase.
		{ const float p = std::atan2(-c.imag(),c.real());
		  return p >= 0.0f? p: p+M_PI; }
	void	scalePower(const unsigned int n, float* a) const			//! scales power spectrum to std::log max.
		{ float m = 0.0f;
		  for (unsigned int t = 0; t < n; t++) m = std::max(m,a[t]);		// find absolute maximum
		  if (m == 0.0f) return;
		  for (unsigned int t = 0; t < n; t++) { const float s = a[t]/m;
			a[t] = s > 0.0f? 10.0f*std::log10(s): -100.0f; } }		// scale std::log to maximum
	fvecD	powerSpectrum(const fvecD& s, const spectrumType st = spectrumType::absPower) const	//! returns power spectrum.
		{ const cvecD c = fft(s); fvecD d(c.N);
		  for (unsigned int t = 0; t < c.N; t++) d[t] = std::abs(c(t));
		  if (st == spectrumType::relPower) scalePower(c.N,d.x);
		  return d; }
	fvecD	phaseSpectrum(const fvecD& s) const					//! returns phase spectrum of channel s in radians.
		{ const cvecD c = fft(s); fvecD d(c.N);
		  for (unsigned int t = 0; t < c.N; t++) d[t] = cphase(c(t));
		  return d; }
	fvecD	autoCorrelation(const fvecD& s) const					//! computes auto-correlation of channel s.
		{ cvecD c = fft(s);
		  for (unsigned int t = 0; t < c.N; t++) c[t] = conj(c(t))*c(t);
		  return ifft(c); }
	fvecD	crossCorrelation(const fvecD& a, const fvecD& b) const			//! computes cross-correlation of channels a and b.
		{ assert(a.N == b.N); cvecD ca = fft(a), cb = fft(b);
		  for (unsigned int t = 0; t < ca.N; t++) ca[t] = conj(ca(t))*cb(t);
		  return ifft(ca); }
	fvecD	center(const fvecD& s) const						//! subtracts mean from channel.
		{ fvecD d = s; d -= s.mean(); return d; }
	fvecD	average(const fvecD& s, const char* stim, const unsigned int len,
			const unsigned int pre,	const char* cond = nullptr,
			const bool center = true) const					//! averages epochs in a channel.
		{ const unsigned int ll = 1024; char line[ll], tc[ll]; 
		  unsigned int cnt = 0, s0, ns = len+pre; 
		  fvecD av(ns); av = 0; FILE *fp = openFile(stim, "r"); 		// open condition file
		  while (fgets(line, ll, fp)) { if (*line == '#') continue;		// skip comment lines
			sscanf(line, "%u %s", &s0, tc); 				// get trigger start & condition
			if (cond && strcmp(cond, tc)) continue;				// check condition
			if (s0 < pre) continue;
			unsigned int off = FTOU(s0-pre); if (off+len >= s.N) continue;	// check range of timesteps to average
			for (unsigned int t = 0; t < ns; t++) av[t] += s(t+off);	// add to average
			cnt++; }
		  if (cnt) av /= cnt;
		  if (center) av -= av.mean();
		  closeFile(fp); return av; }
public:
	signal()									//! allocates an empty signal.
		 : ch(), nc(0), sampleInterval(0.001f), at() { }
	signal(const signal& s)								//! copies from signal s.
		 : ch(s.nc), nc(s.nc), sampleInterval(s.sampleInterval), at(s.at)
		{ for allChannels(c) ch[c] = s.ch[c]; }
	signal(const unsigned int _nc)							//! allocates a signal with nt channels.
		 : ch(_nc), nc(_nc), sampleInterval(0.001f), at()
		{ char buf[40]; sprintf(buf, "%u", nc); at.update("nChannels", buf); }
        fvecD	operator()(const unsigned int c) const					//! returns channel t.
		{ assert(c < nc); return ch[c]; }
        fvecD&	operator()(const unsigned int c)					//! returns reference to channel t.
		{ assert(c < nc); return ch[c]; }
	signal	histogram(const unsigned int nbins = 256) const				//! computes histograms of all channels.
		{ signal s; for allChannels(c) s.ch[c] = ch[c].histogram(nbins); return s; }
	signal	filterFIR(const filterType ft, const float low, const float high) const;
	signal	filterIIR(const filterType ft, const float low, const float high) const;
	signal	filterMA(const unsigned int len) const;
	signal	powerSpectrum(const spectrumType sp) const;
	signal	average(const char *stim, const unsigned int len, const unsigned int pre,
			const char *cond = 0, const bool center = true) const;
	signal	filterNotch(const float low, const float high, const unsigned int q = 4) const;
	signal	pca(const float cutoff) const;
	void	read(FILE *fp);
	void	read(const char *fname);
	void	save(FILE *fp) const;
	void	save(const char *fname) const;
	void	copyAttributes(const attrList &list)					//! copies attributes from list.
		{ at = list; }
};

fmatD asMatrix(const signal& S)
//! converts a signal into a matrix.
{	fmatD M(S.ch.size(),S.ch[0].N); 
	for (unsigned int i = 0; i < M.M; i++) M.setRow(i,S(i));
	return M;
}

signal asSignal(const fmatD& M)
//! converts a matrix into a signal.
{	signal S(M.M);
	for (unsigned int i = 0; i < M.M; i++) S.ch[i] = M.getRow(i);
	return S;
}
#endif
