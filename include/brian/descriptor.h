#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

/*
 *
 * descriptor.h: voxel-wise fiber parameters
 * BRIAN Software Package Version 3.0
 *
 * $Id: descriptor.h 434 2016-11-04 16:22:39Z frithjof $
 *
 * 0.10 (19/09/16): rewritten for BRIAN3
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Definitions for voxel-wise fiber parameters.
*/

//! Implements a fiber descriptor record.

struct fiberDesc {				// 5 floats
	float	pf;				//!< volume fraction
	float	kappa;				//!< concentration
	fvec3	mu;				//!< direction

	fiberDesc(const float _pf = 0.0f, const float _kappa = 0.0f, const fvec3& _mu = 0.0f)
		//! allocates a fiber descriptor for fraction pf, concentration kappa and direction mu.
		 : pf(_pf), kappa(_kappa), mu(_mu) { }
	fiberDesc(const float* v)
		//! allocates a fiber descriptor from a float array.
		 : pf(v[0]), kappa(v[1]), mu(v[2],v[3],v[4]) { }
	void	print()	const								//! prints the record.
		{ printf("pf %e kappa %e mu %e %e %e ", pf, kappa, mu.x, mu.y, mu.z); }
};

#define allDesc(i) (unsigned int i = 0; i < desc.N; i++)

//! Base class for all descriptors.

struct descriptor : public node {
   	fvecD	desc;				//!< fiber descriptors

	descriptor(const unsigned int n = 0)						//! allocates a descriptor with n parameters.
		: node(), desc(n) { desc = 0.0f; }
   	descriptor(const descriptor& d)							//! copies from descriptor d.
		: node(d), desc(d.desc) { }
   	descriptor(descriptor&& d)							//! moves from descriptor d.
		: node(std::move(d)), desc(std::move(d.desc)) { }
	descriptor& operator=(const descriptor& d)					//! assigns from descriptor d including links.
		{ if (this != &d) { node::operator=(d); desc = d.desc; }; return *this; }
	descriptor& operator=(descriptor&& d)						//! move assigns from descriptor d including links.
		{ assert(this != &d); node::operator=(d);
		  desc = std::move(d.desc); return *this; }
	void	read(is& in, const bool hasWeights = false)				//! reads a primitive from stream in.
		{ node::read(in,hasWeights);
		  unsigned int n = 0; in.read(n); desc.resize(n);
		  for allDesc(i) { float d = 0.0f; in.read(d); desc[i] = d; } }
	void	save(os& out, const bool hasWeights = false) const			//! saves a primitive to stream out.
		{ node::save(out,hasWeights); unsigned int n = desc.N; out.save(n);
		  for allDesc(i) { float d = desc(i); out.save(d); } }
	void	print(const bool hasWeights = false) const
		{ printf("["); for (unsigned int i = 0; i < desc.N-1; i++)
			printT(desc(i));
		  printf("%12.5e] : ", desc(desc.N-1)); node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::descriptor; }
	descriptor* clone() const							//! clones this instance.
		{ return new descriptor(*this); }
	virtual float asIntensity() const						//! returns any meaningful intensity for this descriptor.
		{ return 0.0; }
	virtual fvec3 asColor() const							//! returns any meaningful RGB value.
		{ return fvec3(0.0); }
	virtual fmatD asTensor() const							//! returns any meaningful 3x3 tensor.
		{ return fmatD(); }
	virtual fvecD asSpharm() const							//! returns spharm coefficients.
		{ return fvecD(); }
	virtual std::vector<fiberDesc> asFibers() const					//! returns any meaningful fiber descriptors.
		{ return std::vector<fiberDesc>(); }
	virtual renderMethod defaultMode() const					//! display hint.
		{ return renderMethod::off; }
};

//! Implements a streamline record.

struct streamDesc : public descriptor {		// 6 floats
	streamDesc()
		: descriptor(0)	{ }							// no components
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::stream; }
	descriptor* clone() const							//! clones this instance.
		{ return new descriptor(*this); }
	renderMethod defaultMode() const						//! display hint.
		{ return renderMethod::spline; }
};

//! Implements a tensor record.

struct tensorDesc : public descriptor {		// 6 floats
	tensorDesc()
		: descriptor(6)	{ }							// 6 tensor components xx, xy, xz, yy, yz, zz
	fmatD	asTensor() const							//! returns a 3x3 tensor matrix.
		{ fmatD D(3,3);
		  D(0,0) = desc(0); D(0,1) = desc(1); D(0,2) = desc(2);
		  D(1,0) = desc(1); D(1,1) = desc(3); D(1,2) = desc(4);
		  D(2,0) = desc(2); D(2,1) = desc(4); D(2,2) = desc(5);
		  return D; }
	fvecD	decompose(fmatD& R) const
		{ fmatD D = asTensor(); const auto S = sev(D);
		  R = S.U; return S.d/S.d.sum(); }
	fvec3	asColor() const								//! converts largest eigenvector to RGB.
		{ fmatD R; fvecD e = decompose(R);					// convert to tensor & decompose
		  fvec3 c(std::abs(R(0,2)),std::abs(R(1,2)),std::abs(R(2,2)));
		  return c.normalize()*e[2]; }
	fvecD	asSpharm() const							//! converts to spharm (l=2).
		{ fvecD cf(6);								// see Descoteaux, thesis p. 81
		  cf[0] = float((2.0*std::sqrt(M_PI)/3.0)*(desc(0)+desc(3)+desc(5)));
		  cf[1] = float((2.0*std::sqrt(M_PI/15.0))*(desc(0)-desc(3)));
		  cf[2] = float((4.0*std::sqrt(M_PI/15.0))*desc(2));
		  cf[3] = float((-2.0*std::sqrt(M_PI/45.0))*(desc(0)+desc(3)-2.0*desc(5)));
		  cf[4] = float((4.0*std::sqrt(M_PI/15.0))*desc(4));
		  cf[5] = float((4.0*std::sqrt(M_PI/15.0))*desc(1));
		  return cf; }
	std::vector<fiberDesc> asFibers() const						//! returns a fiber descriptor.
		{ fmatD R; fvecD e = decompose(R);
		  fiberDesc fd1(e[2], 1.0f, fvec3(R(0,2),R(1,2),R(2,2))*e[2]);
		  fiberDesc fd2(e[1], 1.0f, fvec3(R(0,1),R(1,1),R(2,1))*e[1]);
		  std::vector<fiberDesc> vfd; vfd.push_back(fd1); vfd.push_back(fd2); 
		  return vfd; }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::tensor; }
	tensorDesc* clone() const							//! clones this instance.
		{ return new tensorDesc(*this); }
	renderMethod defaultMode() const						//! display hint.
		{ return renderMethod::tensor; }
};

//! Implements a spherical harmonics decomposition record of order 4.

struct spharmDesc : public descriptor {
	spharmDesc(const unsigned int n = 0)
		: descriptor(n) { }							// 6 tensor components xx, xy, xz, yy, yz, zz
	fvecD	asSpharm() const							//! returns element vector.
		{ return desc; }
	unsigned int degree() const							//! returns order of SH decomposition with nc coefficients.
		{ return FTOU(std::sqrt(8*desc.N+1)-3)/2; }
	float	generalizedAnisotropy(const unsigned int n) const
		{ float c0 = SQR(desc(0)), cs = 0.0f;
		  for (unsigned int i = 0; i < n; i++) cs += SQR(desc(i));
		  return std::sqrt(1.0f-c0/cs); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::spharm; }
	spharmDesc* clone() const							//! clones this instance.
		{ return new spharmDesc(*this); }
	renderMethod defaultMode() const						//! display hint.
		{ return renderMethod::spharm; }
};

//! Implements a fiber mixture record.

struct fiberMixtureDesc : public descriptor {
	fiberMixtureDesc()								//! allocates an empty fiber mixture record.
		: descriptor(19) { }							// 19 floats
		// 0: isotropic fraction, 1: signal intensity, 2: precision
		// 3: number of fibers, 4-19 up to three fiber descriptors
	std::vector<fiberDesc> asFibers() const						//! returns a vector of fiber descriptors.
		{ std::vector<fiberDesc> vfd; const unsigned int nf = FTOU(desc(3));
		  for (unsigned int i = 0; i < nf; i++)
			vfd.push_back({desc.x+4+5*i});
		  return vfd; }
	void	print(const bool hasWeights = false) const
		{ printf("[ s0 %e si2 %e piso %e\n", desc(1), desc(2), desc(0));
		  const unsigned int nf = FTOU(desc(3));
		  for (unsigned int i = 0; i < nf; i++)
			fiberDesc{desc.x+4+5*i}.print();
		  printf("] : "); node::print(hasWeights); }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::fiberMixture; }
	fiberMixtureDesc* clone() const							//! clones this instance.
		{ return new fiberMixtureDesc(*this); }
	renderMethod defaultMode() const						//! display hint.
		{ return renderMethod::cylinder; }
};

//! Implements a record for the Zeppelin-Cylinder-Dot model.

struct fiberZCDDesc : public descriptor {						// 11 floats
	fiberZCDDesc()									//! allocates an empty ZCD record.
		: descriptor(11) { }							// 11 floats
	void	print(const bool hasWeights = false) const				//! prints the record.
		{ printf("[ s0 %e si2 %e piso %e\n", desc(0), desc(1), desc(10));
		  printf("ic %e rad %e dpar %e mu %e %e %e\n",
			desc(8), desc(4), desc(2), desc(5), desc(6), desc(7));
		  printf("ec %e dper %e ] : ", desc(9), desc(3));
		  node::print(hasWeights); }
	std::vector<fiberDesc> asFibers() const						//! returns a vector of fiber descriptors.
		{ std::vector<fiberDesc> vfd;
		  fvec3 mu{desc(5),desc(6),desc(7)};
		  if (desc(8)) vfd.push_back({desc(8),20.0f,mu});
		  return vfd; }
	nodeID	getID() const								//! identifies this class for serialization.
		{ return nodeID::fiberZCD; }
	fiberZCDDesc* clone() const							//! clones this instance.
		{ return new fiberZCDDesc(*this); }
	renderMethod defaultMode() const						//! display hint.
		{ return renderMethod::cylinder; }
};
#endif

