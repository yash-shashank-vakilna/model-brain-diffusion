#ifndef PARTITION_H
#define PARTITION_H

/*
 *
 * partition.h: partition information for parallel vector and matrices
 * BRIAN Software Package Version 3.0
 *
 * $Id: partition.h 407 2016-10-02 21:36:46Z frithjof $
 *
 * 0.10 (28/11/11): for BRIAN2.2 by FK
 * 0.20 (16/12/13): documented
 * 0.30 (02/11/14): revised
 * 0.50 (17/11/14): released for BRIAN2.7
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Implements partitioning for parallel linear algebra.
*/

#define allProcs(r)	(unsigned int r = 0; r < nprocs; r++)
#define allNParts(i)	(unsigned int i = 0; i < nb.N; i++)

//! Implements a partitioning scheme for parallel vectors and matrices.

struct partition {
	unsigned int nprocs;			//!< number of partitions
	unsigned int myid;			//!< my process id
	bool	initialized;			//!< true if initialized
	ivecD	off;				//!< partition offsets
	uvecD	nb;				//!< number of elements
	uvecD	imi;				//!< import indices
	uvecD	imp;				//!< import pointers
	uvecD	exi;				//!< export indices
	uvecD	exp;				//!< export pointers

	partition()									//! allocates an empty partition.
		 : nprocs(mpiSize()), myid(mpiRank()), initialized(false), off(nprocs+1), nb(), imi(), imp(), exi(), exp()
		{ off = 0; }		
	partition(const unsigned int n)							//! allocates a partition for length n with granulation d.
		 : nprocs(mpiSize()), myid(mpiRank()), initialized(false), off(nprocs+1), nb(), imi(), imp(), exi(), exp()
		{ if (n < nprocs) throw rtException("More processes than rows");
		  const unsigned int sz = n/nprocs; off = 0; 				// determine size of a partition
		  for allProcs(i) off[i+1] = off[i]+sz;
		  off[nprocs] = n; }							// last one is slightly larger
	partition(const ivecD& _off)							//! allocates a partition for length n.
		 : nprocs(mpiSize()), myid(mpiRank()), initialized(false), off(_off), nb(), imi(), imp(), exi(), exp()
		{ if (off.N != nprocs+1) throw rtException("Partition does not match number of processors"); }
	unsigned int st() const								//! returns global start of partition.
		{ return off(myid); };
	unsigned int end() const							//! returns global end+1 of partition.
		{ return off(myid+1); };
	unsigned int len() const							//! returns length of partition.
		{ return off(myid+1)-off(myid); };
	unsigned int tlen() const							//! returns length + imports.
		{ return nprocs > 1? len()+imi.N: len(); };
	unsigned int id() const								//! returns my process id.		
		{ return myid; };
	unsigned int dim() const							//! returns total length.
		{ return off(nprocs); };
	bool	isLocal(unsigned int i) const						//! checks if element i is local.
		{ return i >= st() && i < end(); };
	bool	isRoot() const								//! checks if this is the root process.
		{ return myid == 0; };
	unsigned int partitionOf(const int c)						//! returns partition that contains column c.
		{ for (unsigned int i = 0; i < off.N-1; i++) if (c < off(i+1)) return i; 
		  return 0; }
	void	init(const uvecD& im);
		template <class T>
	void	sync(T* x) const;
};

void partition::init(const uvecD& im)
//! initializes a partition given a vector of import ids.
{	if (initialized) return; initialized = true; imi = im; if (nprocs == 1) return;	// do this just once
	uvecD imnb(nprocs), exnb(nprocs); imnb = 0; exnb = 0;
	for (unsigned int i = 0; i < imi.N; i++) imnb[partitionOf(imi(i))]++;		// find # of imports per partition 
	mpiSync(imnb.x, 1, exnb.x); unsigned int nbt = 0, n = 0;			// sync to receive exports
	for allProcs(i) { if (imnb[i] || exnb[i]) { nbt++; n += exnb[i]; } };		// find # of exports per partition
	nb.resize(nbt); imp.resize(nbt+1); exp.resize(nbt+1); exi.resize(n);		// allocate import and export pointers & indices			
	n = 0; imp[0] = 0; exp[0] = 0;							// set import and export arrays
	for allProcs(i) { if (imnb[i] || exnb[i]) { nb[n] = i; 
		imp[n+1] = imp[n]+imnb[i]; exp[n+1] = exp[n]+exnb[i]; n++; } };																
	uvecD rel(imi.N); for (unsigned int i = 0; i < imi.N; i++)			// for all imports
		rel[i] = imi[i]-off[partitionOf(imi(i))];				// find relative offsets into other partitions
	vecD<MPI_Request> rqe(nb.N); vecD<MPI_Status> ste(nb.N);			// per partition MPI request and status buffer1
	vecD<MPI_Request> rqi(nb.N); vecD<MPI_Status> sti(nb.N);			// per partition MPI request and status buffer1
	for allNParts(i) mpiSend(&rel[imp[i]], imp[i+1]-imp[i], nb[i], &rqe[i]);	// request import indices from other partitions	
	for allNParts(i) mpiRecv(&exi[exp[i]], exp[i+1]-exp[i], nb[i], &rqi[i]);	// receive export indices from other partitions
	mpiWait(rqe.N, rqe.x, ste.x); mpiWait(rqi.N, rqi.x, sti.x);			// wait for completed exports
}

template <class T>
void partition::sync(T* x) const
//! synchronizes x among all partitions.
{	assert(initialized); if (nprocs == 1) return;
	vecD<MPI_Request> rqe(nb.N); vecD<MPI_Status> ste(nb.N); vecD<T> we(exi.N);	// per partition MPI request and status buffer1
	vecD<MPI_Request> rqi(nb.N); vecD<MPI_Status> sti(nb.N); vecD<T> wi(imi.N);	// per partition MPI request and status buffer1
	for allNParts(i) {
		for (unsigned int t = exp(i); t < exp(i+1); t++) we[t] = x[exi(t)];	// fill send buffer
		mpiSend(&we[exp(i)], exp(i+1)-exp(i), nb(i), &rqe[i]); }		// send exports to other partitions
	for allNParts(i) mpiRecv(&wi[imp(i)], imp(i+1)-imp(i), nb(i), &rqi[i]);		// import values from other partitions
	mpiWait(rqi.N, rqi.x, sti.x);							// wait for completed imports
	unsigned int t = len(), l = imp(nb.N);
	for (unsigned int i = 0; i < l; i++) x[t++] = wi[i];				// set imported values in y
	mpiWait(rqe.N, rqe.x, ste.x);							// wait for completed exports
}

//! \relates partition
bool compatible(const partition& p, const partition& q)
//! returns true if matrices partitioned by A and B can be used in binary arithmetic operations.
{	if (&p == &q) return true; if (p.nprocs != q.nprocs) return false;
        if (p.initialized == false || q.initialized == false) return false;
        for (unsigned int i = 0; i <= p.nprocs; i++)
                if (p.off(i) != q.off(i)) return false;
	if (p.exi.N != q.exi.N) return false;
	for (unsigned int i = 0; i < p.exi.N; i++)
                if (p.exi(i) != q.exi(i)) return false;
	if (p.imi.N != q.imi.N) return false;
	for (unsigned int i = 0; i < p.imi.N; i++)
                if (p.imi(i) != q.imi(i)) return false;
	return true;
}

#undef allProcs
#undef allNParts

//! Maps the output of a partitioner to equations (and back).

struct ptmap {
	uvecD	eq;				//!< holds equation number for vertex id.
	uvecD	id;				//!< holds vertex id for equation number.
	ivecD	off;				//!< stores partition offsets.

	ptmap(const limage& pi)								//! allocates a map for partition vector& pi.
		: eq(pi.nel()), id(pi.nel()), off()
		{ unsigned int np = 0;							// determine number of partitions
		  for (unsigned int i = 0; i < id.N; i++) np = std::max(np,pi(i));
		  np++; off.resize(np+1); off[0] = 0; uvecD en(np); en = 0;		// allocate offset array
		  for (unsigned int i = 0; i < id.N; i++) en[pi(i,0,0)]++;		// count vertices per partition
		  for (unsigned int i = 0; i < np; i++) off[i+1] = off[i]+en[i];	// set partition offsets
		  for (unsigned int i = 0; i < np; i++) en[i] = off[i];			// copy offsets
		  for (unsigned int i = 0; i < id.N; i++) {				// for all vertices
			const unsigned int p = pi(i,0,0), e = en[p]++;			// get next eq# for this partition
			eq[i] = e; id[e] = i; } }					// set both maps
	ivecD	offsets() const								//! returns partition offsets.
		{ return off; }
	unsigned int toId(const unsigned int _eq) const					//! returns vertex id for equation eq.
		{ return id(_eq)+1; }							// correct for graph base
	unsigned int toEq(const unsigned int _id) const					//! returns equation eq for vertex id.
		{ return eq(_id-1); }							// correct for graph base
};

#endif
