#ifndef BRIAN_H
#define BRIAN_H

/*
 *
 * brian.h: Common header file
 * BRIAN Software Package Version 3.0
 *
 * $Id: brian.h 509 2017-03-27 20:15:06Z kruggel $
 *
 * 0.10 (11/03/09): initial version
 * 0.20 (01/11/09): introduced new image format
 * 0.30 (31/12/12): released version 2.4
 * 0.40 (16/12/13): documented
 * v406 (28/09/16): bumped to version 3.0
 *
 */

/*! \file
    \brief Common header file and definitions.
*/

/*! \enum connectivity
    \brief Symbolic constants for connectivity in 2D and 3D image lattices.
*/
enum class connectivity {
	c4,					///< c4 edge connectivity, 2D
	c6, 					///< c6 face connectivity, 3D
	c6p, 					///< c6p connectivity, 2D
	c8, 					///< c8 edge and vertex connectivity, 2D
	c18,					///< c18 face and edge connectivity, 2D
	c26,					///< c26 face, edge and vertex connectivity, 2D
	ce					///< unknown connectivity (error)
};

extern connectivity compatibleConnectivity(const connectivity cn);

/*! \enum mcMethod
    \brief Symbolic constants for methods for multiple-comparison correction.
*/
enum class mcMethod {
	none,					///< unknown method (error)
	fw,					///< Friston-Worsley method
	bh					///< Benjamini-Hochberg procedure
};

/*! \enum planeDir
    \brief Symbolic constants for plane directions in 3D images.
*/
enum class planeDir {
	axial = 0,				///< axial plane (xy)
	coronal = 1,				///< coronal plane (xz)
	sagittal = 2,				///< sagittal plane (yz)
	nodir = 3				///< unknown direction (error)
};

/*! \enum wmObject
    \brief Symbolic constants for white matter objects.
*/
enum class wmObject {
	lh,					///< left hemisphere
	rh,					///< right hemisphere
	cb,					///< cerebrum
	brain					///< brain
};

/*! \enum problemType
    \brief Symbolic constants for problem types in elastic FE models.
*/
enum class problemType {
	staticElastic,				///< static linear elastic problem
	linearElastic,				///< dynamic linear elastic problem
	nonlinearElastic1,			///< dynamic non-linear elastic problem, 1-step integration
	nonlinearElastic2,			///< dynamic non-linear elastic problem, 2-step integration
	fsiElastic				///< dynamic fluid-solid non-linear elastic problem
};

/*! \enum renderMethod
    \brief Symbolic constants for displaying graphic primitives.
*/

enum class renderMethod {
	off = 0,				///< suppress display
	color = 1,				///< color-code voxel
	tensor = 2,				///< tensor glyph
	spharm = 3,				///< spharm gylph
	cylinder = 4,				///< cylinder glyph
	line = 5,				///< line glyph
	spline = 6,				///< spline glyph
	base = 7,				///< base glyph 
	surface = 8,				///< surface mesh
	volume = 9,				///< volume mesh 
	none = 10				///< no representation
};

#define SQR(a)		((a)*(a))
#define SIGN(a)		((a>0)-(a<0))
#define CLAMP(a,b,c)	std::min(std::max((a),(b)),(c))					// clamp a in [b,c]
#define ATOU(s)		(static_cast<unsigned int>(abs(atoi(s))))			// use wisely
#define ITOU(i)		(static_cast<unsigned int>(i))					// use wisely
#define FTOU(v)		(static_cast<unsigned int>(floor(v)))				// use wisely
#define LINELEN		1024

#define allMixtures(m)	(unsigned int m = 0; m < nm; m++)
#define allClasses(c)	(unsigned int c = 0; c < nc; c++)
#define allVoxels(i)	(unsigned int i = 0; i < nv; i++)
#define allBins(i)	(unsigned int i = 0; i < nbins; i++)
#define allSites(s, b)	fwIterator fi(ex, b); uvec3 s; \
			while (fi(s))

//----- C includes start here
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <errno.h>
#include <string.h>
//----- C++ includes start here
#include <stdexcept>
#include <exception>
#include <limits>
#include <typeinfo>
#include <complex>
#include <iterator>
#include <vector>
#include <list>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <algorithm>
//----- brian includes start here
#include <brianException.h>			// basic definitions
#include <cmdline.h>
#include <byteStream.h>
#include <file.h>
#include <numbers.h>				// math classes
#include <vector.h>
#include <matD.h>
#include <optimizer.h>
#include <quaternion.h>
#include <matS.h>
#include <image.h>				// image
#include <iterator.h>
#include <imageFuncs.h>
#include <functions.h>
#include <distributions.h>
#include <spline.h>
#include <imagePyramid.h>
#include <histo.h>
#include <shape.h>
#include <dwi.h>
#include <geomTests.h>				// graphs
#include <node.h>
#include <vertex.h>
#include <descriptor.h>
#include <graph.h>
#include <primitive.h>
#include <kdTree.h>
#include <mesh.h>
#include <triangleMesh.h>
#include <sphereMesh.h>
#include <spherePyramid.h>
#include <volumeMesh.h>				// volumes

// non-member entry points
extern fvecD	kmeansBinned(const fvecD& hs, const unsigned int nc);
extern vparam	estimateGMM(const fvecD& hs, const unsigned int nc, const bool verbose = false);
extern void	readHeader(FILE* fp, bundleList& list);
extern void	readMeshes(vmesh& im, const char *fname);
extern void	readMeshes(vmesh& im, FILE *fp);
extern void	saveMeshes(vmesh& im, const char *fname);
extern void	saveMeshes(vmesh& im, FILE *fp);
extern fimage	nlregImageSym(const fimage& src, const fimage& ref, const fimage& pay, 
			const float sigma);
extern void	nlregImageSym(const fimage& src, const fimage& ref, const float sigma, 
			FILE *outf);
extern limage	nlregLabelImage(const fimage& src, const fimage& ref, const limage& lbl,
			const float sigma);
extern fimage	nlregImageDemon(const fimage& src, const fimage& ref, const fimage& pay,
			const float sigma);
extern fimage	nlregImageDemon(const fimage& src, const fimage& ref, const float sigma);
extern fvec3image nlregImageDemonField(const fimage& src, const fimage& ref, const float sigma);
extern fvec3image nlregImageSymField(const fimage& src, const fimage& ref, const float sigma);
extern fimage	nlregJacobian(const fimage& src, const fimage& ref, const float sigma);
extern fimage	linregImage(const fimage& src, const fimage& ref, const fimage& pay, 
			const unsigned int plan, const bool scf, const fvec3& tr, const fvec3& rt,
			const fvec3& sc, const bool pre);
extern limage	linregLabelImage(const fimage& src, const fimage& ref, const limage& pay, 
			const unsigned int plan, const bool scf, const fvec3& tr, const fvec3& rt,
			const fvec3& sc, const bool pre);
extern fmat3	linregImage(const fimage& src, const fimage& ref, const unsigned int plan);
extern fimage	segmentImageFantasm(vfimage& cls, const fimage& src,
			const unsigned int nc, const float beta, const fimage& ref,
			const bool fixCSF);
extern fimage	segmentImageMRF(vfimage& cls, const fimage& src,
			const unsigned int nc, const float sigma, const float beta);
extern bool	segmentImageMC(vfimage& cls, const vfimage& src,
			const unsigned int nc, const float beta);
extern fimage	segmentImageWatershed(const fimage& src, const float rgf, 
			const float imax, const unsigned int size, const float back,
			const float inc);
extern bimage	segmentImageGM(const bimage& src, vfimage& cls, const float bk, const float ba,
			const float thr, const float beta, const char *mf, const bool thk);
extern void	alignDWI(const argVector& inf, const fimage& iso, const char* fname,
			const planeDir pl, const bool reg, const bool verbose);
extern fimage	alignImage(const fimage& src, const fimage& ref, const bool align);
extern bimage	stapleImages(const vbimage& src, const bool verbose);
extern void	evaluateAALRegions(const fimage& src, const vfimage& cls, 
			const unsigned int k, const fimage& ref, const limage& aal, FILE *outf);
extern void	mapAALRegions(triangleMesh& m, const char* fname, const limage& aal);
extern bimage	lesionWM(const vfimage& src, const bimage& m, const float thr);
extern bimage	lesionRG(const vfimage& src, const bimage& m, 
			const uvec3& seed, const unsigned int ov, const unsigned int lim);
extern bimage	lesionSYM(const vfimage& src, const bimage& m,
			const unsigned int win, const float sigma, const float thr, const uvec3& seed);
extern void	solveElasticFEModel(volumeMesh& m, problemType t, const char *bc,
			const bool verbose = false);
extern bimage	genusZeroRG(const bimage& src, const vfimage& cls, const bool verbose = false);
extern bimage	genusZeroLS(const bimage& src, const fimage& wmp, const float thr,
			const float alpha, const bool verbose = false);
extern bimage	genusZeroMSLS(const fimage& src, const float thr, const float alpha,
			const bool verbose = false);
extern descGraph* estimateTensorDWI(const dwImage& src, const unsigned int method,
			const float kappa, const char* fname);
extern descGraph* estimateSpharmDWI(const dwImage& src, const unsigned int deg, const bool sh);
extern descGraph* estimateFiberDWI(const dwImage& src, const unsigned int deg, const bool sh);
extern descGraph* estimateWatsonDWI(const dwImage& src);
extern descGraph* estimateZCDDWI(const dwImage& src);
extern fimage	correctImageN3(const fimage& src, const float lambda);
extern fimage	correctImageHisto(const fimage& src, const unsigned int wd, const float sigma);
extern fmat3	linregSphere(const sphereMesh& src, const sphereMesh& ref, const unsigned int plan, 
			const bool label, const bool verbose = false);
extern vfvec3	nlregSphere(const sphereMesh& src, const sphereMesh& ref, const bool verbose = false);
extern bool	findStreamlines(tensorGraph& tg, streamGraph& sg, const float wt,
			const float len, const float falim);
extern bimage	fillCentralWM(vfimage& cls, const float thr, const wmObject obj,
			const bimage& vent, const limage& lobe);
extern void	lmImages(vfimage& dst, FILE *inp, const char *rname, const bimage& mask,
			const float sigma, const int mt);
extern void	lmSlices(vfimage& dst, FILE *inp, const char *rname, const bimage& mask);
extern void	lmemImages(vfimage& dst, FILE *inp, const char *rname, const char *uname,
			const bimage& mask, const float sigma, const int mt, const int inf);
extern void	lmemSlices(vfimage& dst, FILE *inp, const char *rname, const char *uname,
			const bimage& mask, const int inf);
extern uvecD	segmentBasins(triangleMesh& m, const bimage& wm, const float thr,
			const float lim, const char* id);
extern float	estimateNoise(const fimage& src);
extern fimage	filterMean(const fimage& src, const unsigned int d = 2);
extern fimage	filterVar(const fimage& src, const fimage& mn, const unsigned int d = 2);
extern fimage	filterNLM(const fimage& src, const unsigned int m = 5, const unsigned int d = 2, const bool verbose = false);
extern fimage	filterNLMRician(const fimage& src, const unsigned int m = 5, const unsigned int d = 2, const bool verbose = false);
extern fimage	filterLee(const fimage& src, const unsigned int w, const float d = 2);
extern vfimage	estimateSMT(const vfimage& src);
extern unsigned int fgTN(const bimage& bin, const connectivity cn);
extern unsigned int bgTN(const bimage& bin, const connectivity cn);
extern void	eddyCorrectDWI(vfimage& dst, const argVector& in, const fimage& topup, const planeDir pl);
extern fimage	topupCorrectDWI(const argVector& in);
extern void	correctMesh(triangleMesh& src, const bool verbose = false);
extern fvecD	standardMeans(FILE* inf, const unsigned int nc, const float bg, const bool verbose = false);
extern fimage	standardizeIntensities(const fimage& src, const float bg, const fvecD& st, const unsigned int nm = 1000);
extern fvecD	similarityIndices(const bimage& src, const bimage& ref);
#endif

