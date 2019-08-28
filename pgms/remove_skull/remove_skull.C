

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>
static bimage getmask(const fimage& src, float open_dist)
//! computes a brain mask from T2-weighted image src.
{	const empiricalDist<float> em(src); const float th = em.quantile(0.5f);		// find threshold for 65% 
	bimage mask = src.binarize(th,std::numeric_limits<float>::max());	// classify and postprocess mask
	mask = src.classifyImage(256,3,mask).binarize(1.0f,3.0f).open(open_dist).selectBig(connectivity::c18).fillHoles();
	for (unsigned int z = mask.nz-1; z < mask.nz; z++)					// clear bottom-most slice
		for (unsigned int y = 0; y < mask.ny; y++)
			for (unsigned int x = 0; x < mask.nx; x++) mask(x,y,z) = false;
	if (mask.sum() == 0.0f) throw rtException("Mask image is empty");
	return mask;}

int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		option opts[] = { };			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);	//Parsing cmdline input
		vfimage vin;
		readImages(vin,inf);
		printf("Image loaded\n");
		fimage b0; b0 = vin[vin.size()-1];
		bimage mask; mask = getmask(b0, 2.0f);
		mask.save(outf);
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
