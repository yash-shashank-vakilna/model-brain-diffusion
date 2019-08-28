

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>

static void extract_frm_segim(vfimage& smtimage, const limage& seg_image,
										const unsigned int& tissue_type)
{
	uvec3 s, ex = seg_image.extent();
	fwIterator fi(ex,0); size_t no_im = smtimage.size();
	while(fi(s)){
		if(seg_image(s) != static_cast<long int>(tissue_type))			//matching tissue types
		{			
			for (unsigned int i = 0; i < no_im; i++)
			{
				smtimage[i](s)=0;										//removing tissues which are not required
			}
		}
	}
}

int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		vfimage vin; 
		option opts[] = { };
		std::string clus_fname(argv[4],0,6), append_name = "_seg.v";	//clus_fname = filename of the labeled image
		clus_fname.append(append_name);			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);			//Parsing cmdline input
		readImages(vin,inf);
		printf("Image loaded\n");
		limage seg_image;seg_image.read(clus_fname.c_str());
		unsigned int tissue_type = 2;									//1= CSF, 2=GM, 3=WM
		extract_frm_segim(vin, seg_image, tissue_type); 				//extract tissue type
		saveImages(vin,outf);
		closeFile(inf);closeFile(outf);
		
		
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
