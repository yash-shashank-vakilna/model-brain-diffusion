

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>

template <class T>
static image<T> loadImage(int argc, char** argv, std::string append_name)
//!load another file from the same dir
//name chosen from first patient id of output file_name
//always put before c.parse
{
	image<T> im;
	std::string fname(argv[4],0,6);
	fname.append(append_name);	
	im.read(fname.c_str());
	return im;
}

int main(int argc, char** argv)
{	
//inp = mask file and out = smt file
	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		vfimage vin; 
		option opts[] = { };
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);			//Parsing cmdline input
		printf("Image loaded\n");
		bimage mask; mask.read(inf);
		readImages(vin,outf);
		for (unsigned int i = 0; i < mask.nel(); i++)
		{
			if (!mask(i)) {
				vin[0](i)=0; vin[1](i)=0;}
		}
		saveImages(vin,outf);
		
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
