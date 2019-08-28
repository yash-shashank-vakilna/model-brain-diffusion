

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>

int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		vfimage vin; 
		option opts[] = { };			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);		//Parsing cmdline input
		readImages(vin,inf);
		printf("Image loaded\n");
		
		
		
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
