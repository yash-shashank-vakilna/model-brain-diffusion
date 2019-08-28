

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>

int main(int argc, char** argv)
{	int status = 0;
	try {	
		/* *FILE *inf = nullptr, *outf = nullptr;
		fimage vin; 
		
		fimage vin;
		FILE *clus_f;
		clus_f = fopen("966975_clus_seg.v","r");
		vin.read(clus_f);
		printf("Image loaded\n");
		
		FILE *f = fopen("965367_clus_seg.v","r");
		fimage sum_im;
		sum_im.read(f);
		assert(sum_im.nel() == vin.nel());
		unsigned int nela = vin.nel();
		printf("%d\n",nela);
		for(unsigned int i = 0; i < nela; i++)
		{
			sum_im(i) = sum_im(i) + vin(i);
			printf("%\r",i);
		}
		sum_im.save("test.v");*/
		//sum_im.save(outf);
		fimage im1,im2;
		printf("reading: %s\n",argv[2]);
		im1.read(argv[2]);
		printf("%d %d %d %d\n",im1.nel(),im1.nx,im1.ny,im1.nz);
		printf("reading: %s\n",argv[4]);
		im2.read(argv[4]);
		unsigned int nela=im1.nel();
		fimage sum1(im1.extent());
		
		sum1=1.0f;
		printf("%d\n",sum1.nel());
		for(unsigned int i = 0; i < nela; i++)
		{
			sum1(i)=im1(i)+im2(i);
		}
		sum1.save(argv[4]);
		//closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
