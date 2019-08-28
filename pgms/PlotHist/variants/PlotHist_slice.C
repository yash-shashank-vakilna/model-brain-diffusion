/*
 *
 *Computes histogram inputs from smt images
 *
 */

#include <brian.h>
#include <dwi.h>
#include <math.h>
static vfimage slicing(const vfimage& fim, int z)
{
vfimage vslice(2);
uvec3 ex = fim[0].extent();  fimage slice(ex.x,ex.y,1);
ex.z=z;
for(unsigned int i=0;i<2;i++){
	vslice[i](ex.x,ex.y,1);
		for(unsigned int x=0;x<fim[i].nx;x++){
			for(unsigned int y=0;y<fim[i].ny;y++){
				uvec3 s;
				s.x=x;s.y=y;s.z=0;
				ex.x=x;ex.y=y;
				slice(s)=fim[i](ex);
				printf("inter:%d \n",100000*x+y);			
}}
vslice[i]=slice;}
return vslice;
}

static int GetLinIndex( const float max_ele, const float ele, const int no_bins)
{   
int ind = floor(ele/max_ele*no_bins);
return ind;
}

static void writeText(fimage hist,FILE *file){
    	uvec3 ex = hist.extent(),s;
    	fimage outfile(ex);
	fwIterator fi(ex,0);
    	while(fi(s)){
	    	if (hist(s)==0) {outfile(s)=0;}
	    	else {outfile(s) = log(hist(s));}
	    	}
	outfile.save(file);
    	
}

static void populateHist(fimage& hist, const ivecD& ind)
{
    hist(ind(0),ind(1),0)++;
}

static int GetIndex(const float min_ele, const float max_ele, const float ele, const int no_bins)
{   
    if(ele < 0.000001){return 0;}
    float b[2]; //lin = b[0]*exp(b[1]*log)
    b[1]=log(max_ele/min_ele)/no_bins;
    b[0]=min_ele;
    //b[1]=1;b[0]=0;
    int ind = floor(log(ele/b[0])/b[1]);
    return (ind);
}



fimage PlotHist(const vfimage& fim, const int no_bins)
//Plots 2-D histogram based on no. of bins
//outputs in fimage
//only populateHist dim dependent
{	fimage hist(no_bins,no_bins,1); hist=0.0f; 		//intializing histogram
	int no_elements_im = fim[0].nel(); int no_dim = 2;	
        printf("no of elements: %d \n",no_elements_im);
	ivecD ind(no_dim); ind = 1;
	fvecD max_ele(no_dim),min_ele(no_dim);			//dim indepedndent max,min
    for(int i=0;i<no_dim;i++)					//populating max,min
    	{
	
    	    max_ele[i] = std::numeric_limits<float>::min();	//initializing max to minimum possible
	    min_ele[i] = std::numeric_limits<float>::max();	//initializing max to minimum possible
    	    uvec3 ex = fim[i].extent(),s; 			
    	    fwIterator fi(ex,0);
    	    while (fi(s)) { const float v = fim[i](s);		//finding max & min ele
		if (v < 0.000001) continue;
    	        else if (max_ele[i] < v) { max_ele[i] = v;}
    	        else if (min_ele[i] > v) { min_ele[i] = v;}
    	    }
    	}	
    uvec3 ex = fim[0].extent(),s; 
    fwIterator fi(ex,0);
    while(fi(s)){
            ind[0]=GetLinIndex(max_ele[0],fim[0](s),no_bins);
	    ind[1]=GetIndex(min_ele[1],max_ele[1],fim[1](s),no_bins);//passes appropriate data to Getindex and stores 
	
        populateHist(hist,ind); //increases the counter of hist specified by ind 
    }
    return(hist);
}


	
int main(int argc, char** argv)
{	int status = 0;
	try {	FILE *inf = nullptr, *outf = nullptr;
		option opts[] = { };			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);		//Parsing cmdline input
		vfimage fim,vslice(2); int no_bins=512 ;	
		readImages(fim,inf);						//read from smt file
		vslice = slicing(fim,60);
		fimage hist = PlotHist(vslice, no_bins);    			//Populates Histogram		
		writeText(hist,outf);   //writes the matrix in a text file specified by outf which can be used in matlab
		closeFile(inf); closeFile(outf); }
	catch (...) { brianException(); status = 2; };
	return status;
}
