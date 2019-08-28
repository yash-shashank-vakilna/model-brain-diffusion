/*
 *
 *Computes histogram inputs from smt images
 *Output is an image csv file
 */

#include <brian.h>
#include <dwi.h>
#include <math.h>
static int GetLinIndex( const float max_ele, const float ele, const int no_bins)
//!linear indexing
{   
if(ele > max_ele) {return -1;}
int ind = floor(ele/max_ele*no_bins);
return ind;
}

static void writeText(fimage hist,FILE *file)
//!write image to csv file
{
    	uvec3 ex = hist.extent();
	for (unsigned int i=0;i<ex.y;i++){
	    	for (unsigned int j=0;j<ex.x;j++)
	    	{
	    		fprintf (file, "%.06f, ",hist(j,i,0));
	    	}
		fprintf (file,"\n");
}
}


static void populateHist(fimage& hist, const ivecD& ind)
//!populates the image histogram
{
    if((ind(0) > 0) && (ind(1)>1)){
    hist(ind(0),ind(1),0)++;
    }
}

static int GetLogIndex(const float min_ele, const float max_ele, const float ele, const int no_bins)
//!log binninng algorithm
{   
    if((ele < min_ele) || (ele > max_ele) ){return -1;}
    float b[2]; //lin = b[0]*exp(b[1]*log)
    b[1]=log(max_ele/min_ele)/no_bins;
    b[0]=min_ele;
    //printf("log-const: [%f,%f]\n",b[0],b[1]);
    int ind = floor(log(ele/b[0])/b[1]);
    return (ind);
}

fimage PlotHist(const vfimage& fim, const int no_bins)
//Plots 2-D histogram based on no. of bins
//outputs in csv
//only populateHist dim dependent
{	fimage hist(no_bins,no_bins,1); hist=0.0f; 		//intializing histogram
	int no_elements_im = fim[0].nel(); int no_dim = 2;	
	ivecD ind(no_dim); ind = 1;
	fvecD max_ele(no_dim),min_ele(no_dim);			//dim indepedndent max,min
   /* for(int i=0;i<no_dim;i++)					//populating max,min
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
    	}*/
    	min_ele[0]= 0.001;min_ele[1]=0.01; max_ele[0]=3.05; max_ele[1]=2.5;
    	//min_ele[0]= 0.98;min_ele[1]=0.3; max_ele[0]=1; max_ele[1]=0.5;
    uvec3 ex = fim[0].extent(),s; 
    fwIterator fi(ex,0);
    while(fi(s)){
            ind[0]=GetLinIndex(max_ele[0],fim[0](s),no_bins);
	    ind[1]=GetLogIndex(min_ele[1],max_ele[1],fim[1](s),no_bins);//passes appropriate data to Getindex and stores 
	
        populateHist(hist,ind); //increases the counter of hist specified by ind 
    }
    return(hist);
}

static std::string getSerialID(int argc, char** argv)
//! extract first serial numbers from inp string
{
		unsigned int len = strlen(argv[2]);
		std::string fname(argv[2],0,6);									//clus_fname = filename of the labeled image
		return fname;
}

	
int main(int argc, char** argv)
{	int status = 0;
	try {	FILE *inf = nullptr, *outf = nullptr;
		option opts[] = { };
		std::string serial=getSerialID(argc,argv);
		printf("graph plotted for: %s",serial.c_str());
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);		//Parsing cmdline input
		vfimage fim; int no_bins=512 ;		
		readImages(fim,inf);						//read from smt file
		for (unsigned int i = 0; i < fim.size(); i += 1)
		{
			for (unsigned int j= 0; j < fim[i].nel(); j += 1)
			{
				fim[i](j) = fim[i](j);	//856
			}
		}
		fimage hist = PlotHist(fim, no_bins);    			//Populates Histogram		
		writeText(hist,outf);   //writes the matrix in a text file specified by outf which can be used in matlab
		closeFile(inf); closeFile(outf); }
	catch (...) { brianException(); status = 2; };
	return status;
}
