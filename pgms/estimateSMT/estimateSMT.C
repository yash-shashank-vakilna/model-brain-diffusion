/*
 *
 *V2.0 Find Spherical Mean with multiple b-shells
 *
 */

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>
#define lFree 3.05
class LUT_class
{
	
	dvecD LUT;
	const float max_val = 3.2;
	const float LUT_lc = 1e-6;
	public:
	LUT_class()
	{
		//printf("LUT comstructor called\n");
		unsigned int no_ele = round(max_val/LUT_lc);
		dvecD temp_vec(no_ele);
		LUT(0)=2/std::sqrt(M_PI);
		for (unsigned int i = 1; i < no_ele; i += 1)
		{
			temp_vec(i) = erf(i*LUT_lc)/(i*LUT_lc);;
		}
		LUT = temp_vec;
	}
	
	double lookup(const double& x)
	{	
		unsigned int index = round(x/LUT_lc); 							//round to lc
		if (index >= (max_val/LUT_lc)) throw rtException("max LUT val reached");
		double val = LUT(index);
		return val;
	}
};

static LUT_class erf_lut;

static bimage getmask(const fimage& src)
//! computes a brain mask from T2-weighted image src.
{	const empiricalDist<float> em(src); const float th = em.quantile(0.75f);		// find threshold for 65% 
	bimage mask = src.binarize(th,std::numeric_limits<float>::max());	// classify and postprocess mask
	mask = src.classifyImage(256,3,mask).binarize(1.0f,3.0f).open(2.0f).selectBig(connectivity::c18).fillHoles();
	for (unsigned int z = mask.nz-1; z < mask.nz; z++)					// clear bottom-most slice
		for (unsigned int y = 0; y < mask.ny; y++)
			for (unsigned int x = 0; x < mask.nx; x++) mask(x,y,z) = false;
	if (mask.sum() == 0.0f) throw rtException("Mask image is empty");
	return mask;}


static dvecD modeledSignal(const dvecD& p, const dvecD& b_sh)
//!Calculates signal based on the parametric response model
			{ 
			  dvecD Egh(b_sh.N);
			  const double lp = std::max(p(0)*lFree,0.0);				// safeguard against p(0) = -0
			  const double lt = std::min(p(1)*lp,lp);					// safeguard against p(1) = 1 + eps
			  for (unsigned int i = 0; i < b_sh.N; i++) {				// see Eq. 8
				const double x = std::sqrt(b_sh(i)*(lp-lt));
				const double t = x == 0.0? 1.0:							// use lim erf(x)/x for x -> 0
						0.5*std::sqrt(M_PI)*erf_lut.lookup(x);
			  	 Egh[i] = std::exp(-b_sh(i)*lt)*t;}
			  return Egh; }

static dvecD opt (const dvecD& E, const dvecD& b_sh)
//!Optimize using monte-carlo method and return p(ratios)
{
	unsigned int noit = 1000,  dim = E.N; 
	dvecD p(dim), E_cal(dim), E_opt(dim), popt(dim);
	double max_SSQ = std::numeric_limits<double>::max(); 
	for (unsigned int iter = 0; iter < noit; iter += 1)
	{
		for (unsigned int j = 0; j < dim; j += 1)
		{
			
			p(j)= static_cast<double>(rand())/RAND_MAX;					//sample random number
		}
		E_cal = modeledSignal(p, b_sh);
		double error, SSQ=0; 
		for (unsigned int j = 0; j < dim; j += 1)						//compute error
		{
			error = SQR(E_cal(j)-E(j));
			SSQ += error;
		}											
		if (SSQ < max_SSQ)												//check errors and store minimum error parameter (p)
		{
			popt = p;
			max_SSQ = SSQ;
			
		}
	}
	return popt;

}

static vfimage getPar( const vfimage& meanIm, const dvecD& b_sh ){
//!Finds lp and lt by performing minimization based on parametric response model
	vfimage par(2);
	dvecD E_gen(b_sh.N), prange(2), popt(2); prange = 0.5 ;popt = 0.5;	//declare and intialize optimizer related variable
	uvec3 s,ex=meanIm[0].extent();
	fimage temp(ex); temp = 0;
	par[0]=temp;par[1]=temp;
	fwIterator fi(ex,0);
	unsigned int prog=0, nel=meanIm[0].nel();
	while(fi(s))
	{
	prog++;
	printf("%.2f%% par optimzed\r",(static_cast<float>(100*prog)/nel));
		for (unsigned int i = 0; i < E_gen.N; i += 1)
		{
			E_gen(i)= static_cast<double>(meanIm[i](s));				//loop through image
		} 	
		if (E_gen.sum()==0)
		{
			continue;
		}			
		popt = opt(E_gen,b_sh);											//pass intensities to the optimizer
		par[0](s) = static_cast<float>(popt(0)*lFree); par[1](s)=
						static_cast<float>(popt(1)*popt(0)*lFree); 		//store in par
		
	}
	return(par);
	}

static vfimage normalizeAverageImage(const vfimage& Vin, const std::vector< std::vector< unsigned int>>& sort_ind, const fimage& b0)
//!averages the image (performs SMT) and finds the ratio mean:bo
{
	vfimage meanIm;
	bimage mask = getmask(b0);
	//bimage mask(b0.extent()); mask=true;
	for (unsigned int i = 0; i < sort_ind.size(); i += 1)
	{
		fimage temp(Vin[i].nx,Vin[i].ny,Vin[i].nz);						//Initializing MeanIm 
		temp = 0.0f;
		meanIm.push_back(temp);
	}
	uvec3 ex=b0.extent(),s;
	fwIterator fi(ex,0);
	unsigned int no_b = sort_ind.size();
	unsigned int prog=0, nel=b0.nel();
	
	while(fi(s))														//looping over the image
	{
		bool art=false;	
		for (unsigned int i = 0; i < no_b; i += 1)						//looping over b-val
		{
		if(mask(s)){
		unsigned int size_sort_ind = sort_ind[i].size();
		unsigned int j;
			for ( j = 0; j < size_sort_ind; j += 1)						//looping over images of 1 b-val
			{	
				unsigned int im_index =sort_ind[i][j];
				meanIm[i](s) += Vin[im_index](s);						//averaging
			}
			meanIm[i](s) /= (j*b0(s));
			if ((meanIm[i](s))>1)
			{
				art=true;												//check for artefact that is ratio>1
			} 
		}
		else
		{meanIm[i](s)=0;}
		}
	if (art)
	{
		for (unsigned int i = 0; i < no_b; i += 1)
		{
			meanIm[i](s)=0;												//if artefact 0 the voxel
		}
	}
	prog++;
	printf("%.2f%% image averaged\r",(static_cast<float>(100*prog)/nel));
	}
	return(meanIm);
	
}


static dvecD getSortBsheels(const vfimage& Vin, std::vector< std::vector< unsigned int>>& sort_ind, fimage& b0)
//!looks up gradient strngths from dwi Image and sorts them in sort_b and their initial index in sort_ind
//structure of sort_ind is like cell_array with 1 cell containing all the indices (wrt Vin) of a given B-val
{
	std::vector< std::vector<float>> sort_b;
	unsigned int no_image = Vin.size();
	for (unsigned int i = 0; i < no_image; i += 1)						//looping images
	{	unsigned int cin = Vin[i].at.lookupNumber("gradient_strength");
		float Bin = static_cast<float>(0.001*(cin));
		if (i==0)
		{
			std::vector<float> b_arr; std::vector<unsigned int> i_arr;
			b_arr.push_back(Bin); i_arr.push_back(i);
			sort_b.push_back(b_arr); sort_ind.push_back(i_arr);
			continue;
		}
		bool ifpresent = false;
		unsigned int j = 0;
		for (; j < sort_b.size(); j += 1) 								//looping filled b-val
		{
			if ((Bin>=(0.9*sort_b[j][0])) && (Bin<=(1.1*sort_b[j][0])))  	//if b-val present insert to the vector
																		//10% range of b-val assumed
			{
				ifpresent = true;
				break;
			
			}
		}
				if(ifpresent)
				{
					sort_b[j].push_back(Bin);
					sort_ind[j].push_back(i);
				}
				else
				{														//else create new vector
					std::vector<float> b_arr; std::vector<unsigned int> i_arr;
					b_arr.push_back(Bin); i_arr.push_back(i);
					sort_b.push_back(b_arr); sort_ind.push_back(i_arr);
				}
	}
	unsigned int no_b = sort_b.size();
	if (no_b<3)
	{
		 throw optException("Multi-shell acquisition required");		// must have at least 3 including b0
	}
	dvecD b_sh(no_b-1);
	for (unsigned int i = 0; i < no_b-1; i += 1)
		{
			b_sh(i) = 0.0f;													//finding mean b-val	
		}
	unsigned int i;
	for (i = 0; i < no_b; i += 1)										//finding index of image with b-val=0
	{
		if (sort_b[i][0]==0)
		{
			break;
		}

	}
	fimage temp(Vin[0]); 
	b0=Vin[0]; b0=0.0f;
	unsigned int j = 0;
	for(; j < sort_ind[i].size(); j++)
	{
		int b0_ind = sort_ind[i][j];				
		temp = Vin[b0_ind];													//extracting b0
		for(unsigned int k = 0; k < temp.nel(); k++)
		{
			b0(k)=b0(k)+temp(k);
		}
	}
	for(unsigned int k = 0; k < b0.nel(); k++)
	{
		if(b0(k)==0)
		 {b0(k)=0;} 
		else {b0(k)=b0(k)/j;}
	}
	sort_ind.erase(sort_ind.begin()+i);
	sort_b.erase(sort_b.begin()+i);
	for (unsigned int i = 0; i < no_b-1; i += 1)
	{	j=0;
		for ( j = 0; j < sort_b[i].size(); j += 1)
		{
			b_sh[i] += static_cast<double>(sort_b[i][j]);				//getting mean B-shell
		}
		if ( j > 0 ) {
			b_sh[i] /= (j);
			printf("B-shell values: %f\n",1000*b_sh[i]);
		}
			
	}
	return b_sh;	
}

vfimage estimateSMT(const vfimage& Vin)
//!retunrs voxelwise lp and lt
{	
	vfimage Par,meanIm;
	fimage b0;
	std::vector< std::vector< unsigned int>> sort_ind;
	dvecD b_sh = getSortBsheels(Vin, sort_ind, b0);
	meanIm = normalizeAverageImage(Vin,sort_ind,b0);
	Par = getPar(meanIm, b_sh);
	return Par;
}

int main(int argc, char** argv)
{	int status = 0;
	try {	
		FILE *inf = nullptr, *outf = nullptr;
		vfimage vin, Par; 
		option opts[] = { };			
		cmdline c(argc, argv, 0, opts); c.parse(&inf, &outf);			//Parsing cmdline input
		readImages(vin,inf);
		printf("Image loaded\n");
		Par = estimateSMT(vin);
		saveImages(Par,outf);
		closeFile(inf);closeFile(outf);

	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
