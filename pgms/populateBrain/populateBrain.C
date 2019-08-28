

#include <brian.h>
#include <dwi.h>
#include <math.h>
#include <vector>
#include <distributions.h>
#include <pointset.h>
#include <algorithm>

static std::vector<unsigned int> split(const std::string& s, char delimiter)
{
   std::vector<unsigned int> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
	  tokens.push_back(atoi(token.c_str()));
   }
   return tokens;
}

static std::vector< uvecD > read_csv( FILE *ifp)
{	
	std::vector< uvecD > vec;
	uvecD sam(2);int b,c,a;
	while (fscanf(ifp, "%d %d\n", &a, &b) != EOF) {
  		sam(0)=static_cast<unsigned int>(a), sam(1)=static_cast<unsigned int>(b);
  		//printf("%d %d\n",sam(0),sam(1));
  		vec.push_back(sam);
	}
	return vec;
}


static limage populateImagespace(const std::vector< uvecD>& clusvec, std::vector< unsigned int> clus_no)
{
	//limage clusimage(144,168,111);
	limage clusimage(128,128,64);
	for (unsigned int i = 0; i < clusvec.size(); i += 1)
	{
	//	if(clusvec[i](1)==clus_no)
	if(std::binary_search(clus_no.begin(), clus_no.end(), clusvec[i](1)))
	//check if the cluster assignment no present in clus_no
		{clusimage(clusvec[i](0)) = clusvec[i](1);}
	}
	return clusimage;
}

int main(int argc, char** argv)
//!input: .\populateBrain $1_rclus.csv $1_grp_16 1,2,3
{	int status = 0;
	try {	
		FILE *inf = fopen(argv[1],"r"),
		*outf = fopen(argv[2],"w");
		std::vector< unsigned int> clus_no = split(argv[3],',');
		std::sort(clus_no.begin(),clus_no.end());
		fimage pop_brain;
		std::vector< uvecD> value_vec = read_csv(inf);
		printf("file loaded\n");
		limage clusimage = populateImagespace( value_vec, clus_no);
		clusimage.save(outf);
		closeFile(inf);closeFile(outf);
	}
	
	catch (...) { brianException(); status = 2; }
	return status;
}		
