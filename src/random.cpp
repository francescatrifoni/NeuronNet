#include "random.h"


RandomNumbers::RandomNumbers(unsigned long int s) : seed(s)
{
	if (seed == 0){
		 std::random_device rd;
		 seed = rd();	
	 }
	rng = std::mt19937(seed);
}



void RandomNumbers::uniform_double(std::vector<double>& vect, double lower, double upper)
{
	std::uniform_real_distribution<> uniform_double(lower, upper);
    for (auto I = vect.begin(); I != vect.end(); I++) *I = uniform_double(rng);
}



double RandomNumbers::uniform_double(double lower, double upper)
{	
	std::uniform_real_distribution<>unid(lower,upper);
}  

void RandomNumbers::normal(std::vector<double>& vect, double mean, double sd)
{
	std::normal_distribution<> normal(mean, sd);
    for (auto I = vect.begin(); I != vect.end(); I++) *I = normal(rng);
}



double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<> normd(mean,sd);
}

void RandomNumbers::poisson(std::vector<int>& vect, double mean)
{
	vect.push_back(RandomNumbers::poisson(mean));
	
}



int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<> poisd(mean);
}


    
