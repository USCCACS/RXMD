/*
This code takes several .bnd files and calculates averaged bond-order values of neighbor atoms,
then print them only if their BO values are above a given threshould. 

g++ -std=c++11 BondLifeTime.cpp && ./a.out ../DAT/000*000.bnd 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <stdint.h>
#include <map>
#include <set>
#include <iomanip>

struct Atom
{
	float x, y, z;
	uint32_t id, type, nnbr=0; 
};

struct BondOrder : public Atom
{
	std::set <uint32_t> nbrIDs; 
	std::map <uint32_t, float> nbrBOs; 

	void Print()
	{
		uint intw=6,fltw=8; 
		std::cout.precision(3);
		std::cout << std::setw(intw) << id << 
				std::setw(fltw) << x << std::setw(fltw) << y << std::setw(fltw) << z << 
				std::setw(intw) << type << std::setw(intw) << nnbr << " "; 

		for(auto n : nbrIDs)
			std::cout << std::setw(intw) << n << std::setw(fltw) << nbrBOs[n];

		std::cout << std::endl;
	}
};

std::map<uint32_t, BondOrder> 
GetStableBonds(std::map<uint32_t, BondOrder> & bos, int const& nframes, float const& threshold)
{

    for(auto & b : bos)
    {
        auto & bo = b.second;
        for(auto n : bo.nbrIDs)
		{
            bo.nbrBOs[n]/=nframes;
            if(bo.nbrBOs[n] < threshold) 
			{
				bo.nbrIDs.erase(n);
				bo.nbrBOs.erase(n);
				std::cout << n << "will be removed: AveBo = " << bo.nbrBOs[n] << std::endl;
			}
		}

		// update # of neighbors
		bo.nnbr = bo.nbrIDs.size();
    }

	return bos;
};

int main(int argc, char *argv[])
{
	// use map to store BondOrder objects in case some atoms are missing.
	std::map<uint32_t, BondOrder> bos;

	// go over all bond files
	for(uint32_t i=1; i<argc; i++)
	{
		std::cout << argv[i] << std::endl; 
		std::ifstream input(argv[i]);
		
		std::string line; 
		while(getline(input,line))
		{
			std::stringstream ss(line); 

			// read basic atom info
			Atom a;
			ss >> a.id >> a.x >> a.y >> a.z >> a.type >> a.nnbr; 

			// get neighbor atom info
			auto & bo = bos[a.id];
			uint32_t nbrId; float nbrBo; 

			// read the rest of neighbor ids&bos
			while(ss)
			{
				ss >> nbrId >> nbrBo;
				//a.nbrs.push_back(std::make_pair(id,bo));
				bo.nbrIDs.insert(nbrId);
				bo.nbrBOs[nbrId]+=nbrBo; 
			}

			// keep atom info for stat. note: due to PBC, averaging the coords won't be correct. 
			bo.id = a.id; bo.type = a.type; // non-summations
			bo.x += a.x; bo.y += a.y; bo.z += a.z; bo.nnbr += a.nnbr; // summed up 
		}

	}

	// select neighbor atoms with averaged BO > 0.3
	for(auto & b : GetStableBonds(bos, argc, 0.3))
		b.second.Print();

	return 0;
}
