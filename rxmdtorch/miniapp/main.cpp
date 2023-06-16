#include <iostream>
#include <vector>
#include "rxmdtorch.h"

int main(void)
{
	int nlocal=4, ntotal=4, nbuffer=100;
	double energy = -1.0;

	std::vector<double> pos0 = { 
		0.257, -0.363,  0.000, // N
		0.257,  0.727,  0.000, // H
		0.771, -0.727,  0.890, // H
		0.771, -0.727, -0.890, // H
	};
	std::vector<double> atype0={2.0, 1.0, 1.0, 1.0};

	std::vector<double> atype(nbuffer,0.0);
	std::vector<double> pos(3*nbuffer,0.0);
	std::vector<double> f(3*nbuffer,0.0);
	std::vector<signed int> nbr(nbuffer*10);
	nbr = {3, 2, 3, 4, 1, 1, 1, 1, 1, 1};

	for(int i=0; i<4; i++)
	{
		atype[i] = atype0[i];
		pos[i] = pos0[3*i];
		pos[i+nbuffer] = pos0[3*i+1];
		pos[i+2*nbuffer] = pos0[3*i+2];
	}

	std::cout << "in reproducer main()\n";

	init_rxmdtorch(0);
	get_nn_force_torch(nlocal, ntotal, nbuffer, pos.data(), atype.data(), f.data(), nbr.data(), energy);

	for(int i=0; i<4; i++)
		std::cout << "f: " << i << " " << f[i] << " " << f[i+nbuffer] << " " << f[i+2*nbuffer] << std::endl;
}
