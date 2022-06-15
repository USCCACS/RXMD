#include <iostream>
#include <chrono>
#include <memory>
#include "torch/torch.h"

struct Net : torch::nn::Module {
	Net(int in_dim, int out_dim) {
		fc1 = register_module("fc1", torch::nn::Linear(in_dim, 20));
		fc2 = register_module("fc2", torch::nn::Linear(20, 20));
		fc3 = register_module("fc3", torch::nn::Linear(20, out_dim));
	}

	torch::Tensor forward(torch::Tensor x) {
		x = fc1->forward(x);
		x = fc2->forward(x);
		x = fc3->forward(x);
		return x;
	}

	torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr};
};

#define MAXRC 5.5

void rxmdtorch_featurize_nbrlist(int const natoms, int const maxnbrs, void *nbrdist_voidptr)
{
	float *nbrdist_0 = (float *) nbrdist_voidptr; 
	torch::Tensor nbrdist = torch::zeros({4*natoms*maxnbrs});
	for(int n=0; n<natoms; n++)
		for(int j=0; j<maxnbrs; j++)
			for(int k=4*n*j; k<4*n*j+4; k++)
				nbrdist[k] = nbrdist_0[k];

	std::vector<float> rs = {1.0, 2.0, 3.0};
	std::vector<float> eta = {0.5, 1.5};
	float rc = MAXRC;

	int feature_size = eta.size()*rs.size();

	int eta_size = eta.size();
	int rs_size = rs.size();

	auto net = std::make_shared<Net>(feature_size, 1);

	torch::Tensor total_energy = torch::zeros({1});
	for (int n=0; n<natoms; n++)
	{
		torch::Tensor G2 = torch::zeros({feature_size});

		for (int j=0; j<maxnbrs; j++)
		{
			if (n==j) continue;

			int ii4 = 4*(n*maxnbrs + j);
			auto dr = nbrdist[ii4];
			auto dx = nbrdist[ii4+1]/dr;
			auto dy = nbrdist[ii4+2]/dr;
			auto dz = nbrdist[ii4+3]/dr;

			auto greater = dr > rc;
			if(greater.item<bool>()) continue;

			int idx = 0;
			for (int i1 = 0; i1 < rs.size(); i1++)
			{
				for(int i2 = 0; i2 < eta.size(); i2++)
				{
					auto rij_rs = dr - rs[i1]; 
					auto exp_rij = exp(-eta[i2] * rij_rs * rij_rs);
					auto fc_rij = 0.5*cos(M_PI*dr/rc);

					G2[idx] += exp_rij*fc_rij;
				}
			}
		}

		auto energy = net->forward(G2);
		total_energy += energy;
	}

	total_energy.backward();
}


//==========================================================================
extern "C" void init_rxmdtorch(int natoms)
{
	std::cout << "foo from init\n";
	//test(natoms);
}

extern "C" void predict_rxmdtorch(int natoms, int atom_type, int maxnbrs, void *nbrdist_ptr)
{
	//std::cout << "foo from predict_hybrid\n";
	//predict(natoms, atom_type, maxnbrs, nbrdist_ptr);

	std::cout << "\n\nEntering rxmdtorch_predict()\n\nnatoms,atom_type,maxnbrs: " << 
		natoms << " " << atom_type << " " << maxnbrs << std::endl;

	rxmdtorch_featurize_nbrlist(natoms, maxnbrs, nbrdist_ptr);
}
//==========================================================================


extern "C" void init_rxmdnn_hybrid(int natoms)
{
	std::cout << "foo from init\n";
	//rxmdnn_ptr = std::make_unique<RXMDNN>(natoms,"rxmdnn.in");
}

extern "C" void predict_rxmdnn_hybrid(int natoms, int atom_type, int maxnbrs, void *nbrdist_ptr)
{
	std::cout << "foo from predict_hybrid\n";
	//rxmdnn_ptr->predict(natoms, atom_type, maxnbrs, nbrdist_ptr);
}

/*
extern "C" void init_rxmdnn(void)
{
	std::cout << "foo from init\n";
	//rxmdnn_ptr = std::make_unique<RXMDNN>("pto.xyz","rxmdnn.in");
}
*/

extern "C" void predict_rxmdnn(void)
{
	std::cout << "foo from predict\n";
	int atom_type = 1;
	//rxmdnn_ptr->predict(rxmdnn_ptr->md.natoms, atom_type, rxmdnn_ptr->md.natoms, (void *) rxmdnn_ptr->nbr.nbrdist.data()); 
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = MAXRC;
	//maxrc = rxmdnn_ptr->get_maxrc_rxmdnn(); 
	//std::cout << "foo from maxrc " << maxrc << std::endl;
}
