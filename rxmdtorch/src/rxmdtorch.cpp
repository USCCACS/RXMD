#include <iostream>
#include <chrono>
#include <memory>

#include "torch/torch.h"
#include "torch/script.h"

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
#define DEVICE torch::kCUDA
//#define DEVICE torch::kCPU


struct RXMDTORCH
{
	std::vector<float> rs, eta;
	float rc;

	torch::jit::script::Module net;

	RXMDTORCH(std::string filename)
	{
		torch::set_num_threads(1);

		std::cout << "rxmdtorch init" << std::endl;

		try 
		{
			net = torch::jit::load(filename);
		}
			catch (const c10::Error& e) {
			std::cerr << "error loading the model\n";
		}

		std::cout << filename << " loaded" << std::endl;

		rs = {1.0, 2.0, 3.0};
		eta = {0.5, 1.5};
		rc = MAXRC;

		net.to(DEVICE);

		if (DEVICE == torch::kCPU) torch::set_num_threads(32);
	}

	void featurize_nbrlist(int const natoms, int const maxnbrs, void *nbrtype_voidptr, void *nbrdist_voidptr)
	{
		std::cout << "entered featurize_nbrlist\n";
		float *nbrdist_0 = (float *) nbrdist_voidptr; 
		int *nbrtype = (int *) nbrtype_voidptr; 

		int size = 4*natoms*maxnbrs;

		auto nbrdist = torch::empty({size}, torch::kFloat).to(DEVICE);

		for(int n=0; n<natoms; n++)
			for(int j=0; j<maxnbrs; j++)
				for(int k=4*(n*maxnbrs+j); k<4*(n*maxnbrs+j)+4; k++)
				{
					//std::cout << n << " " << j << " " << k << " " << n*maxnbrs + j << " " << 
					//	nbrtype[n*maxnbrs+j] << " " << nbrdist.sizes() << std::endl;
					nbrdist[k] = nbrdist_0[k];
				}
					
		int feature_size = eta.size()*rs.size();
	
		//auto net = std::make_shared<Net>(feature_size, 1);
	
		torch::Tensor total_energy = torch::zeros({1}, torch::kFloat).to(DEVICE);
		for (int n=0; n<natoms; n++)
		{
			auto G2 = torch::zeros({3*feature_size}, torch::kFloat).to(DEVICE);
	
			for (int j=0; j<maxnbrs; j++)
			{

				int ii4 = 4*(n*maxnbrs + j);
				auto dr = nbrdist[ii4];
				auto dx = nbrdist[ii4+1];
				auto dy = nbrdist[ii4+2];
				auto dz = nbrdist[ii4+3];

				auto jtype = nbrtype[n*maxnbrs+j];

				if(jtype == 0) continue;
				//std::cout << "n,j,jtype " << n << " " << j << " " << jtype << std::endl;
	
				if(dr.item<float>() > rc) continue;
				//if(dr > rc) continue;
	
				int idx = (jtype-1)*feature_size;
				for (unsigned long i1 = 0; i1 < rs.size(); i1++)
				{
					for(unsigned long i2 = 0; i2 < eta.size(); i2++)
					{
						auto rij_rs = dr - rs[i1]; 
						auto exp_rij = exp(-eta[i2] * rij_rs * rij_rs);
						auto fc_rij = 0.5*cos(M_PI*dr/rc);
	
						G2[idx] += exp_rij*fc_rij;
					}
				}
			}
	
			//auto energy = net->forward(G2);
			std::vector<torch::jit::IValue> inputs;
			inputs.push_back(G2);

			//std::cout << "G2 " << n << " " << G2.sizes() << "\n" << G2 << std::endl;
			auto energy = net.forward(inputs).toTensor();
			//std::cout << "energy " << n << " " << energy << std::endl;

			total_energy += energy;
		}
		total_energy.backward();
	}
};

std::unique_ptr<RXMDTORCH> rxmdnn_ptr; 

extern "C" void init_rxmdtorch(int natoms)
{
	rxmdnn_ptr = std::make_unique<RXMDTORCH>("nn.pt");
}

extern "C" void predict_rxmdtorch(int natoms, int maxnbrs, void *nbrtype_ptr, void *nbrdist_ptr)
{
	std::cout << "natoms,maxnbrs " << natoms << " " << maxnbrs << std::endl;

	const auto start = std::chrono::steady_clock::now();
	rxmdnn_ptr->featurize_nbrlist(natoms, maxnbrs, nbrtype_ptr, nbrdist_ptr);
	const auto end = std::chrono::steady_clock::now();
	std::cout << "time(s) " << 
		std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()*1e-6 << std::endl;
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = MAXRC;
	std::cout << "get_maxrc_rxmdnn " << maxrc << std::endl;
}
