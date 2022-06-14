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

void test(int num_atoms=10)
{
	auto start = std::chrono::steady_clock::now();

	std::vector<torch::Tensor> pos;

	for(int i=0; i<num_atoms; i++)
	{
		pos.push_back(torch::randn({3}, torch::requires_grad()));
		std::cout << " === " << i << " ===\n" << pos[i] << std::endl;
	}

	std::vector<float> RS = {1.0, 2.0, 3.0};
	std::vector<float> ETA = {0.5, 1.5};
	int num_features = RS.size()*ETA.size();

	torch::Tensor E = torch::zeros({1});

	auto net = std::make_shared<Net>(num_features, 1);

	for (int i=0; i<num_atoms; i++)
	{

		torch::Tensor G1 = torch::zeros({num_features});
		int idx=0;
		for (auto rs : RS)
		{
			for(auto eta : ETA)
			{
				for(int j=0; j<num_atoms; j++)
				{
					if (i==j) continue; 
					auto dr = torch::norm(pos[j]-pos[i]);
					G1[idx] = G1[idx] + torch::exp(-eta*(rs - dr)*(rs - dr));
				}
				idx += 1;
			}
		}

		auto energy = net->forward(G1);
		energy.backward();
		std::cout << "i,energy : " << i << " " << energy << "\n" << std::endl;
	}


	// Phase 2 : Create Network
	//net->to(torch::kCUDA);

	torch::optim::SGD optimizer(net->parameters(), 1e-6);

	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << std::endl;

}

extern "C" void init_rxmdnn(int natoms)
{
	std::cout << "foo from init\n";
	test(natoms);
}

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
	//maxrc = rxmdnn_ptr->get_maxrc_rxmdnn(); 
	//std::cout << "foo from maxrc " << maxrc << std::endl;
}
