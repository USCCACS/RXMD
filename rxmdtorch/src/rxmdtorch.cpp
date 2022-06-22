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
//#define DEVICE torch::kCUDA
#define DEVICE torch::kCPU


struct RXMDTORCH
{
	std::vector<float> RS, ETA;
	float RC;
	int feature_size;

	std::vector<torch::jit::script::Module> nets;

	RXMDTORCH(std::string filename)
	{
		torch::set_num_threads(1);

		std::cout << "rxmdtorch init" << std::endl;

		try 
		{
			nets.push_back(torch::jit::load("net_Pb.pt"));
			nets.push_back(torch::jit::load("net_Ti.pt"));
			nets.push_back(torch::jit::load("net_O.pt"));
		}
			catch (const c10::Error& e) {
			std::cerr << "error loading the model\n";
		}

		std::cout << filename << " loaded" << std::endl;

		RS = {1.0, 2.0, 3.0, 4.0, 5.0};
		ETA = {1.0, 1.5, 2.0};
		RC = MAXRC;

		feature_size = RS.size()*ETA.size();

		//net.to(DEVICE);

		if (DEVICE == torch::kCPU) torch::set_num_threads(32);
	}

	void get_nn_force(int const natoms, int const nbuffer, int const maxnbrs, 
			void *pos_voidptr, void *type_voidptr, void *force_voidptr)
	{
		std::cout << "natoms,nbuffer,maxnbrs " << natoms << " " << nbuffer << " " << maxnbrs << std::endl;
		double *pos_vec0 = (double *) pos_voidptr; 
		double *type_vec = (double *) type_voidptr; 
		double *force_vec = (double *) force_voidptr; 

		//FIXME copy & align all position data 
		std::vector<float> pos_vec(3*nbuffer);
		for(int i=0; i<nbuffer; i++)
		{
			pos_vec[3*i] = (float)pos_vec0[i];
			pos_vec[3*i+1] = (float)pos_vec0[nbuffer+i];
			pos_vec[3*i+2] = (float)pos_vec0[2*nbuffer+i];
		}

		std::vector<torch::Tensor> pos(nbuffer); 
		for(int i=0; i<nbuffer; i++)
			pos[i] = torch::from_blob(&pos_vec[3*i], {3}, torch::requires_grad());

		/*
		for(int i=0; i<nbuffer; i++)
		{
			std::cout << "===== " << i << " " << type_vec[i] << " =======\n";
			std::cout << pos[i] << std::endl;
		}
		*/

		for(int i=0; i<natoms; i++)
		{
			int itype = type_vec[i];

			auto g2 = torch::zeros({3*feature_size});
			for(int j=0; j<nbuffer; j++)
			{
				if(i==j) continue;
				if(type_vec[j] == 0) continue;

				auto dr = torch::linalg_norm(pos[j] - pos[i]);
				//std::cout << " pos[i] ========= \n";
				//std::cout << pos[i] << std::endl;
				//std::cout << " pos[j] ========= \n";
				//std::cout << pos[j] << std::endl;
				//std::cout << i << " " << ia << " " << ib << " " << j << " " << dr << std::endl;
				if(dr.item<float>()>RC) continue;

				int stride = (type_vec[j]-1)*feature_size;

				int idx=0;
				for(int ia=0; ia < (int) RS.size(); ia++)
				{
					auto rij_rs = dr - RS[ia]; 
					auto fc_rij = 0.5*torch::cos(M_PI*dr/RC);

					for(int ib=0; ib < (int) ETA.size(); ib++)
					{	
						auto exp_rij = torch::exp(-ETA[ib] * rij_rs * rij_rs);
	
						g2[idx+stride] += exp_rij*fc_rij;
					}
					idx++;
				}
			}

			std::vector<torch::jit::IValue> inputs;
			inputs.push_back(g2);
			auto energy = nets[itype-1].forward(inputs).toTensor();

			energy.backward();

			//std::cout << "====== " << i << " " << itype << " ======" << std::endl;
			//std::cout << energy << std::endl;
			//std::cout << g2 << std::endl;
			//std::cout << pos[i].grad() << std::endl;

			auto force = pos[i].grad();
			force_vec[i] = force[0].item<float>();
			force_vec[i+nbuffer] = force[1].item<float>();
			force_vec[i+nbuffer*2] = force[2].item<float>();
		}
	}
};

std::unique_ptr<RXMDTORCH> rxmdnn_ptr; 

extern "C" void init_rxmdtorch()
{
	rxmdnn_ptr = std::make_unique<RXMDTORCH>("nn.pt");
}

extern "C" void get_nn_force_torch(int natoms, int nbuffer, int maxnbrs, 
		void *pos_ptr, void *type_ptr, void *force_ptr)
{
	std::cout << "natoms,maxnbrs " << natoms << " " << maxnbrs << std::endl;
	const auto start = std::chrono::steady_clock::now();
	rxmdnn_ptr->get_nn_force(natoms, nbuffer, maxnbrs, pos_ptr, type_ptr, force_ptr);
	const auto end = std::chrono::steady_clock::now();
	std::cout << "time(s) " << 
		std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()*1e-6 << std::endl;
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = MAXRC;
	std::cout << "get_maxrc_rxmdnn " << maxrc << std::endl;
}
