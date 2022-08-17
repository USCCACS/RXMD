#include <iostream>
#include <chrono>
#include <memory>

#include "torch/torch.h"
#include "torch/script.h"

using namespace torch::indexing;

#define MAXRC 7.5
//#define DEVICE torch::kCUDA
#define DEVICE torch::kCPU

struct RXMDTORCH
{
	std::vector<float> RS, ETA;
	float RC;
	int num_hparams, feature_size;
	std::map<std::string, int> type_numeric;

	std::vector<torch::jit::script::Module> nets;

	RXMDTORCH(std::string dirname = "")
	{
		torch::set_num_threads(1);

		std::cout << "rxmdtorch init" << std::endl;

		//std::vector<std::string> filenames = {"H.model_scripted.pt", "N.model_scripted.pt"};
		std::vector<std::string> filenames = {"Pb.model_scripted.pt", "Ti.model_scripted.pt", "O.model_scripted.pt"};

		for (std::string filename : filenames)
		{
			filename = dirname + filename;  
			try 
			{
				nets.push_back(torch::jit::load(filename));
			}
			catch (const c10::Error& e) 
			{
				std::cerr << "error loading the model : " + filename << std::endl;
			}
		}

		RS = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
		ETA = {0.5,1.0,3.0};
		RC = MAXRC;

		type_numeric = {{"Pb",1}, {"Ti",2}, {"O",3}};

		num_hparams = RS.size()*ETA.size();
		feature_size = 3*num_hparams; 

		for(auto net : nets) net.to(DEVICE);
		for(auto net : nets) net.eval();

		//if (DEVICE == torch::kCPU) torch::set_num_threads(32);
	}

	void get_nn_force(int const natoms, int const nbuffer, int const maxnbrs, 
			void *pos_voidptr, void *type_voidptr, void *force_voidptr)
	{
		std::cout << "natoms,nbuffer,maxnbrs " << natoms << " " << nbuffer << " " << maxnbrs << std::endl;
		double *pos_vec0 = (double *) pos_voidptr; 
		double *type_vec = (double *) type_voidptr; 
		double *force_vec = (double *) force_voidptr; 

		//FIXME copy & align all position data 
		//std::vector<float> pos_vec(3*nbuffer);
		std::vector<float> pos(3*nbuffer);
		for(int i=0; i<nbuffer; i++)
		{
			pos[3*i] = (float)pos_vec0[i];
			pos[3*i+1] = (float)pos_vec0[nbuffer+i];
			pos[3*i+2] = (float)pos_vec0[2*nbuffer+i];
		}

		float total_energy=0.0;

		int counter = 0;
		for(int i=0; i<natoms; i++)
		{
			//std::cout << " counter " << counter << std::endl;
			//std::cout << " ====== " << i << " =============\n";
			int itype = type_vec[i];

			std::vector<float> g2_vec(feature_size);

			auto feature_jacobian_self=torch::zeros({feature_size,3});
			auto feature_jacobian_neighbor=torch::zeros({nbuffer,feature_size,3});

auto start = std::chrono::steady_clock::now();

			for(int j=0; j<nbuffer; j++)
			{
				if(i==j) continue;
				if(type_vec[j] == 0) continue;
				//auto rij = pos[i] - pos[j];
				//auto dr = torch::linalg_norm(rij)


				const int jtype = type_vec[j];
				const int stride = (jtype-1)*num_hparams; 

				const float dx = pos[3*i] - pos[3*j];
				const float dy = pos[3*i+1] - pos[3*j+1];
				const float dz = pos[3*i+2] - pos[3*j+2];
				const float dr = std::sqrt(dx*dx + dy*dy + dz*dz);

				//std::cout << "dx,dy,dz " << dx << " " << dy << " " << dz << std::endl;
				auto rij_ = torch::tensor({dx,dy,dz}, torch::requires_grad());
				//std::cout << "rji " << rij_ << std::endl;
				auto dr_ = torch::linalg_norm(rij_);

				//if(dr.item<float>()>RC) continue;
				if(dr>RC) continue;

				auto fc_rij = 0.5*(std::cos(M_PI*dr/RC)+1.0);
				auto dfc_rij = -0.5*(M_PI/RC)*std::sin(M_PI*dr/RC);

				//std::cout << dr_.item<float>() << " " << jtype << "\n";
				//counter++; 

				//int idx=0;
				int idx = stride;
				for(int ia=0; ia < (int) ETA.size(); ia++)
				{
					for(int ib=0; ib < (int) RS.size(); ib++)
					{	
						/*
						auto rij = torch::tensor({dx,dy,dz}, torch::requires_grad());
						auto dr = torch::linalg_norm(rij);

						auto rij_rs = dr - RS[ia];
						auto fc_rij = 0.5*(torch::cos(M_PI*dr/RC)+1.0);
						auto exp_rij = torch::exp(-ETA[ib] * rij_rs * rij_rs);
						//std::cout << "==========================\n";

						auto feature = exp_rij*fc_rij;
						//std::cout << "feature " << feature << std::endl;

						feature.backward();
						//std::cout << "rij.grad() " << rij.grad() << std::endl;

						//auto feat_d=dr.grad()*rij/dr;
						auto feat_d=rij.grad();
						//std::cout << "feat_d " << feat_d << std::endl;

						for (int ia=0; ia<3; ia++)
						{
							feature_jacobian_self.index({idx,ia}) = 
								feature_jacobian_self.index({idx,ia}) + feat_d[ia];
							feature_jacobian_neighbor.index({j,idx,ia}) = 
								feature_jacobian_neighbor.index({j,idx,ia}) - feat_d[ia];
						}

						//g2[idx]+=feature;
						//std::cout << "feature " << feature << std::endl;
						g2_vec[idx] += feature.item<float>();
						*/

						auto rij_rs = dr - RS[ib];
						auto exp_rij = std::exp(-ETA[ia] * rij_rs * rij_rs);
						auto dexp_rij = -2.0*ETA[ia]*exp_rij*rij_rs;

						auto feature = exp_rij*fc_rij;
						g2_vec[idx] += feature;

						/*
						if(idx==0) 
						std::cout << "jtype,idx,eta,rs,feature,fc_rij,exp_rij: " << 
						jtype << " " << idx << " " << ETA[ia] << " " << RS[ib] << " " << 
						feature << " " << fc_rij << " " << exp_rij << std::endl;
						*/

						auto coeff = (exp_rij*dfc_rij + dexp_rij*fc_rij)/dr;

						auto feat_d = torch::tensor({coeff*dx,coeff*dy,coeff*dz}, torch::requires_grad());

						//for (int ia=0; ia<3; ia++)
						//{
						//	feature_jacobian_self.index({idx,ia}) = feature_jacobian_self.index({idx,ia}) + feat_d[ia];
						//	feature_jacobian_neighbor.index({idx,ia}) = feature_jacobian_neighbor.index({idx,ia}) - feat_d[ia];
						//}
						
						feature_jacobian_self.index({idx,Slice()}) = feature_jacobian_self.index({idx,Slice()}) + feat_d;
						feature_jacobian_neighbor.index({j,idx,Slice()}) = feature_jacobian_neighbor.index({j,idx,Slice()}) - feat_d;

						idx++;
					}
				}
			}
			//std::cout << "jacobian_self\n" << feature_jacobian_self << std::endl;

auto end = std::chrono::steady_clock::now();
std::chrono::duration<double> es1 = end-start;

auto start2 = std::chrono::steady_clock::now();

			auto g2 = torch::from_blob(g2_vec.data(),{1,feature_size}).requires_grad_();
			//std::cout << "g2 " << g2 << std::endl;

			std::vector<torch::jit::IValue> inputs;

			inputs.push_back(g2);
			//std::cout << "inputs: " << inputs << std::endl;

			auto atomic_energy = nets[itype-1].forward(inputs).toTensor();
			//std::cout << "atomic_energy: " << atomic_energy << std::endl;

			total_energy += atomic_energy.item<float>();
			//std::cout << "total_energy: " << total_energy << std::endl;

			auto grad_output = torch::ones_like(atomic_energy);
			//std::cout << "grad_output: " << grad_output << std::endl;

			auto g2_ = g2.index({0,Slice()});
			//std::cout << "g2_ " << g2_ << std::endl;

			//auto atomic_energy_ = atomic_energy.index({0,Slice()});
			//std::cout << "atomic_energy_: " << atomic_energy_ << std::endl;

			// torch::autograd::grad returns a std::vector where the first element is the gradient tensor
 		        auto gradient = torch::autograd::grad({atomic_energy}, {g2})[0];
			//std::cout << "gradient " << gradient << std::endl;
			//std::cout << "feature_jacobian_self " << feature_jacobian_self << std::endl;

			// i-atom force
			auto force_i_ = torch::matmul(gradient, feature_jacobian_self);
			auto force_i = force_i_.index({0,Slice()});
			//std::cout << "force_i: " << force_i << std::endl;
			//std::cout << "force_i_: " << force_i_ << std::endl;

			force_vec[i] += force_i[0].item<float>();
			force_vec[i+nbuffer] += force_i[1].item<float>();
			force_vec[i+nbuffer*2] += force_i[2].item<float>();

auto end2 = std::chrono::steady_clock::now();
std::chrono::duration<double> es2 = end2-start2;

auto start3 = std::chrono::steady_clock::now();
			// Reaction force from neighbor atoms
			for(int j=0; j<nbuffer; j++)
			{	
				//auto feature_jacobian_j=feature_jacobian_neighbor[j,:,:];
				auto feature_jacobian_j = feature_jacobian_neighbor.index({j,Slice(None,None)});
				auto force_j_ = torch::matmul(gradient, feature_jacobian_j);
				auto force_j = force_j_.index({0,Slice()});
				//std::cout << "force_j: " << force_j << std::endl;

				force_vec[j] += force_j[0].item<float>();
				force_vec[j+nbuffer] += force_j[1].item<float>();
				force_vec[j+nbuffer*2] += force_j[2].item<float>();
			}

auto end3 = std::chrono::steady_clock::now();
std::chrono::duration<double> es3 = end3-start3;

//std::cout << "elapsed time: " << es1.count() << " " << es2.count() << " " << es3.count()<< "s\n";
		}

	}
};

std::unique_ptr<RXMDTORCH> rxmdnn_ptr; 

extern "C" void init_rxmdtorch()
{
	rxmdnn_ptr = std::make_unique<RXMDTORCH>("");
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
