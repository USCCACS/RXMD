#include <iostream>
#include <chrono>
#include <memory>

#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>
#include <numeric>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <torch/torch.h>
#include <torch/script.h>
#include <torch/csrc/jit/runtime/graph_executor.h>

struct RXMDNN
{
	int myrank, nlocal, ntotal; 
	int deviceidx, devicecount; 
	torch::Device device = torch::kCPU;

	double cutoff; 

	torch::jit::script::Module model;

	RXMDNN(int myrank, std::string modelpath)
	{
		if(torch::cuda::is_available())
		{

			deviceidx = myrank;

			if(deviceidx >= 0) 
			{
				devicecount = torch::cuda::device_count();

				if(deviceidx >= devicecount) 
					deviceidx = deviceidx % devicecount;
			}

			device = c10::Device(torch::kCUDA,deviceidx);

  		} else {
			device = torch::kCPU;
		};

		std::unordered_map<std::string, std::string> metadata = {
			{"config", ""},
			{"nequip_version", ""},
			{"r_max", ""},
			{"n_species", ""},
			{"type_names", ""},
			{"_jit_bailout_depth", ""},
			{"_jit_fusion_strategy", ""},
			{"allow_tf32", ""}
		};

		model = torch::jit::load(modelpath, device, metadata);
		model.eval();

		cutoff = std::stod(metadata["r_max"]);

		if(myrank == 0) 
		{
			std::cout << "Allegro is using device " << device << "\n";
			std::cout << "modelpath: " << modelpath << std::endl;
			std::cout << "n_species: " << metadata["n_species"] << std::endl;
			std::cout << "r_max: " << metadata["r_max"] << std::endl;
		}
	}

	double get_maxrc()
	{
		return cutoff; 
	}

	void get_nn_force(int const nlocal, int const ntotal, int const nbuffer, 
			void *pos_voidptr, void *type_voidptr, void *force_voidptr, void *nbrlist_voidptr, double &energy)
	{
		//std::cout << "nlocal,ntotal,nbuffer" << nlocal << " " << ntotal << " " << nbuffer<< std::endl;
		double *pos_vec0 = (double *) pos_voidptr;
		double *type_vec = (double *) type_voidptr;
		double *force_vec = (double *) force_voidptr;
		signed int *nbrlist_vec = (signed int *) nbrlist_voidptr;

		std::vector<float> pos_(3*ntotal);
		for(int i=0; i<ntotal; i++)
		{
			pos_[3*i] = (float)pos_vec0[i];
			pos_[3*i+1] = (float)pos_vec0[nbuffer+i];
			pos_[3*i+2] = (float)pos_vec0[2*nbuffer+i];
		}

		// Total number of bonds (sum of number of neighbors)
		int nedges = 0;

		// Number of bonds per atom
		std::vector<int> neigh_per_atom(nlocal, 0);
		std::vector<std::vector<int>> nbrlist(ntotal);

		int c = 0; 
		for(int i=0; i<nlocal; i++)
                {
			signed int ni = nbrlist_vec[c++];
			for(int j1=0; j1<ni; j1++)
			{
				signed int j = nbrlist_vec[c++]-1; // zero-indexed
				nbrlist[i].push_back(j);
				neigh_per_atom[i]++;
				nedges++;
			}
		}

		// Cumulative sum of neighbors, for knowing where to fill in the edges tensor
		std::vector<int> cumsum_neigh_per_atom(nlocal);

		for(int ii = 1; ii < nlocal; ii++)
			cumsum_neigh_per_atom[ii] = cumsum_neigh_per_atom[ii-1] + neigh_per_atom[ii-1];

		torch::Tensor pos_tensor = torch::zeros({ntotal, 3});
		torch::Tensor edges_tensor = torch::zeros({2,nedges}, torch::TensorOptions().dtype(torch::kInt64));
		torch::Tensor ij2type_tensor = torch::zeros({ntotal}, torch::TensorOptions().dtype(torch::kInt64));

		auto pos = pos_tensor.accessor<float, 2>();
		auto edges = edges_tensor.accessor<long, 2>();
		auto ij2type = ij2type_tensor.accessor<long, 1>();

		// set position and type for all atoms including buffer atoms, 
		// but set complete edge information only for resident atoms. 
		for(int i = 0; i < ntotal; i++)
		{
			//ij2type[i] = type_mapper[itype - 1];
			ij2type[i] = (int)type_vec[i] - 1;

			pos[i][0] = pos_[3*i];
			pos[i][1] = pos_[3*i+1];
			pos[i][2] = pos_[3*i+2];

			//std::cout << ij2type[i] << " " << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << std::endl;

			if(i >= nlocal) continue;

			int edge_counter = cumsum_neigh_per_atom[i];
			for(int j1 = 0; j1 < (int)nbrlist[i].size(); j1++)
			{
				int j = nbrlist[i][j1];

				// TODO: double check order
				edges[0][edge_counter] = i;
				edges[1][edge_counter] = j;

				edge_counter++;
			}
    		}

		c10::Dict<std::string, torch::Tensor> input;
		input.insert("pos", pos_tensor.to(device));
		input.insert("edge_index", edges_tensor.to(device));
		input.insert("atom_types", ij2type_tensor.to(device));
		std::vector<torch::IValue> input_vector(1, input);

		auto output = model.forward(input_vector).toGenericDict();

		torch::Tensor forces_tensor = output.at("forces").toTensor().cpu();
		auto forces = forces_tensor.accessor<float, 2>();

		//torch::Tensor total_energy_tensor = output.at("total_energy").toTensor().cpu(); WRONG WITH MPI

		torch::Tensor atomic_energy_tensor = output.at("atomic_energy").toTensor().cpu();
		auto atomic_energies = atomic_energy_tensor.accessor<float, 2>();
		float atomic_energy_sum = atomic_energy_tensor.sum().data_ptr<float>()[0];

		/*
		std::cout << "atomic energy sum: " << atomic_energy_sum << std::endl;
		//std::cout << "Total energy: " << total_energy_tensor << "\n";
		std::cout << "atomic energy shape: " << atomic_energy_tensor.sizes()[0] << "," << atomic_energy_tensor.sizes()[1] << std::endl;
		std::cout << "atomic energies: " << atomic_energy_tensor << std::endl;
		*/

		// Write forces and per-atom energies (0-based tags here)
		double eng_vdwl = 0.0;

		for(int i = 0; i < ntotal; i++)
		{
			force_vec[i] = forces[i][0];
			force_vec[nbuffer+i] = forces[i][1];
			force_vec[2*nbuffer+i] = forces[i][2];
			//if (eflag_atom && i < nlocal) eatom[i] = atomic_energies[i][0];
			if(i < nlocal) eng_vdwl += atomic_energies[i][0];
  		}

		//std::cout << "total energy eng_vdwl: " << eng_vdwl << std::endl;
		energy = eng_vdwl;
	}
};

std::unique_ptr<RXMDNN> rxmdnn_ptr; 

extern "C" void init_rxmdtorch(int myrank)
{
        std::ifstream in("rxmdnn.in");
        std::string line, word, modelpath;
        while (std::getline(in,line))
        {
                if (line[0]=='#') continue;

                std::stringstream ss(line);

                int n_spieces=0;

                ss >> word; // model name
                ss >> modelpath;
                ss >> n_spieces;

                std::vector<std::string> element(n_spieces);
                std::vector<double> mass(n_spieces);

                for(int i=0; i<n_spieces; i++)
                {
                        ss >> element[i];
                        ss >> mass[i];
                }

		if(myrank == 0)
		{
                	std::cout << "modelpath: " << modelpath << std::endl;
                	std::cout << "n_spieces: " << n_spieces<< std::endl;
                	for(int i=0; i<n_spieces; i++) std::cout << element[i] << " " << mass[i] << ", ";
                	std::cout << std::endl;
		}
        }

	//std::cout << "init_rxmdtorch : myrank " << myrank << std::endl;
	rxmdnn_ptr = std::make_unique<RXMDNN>(myrank, modelpath);
}

extern "C" void get_nn_force_torch(int nlocal, int ntotal, int nbuffer, 
		void *pos_ptr, void *type_ptr, void *force_ptr, void *nbrlist_ptr, double &energy)
{
	//std::cout << "nlocal,nbuffer " << nlocal << " " << nbuffer << std::endl;
	const auto start = std::chrono::steady_clock::now();
	rxmdnn_ptr->get_nn_force(nlocal, ntotal, nbuffer, pos_ptr, type_ptr, force_ptr, nbrlist_ptr, energy);
	const auto end = std::chrono::steady_clock::now();
	//std::cout << "time(s) " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()*1e-6 << std::endl;
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = rxmdnn_ptr->get_maxrc();
	//std::cout << "get_maxrc_rxmdnn " << maxrc << std::endl;
}
