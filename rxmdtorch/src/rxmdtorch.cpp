#include <iostream>
#include <chrono>
#include <memory>

#include <tuple>
#include <future>
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

struct model_spec
{
	std::string modelpath;
	int n_spieces;
	std::vector<std::string> element;
	std::vector<double> mass;
};

struct RXMDNN
{
	int myrank, nlocal, ntotal; 
	int deviceidx, devicecount; 
	torch::Device device = torch::kCPU;
	std::vector<torch::Device> devices; 

	double cutoff; 

	std::vector<torch::jit::script::Module> models;
	torch::jit::script::Module model;
	int current_model_id = 0; 

	std::vector<model_spec> model_specs;

	RXMDNN(int _myrank, std::vector<model_spec> _model_specs) : model_specs(_model_specs)
	{
		myrank = _myrank; 

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
			for (int deviceidx = 0; deviceidx < torch::cuda::device_count(); deviceidx++)
				devices.push_back(c10::Device(torch::kCUDA,deviceidx));

  		} else if(torch::xpu::is_available()) {

			device = c10::DeviceType::XPU;

		} else{

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

		for (const auto & mspec : model_specs)
		{
			for (auto & dev : devices)
			{
				std::cout << "loading model.. " << mspec.modelpath << std::endl;
				models.push_back(torch::jit::load(mspec.modelpath, dev, metadata));
			}
		}
		
		for (auto & m : models) m.eval();
		for (auto & m : models) torch::jit::optimize_for_inference(m);


		current_model_id = 0; 
		update_current_model(current_model_id);

		cutoff = std::stod(metadata["r_max"]);

		if(myrank == 0) 
		{
			std::cout << "Allegro is using device " << device << "\n";
			for (const auto & mspec : model_specs)
				std::cout << "modelpath: " << mspec.modelpath << std::endl;
			std::cout << "n_species: " << metadata["n_species"] << std::endl;
			std::cout << "r_max: " << metadata["r_max"] << std::endl;
			std::cout << "model.size(): " << model_specs.size() << std::endl;

			std::cout << "PyTorch version from parts: "
				<< TORCH_VERSION_MAJOR << "."
				<< TORCH_VERSION_MINOR << "."
				<< TORCH_VERSION_PATCH << std::endl;
				std::cout << "PyTorch version: " << TORCH_VERSION << std::endl;

		}

	}

	int get_num_models()
	{
		return model_specs.size();
	}

	void update_current_model(int const id)
	{
		if (id >= models.size()) 
			std::cout << "ERROR: model id> model.size():  " 
				<< id << " " << models.size() << std::endl;

		model = models[id];
	}

	double get_maxrc()
	{
		return cutoff; 
	}

	void get_nn_force(int const nlocal, int const ntotal, int const nbuffer, 
			void *pos_voidptr, void *type_voidptr, void *force_voidptr, void *nbrlist_voidptr, 
			double &energy)
	{
		//std::cout << "myrank,nlocal,ntotal,nbuffer: " << 
		//	myrank << " " << nlocal << " " << ntotal << " " << nbuffer<< std::endl;

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
		std::vector<int> cumsum_neigh_per_atom(nlocal, 0);
		for(int ii = 1; ii < nlocal; ii++)
			cumsum_neigh_per_atom[ii] = cumsum_neigh_per_atom[ii-1] + neigh_per_atom[ii-1];

		std::vector<torch::Tensor> pos_tensors; 
		std::vector<torch::Tensor> ij2type_tensors;

		for (auto & dev : devices)
		{
			torch::Tensor pos_tensor = torch::zeros({ntotal, 3});
			torch::Tensor ij2type_tensor = torch::zeros({ntotal}, torch::TensorOptions().dtype(torch::kInt64));

			auto pos = pos_tensor.accessor<float, 2>();
			auto ij2type = ij2type_tensor.accessor<long, 1>();

		// set position and type for all atoms including buffer atoms, 
		// but set complete edge information only for resident atoms. 
			for(int i = 0; i < ntotal; i++)
			{
				ij2type[i] = (int)type_vec[i] - 1;

				pos[i][0] = pos_[3*i];
				pos[i][1] = pos_[3*i+1];
				pos[i][2] = pos_[3*i+2];
    			}

			pos_tensors.push_back(pos_tensor.to(dev));
			ij2type_tensors.push_back(ij2type_tensor.to(dev));
		}


		int batch_size = std::ceil(nlocal/devices.size()); 
		int n_batches = nlocal/batch_size; 
		std::cout << "batch_size,n_batches: " << batch_size << " " << n_batches << std::endl;

		typedef torch::Dict<c10::IValue,c10::IValue> output_t;
		std::vector<std::future<std::tuple<output_t, int, int>>> outputs;

		//bool model_parallel = false; 
		bool model_parallel = true; 

		for (int nb = 0; nb <= n_batches; nb++)
		{
			int n_start = nb*batch_size;
			int n_end = std::min((nb+1)*batch_size, nlocal); 

			//std::cout << "n_start,n_end: " << n_start << " " << n_end << std::endl;

			if (n_start == n_end) break;

			outputs.push_back(std::async(
				model_parallel ? std::launch::async : std::launch::deferred, [&, nb, n_start, n_end] {

				//std::cout << "nb,n_start,n_end,nlocal: " << nb << " " 
				//<< n_start << " " << n_end << " " << nlocal << std::endl;

				int nedges = 0;
				for(int i = n_start; i < n_end; i++) nedges += neigh_per_atom[i];

				torch::Tensor edges_tensor = torch::zeros({2,nedges}, torch::TensorOptions().dtype(torch::kInt64));
				auto edges = edges_tensor.accessor<long, 2>();

				int iedges = 0; 
				for(int i = n_start; i < n_end; i++)
				{
					for(int j1 = 0; j1 < (int)nbrlist[i].size(); j1++)
					{
						int j = nbrlist[i][j1];
						edges[0][iedges] = i;
						edges[1][iedges] = j;
						iedges++;
					}
				}

				//std::cout << "myrank,nb/n_batches,iedges : " << myrank << 
				//	" " << nb << " / " << n_batches << " : " << iedges << std::endl;

				c10::Dict<std::string, torch::Tensor> input;

				input.insert("pos", pos_tensors[nb]);
				input.insert("atom_types", ij2type_tensors[nb]);

				auto dev = c10::Device(torch::kCUDA,nb);
				input.insert("edge_index", edges_tensor.to(devices[nb]));
				std::vector<torch::IValue> input_vector(1, input);

				//std::cout << "before calling m.forward.. " << std::endl;
				//auto output = model.forward(input_vector).toGenericDict();
				return std::make_tuple(models[nb].forward(input_vector).toGenericDict(), n_start, n_end);
				}
			));
		}

		double eng_vdwl = 0.0;
		for (auto &out_future: outputs) 
		{
			try { 
				//out_future.wait();
				auto out_tuple = out_future.get(); 

			       	output_t output = std::get<0>(out_tuple);
				int n_start = std::get<1>(out_tuple);
				int n_end = std::get<2>(out_tuple);

				//std::cout << "n_start,n_end: " << n_start << " " << n_end << std::endl;

				torch::Tensor atomic_energy_tensor = output.at("atomic_energy").toTensor().cpu();
				auto atomic_energies = atomic_energy_tensor.accessor<float, 2>();
				float atomic_energy_sum = atomic_energy_tensor.sum().data_ptr<float>()[0];

				torch::Tensor forces_tensor = output.at("forces").toTensor().cpu();
				auto forces = forces_tensor.accessor<float, 2>();

				for(int i = n_start; i < n_end; i++) eng_vdwl += atomic_energies[i][0];

				for (int ii=0; ii<ntotal; ii++)
				{
					force_vec[ii] += forces[ii][0];
					force_vec[nbuffer+ii] += forces[ii][1];
					force_vec[2*nbuffer+ii] += forces[ii][2];
				}

			} catch (const std::exception &e) {
				std::cout << e.what() << std::endl;
				exit(1);
			}

		}

		//std::cout << "total energy eng_vdwl: " << eng_vdwl << std::endl;
		energy = eng_vdwl;
	}
};

std::unique_ptr<RXMDNN> rxmdnn_ptr; 


extern "C" void init_rxmdtorch(int myrank)
{
	std::vector<model_spec> model_specs; 

        std::ifstream in("rxmdnn.in");
        std::string line, word, modelpath;
        while (std::getline(in,line))
        {
                if (line[0]=='#') continue;

                std::stringstream ss(line);

		model_spec mspec; 

                int n_spieces;

                ss >> word; // model name
                ss >> modelpath;
                ss >> n_spieces;

		mspec.modelpath = modelpath;
		mspec.n_spieces = n_spieces;

                std::string element;
                double mass;

                for(int i=0; i<n_spieces; i++)
                {
                        ss >> element; ss >> mass;

			mspec.element.push_back(element);
			mspec.mass.push_back(mass);
                }

		if(myrank == 0)
		{
                	std::cout << "modelpath: " << mspec.modelpath << std::endl;
                	std::cout << "n_spieces: " << mspec.n_spieces<< std::endl;
                	for(int i=0; i<n_spieces; i++) std::cout << mspec.element[i] << " " << mspec.mass[i] << ", ";
                	std::cout << std::endl;
		}

		model_specs.push_back(mspec);
        }

	//std::cout << "init_rxmdtorch : myrank " << myrank << std::endl;
	rxmdnn_ptr = std::make_unique<RXMDNN>(myrank, model_specs);
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

extern "C" void get_num_models_rxmdnn(int & n_models)
{
	n_models = rxmdnn_ptr->get_num_models();
	//std::cout << "get_num_models " << n_models << std::endl;
}

extern "C" void update_current_model_rxmdnn(int const id)
{
	rxmdnn_ptr->update_current_model(id);
	//std::cout << "update_current_model " << id << std::endl;
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = rxmdnn_ptr->get_maxrc();
	//std::cout << "get_maxrc_rxmdnn " << maxrc << std::endl;
}
