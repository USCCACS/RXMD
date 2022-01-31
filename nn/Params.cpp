#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "json.hpp"

std::stringstream getss(std::ifstream & fin)
{
	std::string line;
	std::getline(fin,line); 
	//std::cout << line << std::endl;
	std::stringstream ss(line);
	return ss; 
}

template <typename T>
void read_params(std::ifstream & fin, std::vector<T> & v)
{
	auto ss = getss(fin);
	for(int i = 0; i<v.size(); i++)
	{
		if (ss.eof()) ss = getss(fin);
		T data;

		ss >> data; v[i] = data;
	}
}

template<typename T>
void read_G2(std::ifstream & fin, const std::vector<int> &sf, std::vector<std::vector<T>> & v)
{
	assert(sf.size() == v.size());

	auto ss = getss(fin);
	for(int i=0; i<v.size(); i++)
	{
		assert(sf[i]==2); //support only G2 symm function

		if(ss.eof()) ss = getss(fin);
		for(int j=0; j<v[i].size(); j++) ss >> v[i][j];
		//for(int j=0; j<v[i].size(); j++) std::cout << i << " " << j << " " << v[i][j] << std::endl;
	}
}


struct AenetParams
{
	std::string element, filename; 
	float mass; 

	// load_Network
	int nlayers, nnodes_max, Wsize, nvalues;
	std::vector<int> nnodes, f_a, iw, iv;
	std::vector<float> W;
	bool init, evaluated, derivatives; 
	
	// load_Setup
	std::string description; // unused
	std::string atomtype, sftype;
	int nenv, nsf, nsfparam; 
	std::vector<std::string> envtypes;
	float Rc_min, Rc_max;

	std::vector<int> sf;
	std::vector<std::vector<float>> sfparam;
	std::vector<std::vector<int>> sfenv;
	int neval;
	std::vector<float> sfval_min, sfval_max, sfval_avg, sfval_cov;

	std::vector<float> D, value, deriv, work, work2, work3, work4;

	AenetParams(const std::string & parminfo)
	{
		std::stringstream ss(parminfo);

		std::string token; 
		ss >> token >> element >> mass >> filename; 
		std::cout << token << " " <<  element << " " << mass << " " << filename << std::endl; 

		std::ifstream fin = std::ifstream(filename);

		load_Network(fin);
		load_Setup(fin);

		fin.close();
	};

	void load_Network(std::ifstream & fin)
	{
		getss(fin) >> nlayers; 
		getss(fin) >> nnodes_max; 
		getss(fin) >> Wsize; 
		getss(fin) >> nvalues; 

		nnodes.resize(nlayers);
		f_a.resize(nlayers-1);
		iw.resize(nlayers);
		iv.resize(nlayers);
		W.resize(Wsize);

		read_params(fin, nnodes);
		read_params(fin, f_a);
		read_params(fin, iw);
		read_params(fin, iv);
		read_params(fin, W);
	}

	void load_Setup(std::ifstream &fin)
	{
		// getss(fin) >> description; // unsed
		getss(fin) >> atomtype;
		getss(fin) >> nenv; 

		envtypes.resize(nenv); 
		for (int i=0; i<envtypes.size(); i++) getss(fin) >> envtypes[i]; 

		getss(fin) >> Rc_min; 
		getss(fin) >> Rc_max; 
		getss(fin) >> sftype; 
		getss(fin) >> nsf; 
		getss(fin) >> nsfparam; 

		sf.resize(nsf);
		read_params(fin, sf);

		sfparam.resize(nsf);
		for(int i=0;i<sfparam.size();i++) sfparam[i].resize(nsfparam);
		read_G2(fin, sf, sfparam);

		const int NENV_MAX=2;
		sfenv.resize(nsf);
		for(int i=0;i<sfenv.size();i++) sfenv[i].resize(NENV_MAX);
		read_G2(fin, sf, sfenv);
		
		getss(fin) >> neval; 

		sfval_min.resize(nsf); read_params(fin, sfval_min);
		sfval_max.resize(nsf); read_params(fin, sfval_max);
		sfval_avg.resize(nsf); read_params(fin, sfval_avg);
		sfval_cov.resize(nsf); read_params(fin, sfval_cov);
	}

	std::vector<int> get_nodes(void)
	{
		return nnodes;
	};

	void show(void)
	{
		std::cout << "=============================================\n";

		std::cout << "nlayers : " << nlayers << std::endl;
		std::cout << "nnodes_max : " << nnodes_max << std::endl;
		std::cout << "Wsize : " << Wsize << std::endl;
		std::cout << "nvalues : " << nvalues<< std::endl << std::endl;

		std::cout << "nnodes : "; for (auto i : nnodes) { std::cout << i << " ";}; std::cout << std::endl;
		std::cout << "f_a : "; for (auto i : f_a) { std::cout << i << " ";}; std::cout << std::endl;
		std::cout << "iw : "; for (auto i : iw) { std::cout << i << " ";}; std::cout << std::endl;
		std::cout << "iv : "; for (auto i : iv) { std::cout << i << " ";}; std::cout << std::endl; 
		std::cout << "W : "; 
		for (int i = 0; i<W.size(); i++) 
		{ 
			if(i%16 == 0) std::cout << std::endl;
			std::cout << std::scientific << std::setprecision(3) << std::setw(11) << W[i];
		}; 
		std::cout << std::endl;

		std::cout << "=============================================\n";

		std::cout << "atomtype,nenv,Rc_min,Rc_max: " 
			<< atomtype << " " << nenv << " " << Rc_min << " " << Rc_max << std::endl;
		std::cout << "sftype,nsf,nsfparam: " << sftype << " " << nsf << " " << nsfparam << std::endl;
		std::cout << "envtypes : "; for (auto e : envtypes) std::cout << e << " "; std::cout << std::endl; 

		std::cout << "sf : "; 
		for (int i = 0; i<sf.size(); i++) std::cout << sf[i] << " "; 
		std::cout << std::endl;

		std::cout << "sfparam : ";
		for (int i = 0; i<sfparam.size(); i++) 
		{ 
			if(i%4==0) std::cout << std::endl;
			for(auto a : sfparam[i]) std::cout << a << " "; std::cout << ", "; 
		}; 
		std::cout << std::endl;

		std::cout << "sfenv: " << std::endl;
		for (int i = 0; i<sfenv.size(); i++) 
		{ 
			if(i%16==0) std::cout << std::endl;
			for(auto a : sfenv[i]) std::cout << a << " "; std::cout << "  ";
		}; 
		std::cout << std::endl;

		std::cout << "=============================================\n";
	};

};

struct JAXMDParams
{
	nlohmann::json json;

	std::string element, filename; 
	float mass; 

	int feature_size; 

	std::vector<int> nnodes;

	JAXMDParams(const std::string & parminfo)
	{
		std::stringstream ss(parminfo);

		std::string token; 
		ss >> token >> element >> mass >> filename; 
		std::cout << token << " " <<  element << " " << mass << " " << filename << std::endl; 

		std::ifstream fin = std::ifstream(filename);
		fin >> json;

		const auto j = json["sym_hyper_parameters"];

		// add feature size to nnodes first
		//feature_size = j["radial_etas"].size() + j["angular_etas"].size();
		
		feature_size = j["radial_etas"].size(); //FIXME only radial feature is supposed now. 

		nnodes.push_back(feature_size);

		for(auto const & d : json.items())
		{
			auto const & key = d.key();
			auto const & data = d.value();
			if (key.find("Layer") != std::string::npos) nnodes.push_back(data["Bias"].size());
		};


		fin.close();
	}

	void show(void)
	{
		std::cout << "Activation : " << json["Activation"] << std::endl;
		std::cout << "sym hyperparams : " << std::endl;

		const auto j = json["sym_hyper_parameters"];
		for (auto const & d : j.items())
		{
			std::cout << d.key() << " " << d.value().size() << std::endl;
		}

		std::cout << "radial features" << std::endl;
		for(int i=0; i<j["radial_etas"].size(); i++)
			std::cout << j["radial_etas"][i] << " " << 0.0 << " " << j["cutoff_distance"] << std::endl;

		std::cout << "angular features" << std::endl;
		for(int i=0; i<j["angular_etas"].size(); i++)
			std::cout << j["angular_etas"][i] << " " << j["lambdas"][i] << " " << j["zetas"][i] << std::endl;

		for(auto const & d : json.items())
		{
			auto const & key = d.key();
			auto const & data = d.value();
			if (key.find("Layer") != std::string::npos) 
			{
				std::cout << key << ", " << data["Bias"].size() << " " << data["Weight"].size() << std::endl;
				//std::cout << data["Weight"] << std::endl;
			}	
		}
	};

	std::vector<int> get_nodes()
	{
		std::cout << "in get_nodes() nnodes: "; 
		for(auto const & d : nnodes) std::cout << d << " ";
		std::cout << std::endl;

		return nnodes;
	}
};


struct Params
{
	std::vector<float> eta, eta0;
	std::vector<float> rs, rs0;
	std::vector<float> rc;
	const float rc0 = 4.5;
	int num_elems; 

	int feature_size; 

	Params() : num_elems{3}, eta0{0.5, 1.0, 3.0}, rs0{1.0, 2.0, 3.0, 4.0}
	{
		for(int k=0; k<num_elems; k++)
			for(int i=0; i<eta0.size(); i++)
				for(int j=0; j<rs0.size(); j++)
				{
					eta.push_back(eta0[i]);
					rs.push_back(rs0[j]);
					rc.push_back(rc0);
				}

		feature_size = num_elems*eta0.size()*rs0.size();
		assert(feature_size > 0);

		assert(eta.size() == feature_size);
		assert(rs.size() == feature_size);
	}

	Params(AenetParams const & ap)
	{
		//ap.show();
		assert(ap.sfparam.size() == ap.sfenv.size());
		feature_size = ap.sfparam.size();

		for (int i = 0; i<ap.sfparam.size(); i++) 
		{ 
			rc.push_back(ap.sfparam[i][0]);
			eta.push_back(ap.sfparam[i][1]);
			rs.push_back(ap.sfparam[i][2]);
		}; 
	}

	Params(JAXMDParams const & jp) // FIXME currently only radial feature is supposed. 
	{
		//jp.show();
		
		feature_size = jp.feature_size;

		const auto j = jp.json["sym_hyper_parameters"];

		for (int i = 0; i<j["radial_etas"].size(); i++) 
		{ 
			rc.push_back(j["cutoff_distance"]); 
			eta.push_back(j["radial_etas"][i]); 
			rs.push_back(0.0);  // rs is not used in JaxMD?
		}; 
	};
};
