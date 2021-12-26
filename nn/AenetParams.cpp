#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cassert>

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
	double mass; 

	// load_Network
	int nlayers, nnodes_max, Wsize, nvalues;
	std::vector<int> nnodes, f_a, iw, iv;
	std::vector<double> W;
	bool init, evaluated, derivatives; 
	
	// load_Setup
	std::string description; // unused
	std::string atomtype, sftype;
	int nenv, nsf, nsfparam; 
	std::vector<std::string> envtypes;
	double Rc_min, Rc_max;

	std::vector<int> sf;
	std::vector<std::vector<double>> sfparam;
	std::vector<std::vector<int>> sfenv;
	int neval;
	std::vector<double> sfval_min, sfval_max, sfval_avg, sfval_cov;

	std::vector<double> D, value, deriv, work, work2, work3, work4;

	AenetParams(const std::string & parminfo)
	{
		std::stringstream ss(parminfo);
		ss >> element >> filename >> mass; 

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

/*
int main(int argc, char *argv[])
{
	std::vector<AenetParams> p; 
	for(int i=1; i<argc; i++)
		p.push_back(AenetParams(argv[i]));


	for(int i=0; i<p.size(); i++) p[i].show();

}
*/
