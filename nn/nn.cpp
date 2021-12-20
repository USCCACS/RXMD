#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <set>
#include <regex>

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "AenetParams.cpp"

#define MAX_NEIGHBOR 100

struct dist 
{
	double dx,dy,dz;
	dist(double _dx, double _dy, double _dz) : dx(_dx), dy(_dy), dz(_dz){};
};

double apply_pbc(double x, double lattice)
{
    if(x>=0.5*lattice) x-=lattice;
    if(x<-0.5*lattice) x+=lattice;
    return x;
};

struct Params
{
	std::vector<double> eta, eta0;
	std::vector<double> rs, rs0;
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
				}

		feature_size = num_elems*eta0.size()*rs0.size();
		assert(feature_size > 0);

		assert(eta.size() == feature_size);
		assert(rs.size() == feature_size);
	}
};

struct MDFrame
{
    std::string filename;

    int natoms;
    double lattice[6];

    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<int> mol_id;  // molecular Id
    std::vector<std::string> name;
    std::map<std::string,int> elems;

    void print()
    {
      std::cout << " # of Atoms : " << natoms << std::endl;
      std::cout << " Lattice Consts : " <<
              lattice[0] << " " << lattice[1] << " " << lattice[2] << " " <<
              lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;

      for (auto e : elems)
      	std::cout << e.first << " - " << e.second << " " ;
      std::cout << std::endl;

    }
};

typedef std::tuple<double, int, std::string, dist, bool, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;
typedef std::vector<double> nbrdist_t;

struct NeighborList
{
	int max_neighbors;
	int num_atoms; 

	nbrlist_t nbrlist;
	nbrdist_t nbrdist;

	NeighborList() {}; 

	NeighborList(MDFrame & mdframe, int _max=MAX_NEIGHBOR) : max_neighbors(_max), num_atoms(mdframe.natoms)
	{
		nbrlist.resize(mdframe.natoms);
		nbrdist.resize(mdframe.natoms*mdframe.natoms);

		for (int i=0; i<mdframe.natoms; i++)
		{
			for (int j=0; j<mdframe.natoms; j++)
			{

				double dx = mdframe.x[j] - mdframe.x[i];
				double dy = mdframe.y[j] - mdframe.y[i];
				double dz = mdframe.z[j] - mdframe.z[i];

				dx = apply_pbc(dx, mdframe.lattice[0]);
				dy = apply_pbc(dy, mdframe.lattice[1]);
				dz = apply_pbc(dz, mdframe.lattice[2]);

				double dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				auto data = std::make_tuple(dr, j, mdframe.name[j], dist(dx,dy,dz), true, i, mdframe.name[i]);

				nbrlist[i].push_back(data);
				nbrdist[i*mdframe.natoms + j] = dr;
				//std::cout << i << " " << j << " " << dr << std::endl;
			}
		}
	}
};

std::vector<double> featurize_nbrlist(MDFrame &mdframe, NeighborList &nbr, Params &p)
{
	auto natoms = mdframe.natoms;

	std::vector<double> G2;

	G2.resize(p.feature_size,0.0);

	std::cout << "eta.rs,G2 size : " << G2.size() << " " << p.eta.size() << " " << p.rs.size() << std::endl;

	int G2_size = G2.size();
	int eta_size = p.eta.size();
	int rs_size = p.rs.size();
	int feature_size = p.feature_size; 
	int nbrdist_size = nbr.nbrdist.size();
	std::cout << "size G2, eta, rs, feature_size, nbrdist : " << G2_size << " " 
		<< eta_size << " " << rs_size << " " << feature_size << " " << nbrdist_size << std::endl;

	auto G2_ptr = G2.data();
	auto eta_ptr = p.eta.data();
	auto rs_ptr = p.rs.data();
	auto nbrdist_ptr = nbr.nbrdist.data();

	#pragma omp target data map(tofrom: G2_ptr[0:G2_size]), \
	map(to: eta_ptr[0:feature_size]), \
       	map(to: rs_ptr[0:feature_size]), \
       	map(to: nbrdist_ptr[0:nbrdist_size]), \
	map(to: natoms, feature_size, nbrdist_size)
	for (int i=0; i<natoms; i++)
	{
		for (int ii = 0; ii < feature_size; ii++)
		{
			auto rs_val = rs_ptr[ii];
			auto eta_val = eta_ptr[ii];

			for (int j=0; j<natoms; j++)
			{
				auto dr = nbrdist_ptr[i*natoms + j] - rs_val; 
				G2_ptr[ii] += exp(-eta_val*dr*dr);
			}
		}
	}

	return G2;
}

struct Dense
{
	int width, height, n; 
	int nodes; // nodes per layer. height == nodes

	std::vector<double> w, b; 

	Dense(int _width, int _height, bool randomize = true)
		: width(_width), height(_height), nodes(_height)
	{
		assert(height>0); assert(width>0); 

		n = height*width;

		if(randomize)
		{
			b.resize(height);
			w.resize(height*width);
			for(int i=0; i < b.size(); i++) b[i] = 2.0*((double)rand()/RAND_MAX-0.5);
			for(int i=0; i < w.size(); i++) w[i] = 2.0*((double)rand()/RAND_MAX-0.5);
		}
	}

	void set_w(std::vector<double> & _w)
	{
		w.resize(0);
		std::copy(_w.begin(), _w.end(), back_inserter(w));
	}

	void set_b(std::vector<double> & _b)
	{
		b.resize(0);
		std::copy(_b.begin(), _b.end(), back_inserter(b));
	}

};

double relu(double x)
{
	return x < 0.0 ? 0.0 : x;
};

double drelu(double x)
{
	return x < 0.0 ? 0.0 : 1.0;
};

struct Net
{
	std::vector<Dense> nn; 

	Net(const std::vector<int> & dims)
	{
		for(int i=0; i<dims.size()-1; i++)
		{
			assert(dims[i]>0); 
			assert(dims[i+1]>0); 

			nn.push_back(Dense(dims[i],dims[i+1],true));
			std::cout << "added layer " << i << " : " << dims[i] << " x " << dims[i+1] << "   " <<
				"wsize " << nn[i].w.size()<< " bsize " << nn[i].b.size() << std::endl;
		}
	}

	void set_wb_aenet(const std::vector<double> & W, const std::vector<int> & nnodes)
	{
		std::cout << nnodes.size() << " " << nn.size() << std::endl;
		assert(nnodes.size()-1 == nn.size());

		int iw=0;
		for(int i=0; i<nn.size()-1; i++)
		{
			int width = nnodes[i];
			int height = nnodes[i+1];

			nn[i].w.resize(width*height);
			nn[i].b.resize(height);

			for(int n2=0; n2<height; n2++)
			{
				for(int n1=0; n1<width; n1++) nn[i].w[n2*width + n1] = W[iw++];
				nn[i].b[n2] = W[iw++];
			}
		}
	}

	std::vector<double> predict(const std::vector<double> & feature, std::vector<double> & dfeature)
	{
		std::cout << feature.size() << " " << nn.size() << " " << nn[0].width << std::endl; 

		//assert(feature.size() == nn[nn.size()-1].width); 
		assert(feature.size() == nn[0].width); 

		std::vector<double> vin,vout; 
		std::vector<double> dvin,dvout,dsigma; 

		std::copy(feature.begin(), feature.end(), back_inserter(vin));
		std::copy(dfeature.begin(), dfeature.end(), back_inserter(dvin));

		/*
		std::cout << "vin ===============================\n";
		for(int ih=0; ih<vin.size(); ih++) std::cout << vin[ih] << " "; std::cout << std::endl;
		std::cout << "===============================\n";
		std::cout << "vout ===============================\n";
		for(int ih=0; ih<dvin.size(); ih++) std::cout << dvin[ih] << " "; std::cout << std::endl;
		std::cout << "===============================\n";
		*/

		//for(int i=nn.size()-1; 0<=i; i--)
		for(int i=0; i<nn.size(); i++)
		{
			const int height = nn[i].height;
			const int width = nn[i].width;

			//vout.resize(height);
			//for(int ih=0; ih<height; ih++) vout[ih]=nn[i].b[ih];

			vout.resize(0);
			std::copy(nn[i].b.begin(), nn[i].b.end(), back_inserter(vout));

			std::cout << "ilayer,vin,vout: " 
				<< i << " " << vin.size() << " " << vout.size() << std::endl;

			/*
			std::cout << "vin ===============================\n";
			for(int ih=0; ih<vin.size(); ih++) std::cout << vin[ih] << " "; std::cout << std::endl;
			std::cout << "===============================\n";

			std::cout << "vout before ===============================\n";
			for(int ih=0; ih<vout.size(); ih++) std::cout << vout[ih] << " "; std::cout << std::endl;
			std::cout << "===============================\n";
			*/


			for(int ih=0; ih<height; ih++)
			{
			for(int iw=0; iw<width; iw++)
			{
				vout[ih] += nn[i].w[ih*width+iw]*vin[iw];
			        //std::cout << ih << " " << iw << " " << 
				//	vout[ih] << " " << nn[i].w[ih*width+iw] << " " << vin[iw] << std::endl;
			}
			}
			//std::cout << std::endl;

			/*
			std::cout << "nn.w ===============================\n";
			for(int ih=0; ih<height; ih++)
			for(int iw=0; iw<width; iw++)
				std::cout <<  nn[i].w[ih*width+iw] << " ";
			std::cout << std::endl;
			std::cout << "===============================\n";

			std::cout << "vout after ===============================\n";
			for(int ih=0; ih<height; ih++) std::cout << vout[ih] << " "; std::cout << std::endl;
			std::cout << "===============================\n";
			*/

			// compute sigma_deriv(z) first, then appy activation function sigma(z)
			dsigma.resize(height);
			for(int ih=0; ih<height; ih++) dsigma[ih]=drelu(vout[ih]);
			for(int ih=0; ih<height; ih++) vout[ih]=relu(vout[ih]);

			std::cout << "=============================== vout after activation \n";
			for(int ih=0; ih<height; ih++) std::cout << vout[ih] << " "; std::cout << std::endl;
			std::cout << "===============================\n";

			dvout.resize(3*height,0.0);

			for(int ih=0; ih<height; ih++)
			for(int iw=0; iw<width; iw++)
			for(int ii=0; ii<3; ii++) 
				dvout[3*ih+ii] += nn[i].w[ih*width + iw]*dvin[3*iw+ii];

			for(int ih=0; ih<height; ih++) 
			for(int ii=0; ii<3; ii++) dvout[3*ih+ii]*=dsigma[ih];

			vin.resize(0);
			std::copy(vout.begin(), vout.end(), back_inserter(vin));

			dvin.resize(0);
			std::copy(dvout.begin(), dvout.end(), back_inserter(dvin));

		}

		return vin;
	}
};



MDFrame read_single_mdframe(std::ifstream &in, std::string _filename="NA")
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
    mdframe.filename = _filename;
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;

    ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2];
    mdframe.lattice[3] = mdframe.lattice[4] = mdframe.lattice[5] = 90.0;

    mdframe.name.resize(mdframe.natoms);
    mdframe.x.resize(mdframe.natoms);
    mdframe.y.resize(mdframe.natoms);
    mdframe.z.resize(mdframe.natoms);
    mdframe.mol_id.resize(mdframe.natoms);

    for (int i=0; i<mdframe.natoms; i++)
    {
		std::string name;
		float x,y,z,vx,vy,vz;
		int id;

		std::stringstream ss;
		std::getline(in,str);
		ss << str;

        //ss >> name >> x >> y >> z >> dummy >> id;  
		ss >> name >> x >> y >> z; 
		id = i+1;

		mdframe.name[id-1] = name;
		mdframe.x[id-1] = x;
		mdframe.y[id-1] = y;
		mdframe.z[id-1] = z;
		mdframe.mol_id[id-1] = id;

		if (mdframe.elems.count(name) == 0)
		{
			mdframe.elems[name]=1;
		} else {
			mdframe.elems[name]++;
		}
	}

    return mdframe;
};

int main(int argc, char* argv[])
{
	std::string filename(argv[1]);
	std::ifstream fin(filename);

	AenetParams ap =  AenetParams(argv[2]);
	ap.show();

	std::vector<int> nodes = ap.get_nodes();

	auto nn = Net(nodes);
	nn.set_wb_aenet(ap.W, ap.nnodes);

	auto md = read_single_mdframe(fin);
	md.print();
	auto param = Params();
	auto nbr = NeighborList(md);
	auto G2 = featurize_nbrlist(md, nbr, param);
	std::vector<double> dG2;
	dG2.resize(3*G2.size(),0.1); // 3 derivatives dE/dxyz
	auto result = nn.predict(G2, dG2); 

	std::cout << "result.size() " << result.size() << std::endl;
	for(auto r : result) std::cout << r << " ";
	std::cout << std::endl;

}
