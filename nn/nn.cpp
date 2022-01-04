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
#include <tuple>
#include <iomanip>
#include <memory>

#include <cstdlib>
#include <cmath>
#include <cassert>

#include "AenetParams.cpp"

#define MAX_NEIGHBOR 100

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

	Params(AenetParams const &ap)
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
};

struct dist 
{
	float dx,dy,dz;
	dist(float _dx, float _dy, float _dz) : dx(_dx), dy(_dy), dz(_dz){};
};

float apply_pbc(float x, float lattice)
{
    if(x>=0.5*lattice) x-=lattice;
    if(x<-0.5*lattice) x+=lattice;
    return x;
};

struct MDFrame
{
    std::string filename;

    int natoms;
    float lattice[6];

    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;
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

MDFrame read_single_mdframe(std::ifstream &in, std::string _filename="NA")
{
	MDFrame mdframe;
	std::string str;

	mdframe.filename = _filename;

	getss(in) >> mdframe.natoms; 
	getss(in) >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2] >> 
		mdframe.lattice[3] >>  mdframe.lattice[4] >>  mdframe.lattice[5];

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

	getss(in) >> name >> x >> y >> z; 

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

typedef std::tuple<float, int, std::string, dist, bool, int, std::string> tuple_t;
typedef std::vector<std::vector<tuple_t>> nbrlist_t;
typedef std::vector<float> nbrdist_t;

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
		nbrdist.resize(4*mdframe.natoms*mdframe.natoms);

		for (int i=0; i<mdframe.natoms; i++)
		{
			for (int j=0; j<mdframe.natoms; j++)
			{

				float dx = mdframe.x[j] - mdframe.x[i];
				float dy = mdframe.y[j] - mdframe.y[i];
				float dz = mdframe.z[j] - mdframe.z[i];

				dx = apply_pbc(dx, mdframe.lattice[0]);
				dy = apply_pbc(dy, mdframe.lattice[1]);
				dz = apply_pbc(dz, mdframe.lattice[2]);

				float dr = sqrt(dx*dx + dy*dy + dz*dz); 
				if (dr == 0.0) continue;

				auto data = std::make_tuple(dr, j, mdframe.name[j], dist(dx,dy,dz), true, i, mdframe.name[i]);

				nbrlist[i].push_back(data);

				int idx = 4*(i*mdframe.natoms + j); 
				nbrdist[idx]   = dr;
				nbrdist[idx+1] = dx;
				nbrdist[idx+2] = dy;
				nbrdist[idx+3] = dz;
				//std::cout << i << " " << j << " " << dr << " " << dx << " " << dy << " " << dz << std::endl;
			}
		}
	}
};

std::tuple<std::vector<float>,std::vector<float>> 
featurize_nbrlist(int const natoms, int const maxnbrs, std::vector<float> const & nbrdist, Params const &p)
{
	//auto natoms = mdframe.natoms;

	std::vector<float> G2, dG2;

	G2.resize(natoms*p.feature_size,0.0);
	dG2.resize(3*natoms*p.feature_size,0.0);

	std::cout << "natoms,maxnbrs,size(G2,dG2,eta,rs) : " << natoms << " " << maxnbrs << " " << 
		G2.size() << " " << dG2.size() << " " << p.eta.size() << " " << p.rs.size() << std::endl;

	int G2_size = G2.size();
	int dG2_size = dG2.size();
	int eta_size = p.eta.size();
	int rs_size = p.rs.size();
	int rc_size = p.rc.size();
	int feature_size = p.feature_size; 
	int nbrdist_size = nbrdist.size();
	std::cout << "size(G2,dG2,eta,rs,rc,feature_size,nbrdist) : " << G2_size << " " << dG2_size << " " 
		<< eta_size << " " << rs_size << " " << rc_size << " " 
		<< feature_size << " " << nbrdist_size << std::endl;

	auto G2_ptr = G2.data();
	auto dG2_ptr = dG2.data();
	auto eta_ptr = p.eta.data();
	auto rs_ptr = p.rs.data();
	auto rc_ptr = p.rc.data();
	auto nbrdist_ptr = nbrdist.data();

	#pragma omp target data map(tofrom: G2_ptr[0:G2_size], tofrom: dG2_ptr[0:dG2_size]), \
		map(to: eta_ptr[0:feature_size], rs_ptr[0:feature_size], rc_ptr[0:feature_size]), \
		map(to: nbrdist_ptr[0:nbrdist_size]),  map(to: natoms, feature_size, nbrdist_size)
	for (int n=0; n<natoms; n++)
	{
		for (int ii = 0; ii < feature_size; ii++)
		{
			auto rs_val = rs_ptr[ii];
			auto eta_val = eta_ptr[ii];
			auto rc_val = rc_ptr[ii];

			for (int j=0; j<maxnbrs; j++)
			{
				if (n==j) continue;

				int ii4 = 4*(n*maxnbrs + j);
				auto dr = nbrdist_ptr[ii4];
				auto dx = nbrdist_ptr[ii4+1]/dr;
				auto dy = nbrdist_ptr[ii4+2]/dr;
				auto dz = nbrdist_ptr[ii4+3]/dr;

				auto rij_rs = dr - rs_val; 
				auto exp_rij = exp(-eta_val * rij_rs * rij_rs);
				auto fc_rij = 0.5*cos(M_PI*dr/rc_val);

				auto G2_val = exp_rij*fc_rij;
				if (dr > rc_val) G2_val = 0.0;

				int idx = n*feature_size + ii;
				G2_ptr[idx] += G2_val;

				auto G2_deriv = exp_rij*(-2.0*eta_val*fc_rij*rij_rs - M_PI/(2.0*rc_val)*sin(M_PI*dr/rc_val));
				if (dr > rc_val) G2_deriv= 0.0;

				int idx3 = 3*(n*feature_size + ii);
				//std::cout << n << " " << j << " " << dr << " " << rc_val << " " << G2_val << " " << G2_deriv << std::endl;
				dG2_ptr[idx3  ] += G2_deriv*dx;
				dG2_ptr[idx3+1] += G2_deriv*dy;
				dG2_ptr[idx3+2] += G2_deriv*dz;
				/*
				std::cout << n << " " << j << " " << dr << " " << rc_val << " " << 
						ii4 << " " << G2_deriv << " " << dx << " " << dy << " " << dz << " " << 
						G2_ptr[idx] << " " << dG2_ptr[idx3  ] << " " << dG2_ptr[idx3+1] << " " << dG2_ptr[idx3+2] << std::endl;
				*/

			}
		}
	}

	return std::make_tuple(G2,dG2);
}


std::tuple<std::vector<float>,std::vector<float>> 
featurize_nbrlist(int const natoms, int const maxnbrs, void *nbrdist_voidptr, Params const &p)
{
	float *nbrdist = (float *) nbrdist_voidptr; 

	//auto natoms = mdframe.natoms;

	std::vector<float> G2, dG2;

	G2.resize(natoms*p.feature_size,0.0);
	dG2.resize(3*natoms*p.feature_size,0.0);

	std::cout << "natoms,maxnbrs,size(G2,dG2,eta,rs) : " << natoms << " " << maxnbrs << " " 
		<< G2.size() << " " << dG2.size() << " " << p.eta.size() << " " << p.rs.size() << std::endl;

	int G2_size = G2.size();
	int dG2_size = dG2.size();
	int eta_size = p.eta.size();
	int rs_size = p.rs.size();
	int rc_size = p.rc.size();
	int feature_size = p.feature_size; 
	int nbrdist_size = natoms*maxnbrs;

	std::cout << "size(G2,dG2,eta,rs,rc,feature_size,nbrdist) : " << G2_size << " " << dG2_size << " " 
		<< eta_size << " " << rs_size << " " << rc_size << " " << feature_size << " " << nbrdist_size << std::endl;

	auto G2_ptr = G2.data();
	auto dG2_ptr = dG2.data();
	auto eta_ptr = p.eta.data();
	auto rs_ptr = p.rs.data();
	auto rc_ptr = p.rc.data();
	auto nbrdist_ptr = nbrdist;

	#pragma omp target data map(tofrom: G2_ptr[0:G2_size], dG2_ptr[0:dG2_size]), \
		map(to: eta_ptr[0:feature_size], rs_ptr[0:feature_size], rc_ptr[0:feature_size]), \
		map(to: nbrdist_ptr[0:nbrdist_size]),  map(to: natoms, feature_size, nbrdist_size)
	for (int n=0; n<natoms; n++)
	{
		for (int ii = 0; ii < feature_size; ii++)
		{
			auto rs_val = rs_ptr[ii];
			auto eta_val = eta_ptr[ii];
			auto rc_val = rc_ptr[ii];

			for (int j=0; j<maxnbrs; j++)
			{
				int ii4 = 4*(n*maxnbrs + j);

				auto dr = nbrdist_ptr[ii4];
				if(dr == 0.0) continue;

				auto dx = nbrdist_ptr[ii4+1]/dr;
				auto dy = nbrdist_ptr[ii4+2]/dr;
				auto dz = nbrdist_ptr[ii4+3]/dr;

				auto rij_rs = dr - rs_val; 
				auto exp_rij = exp(-eta_val * rij_rs * rij_rs);
				auto fc_rij = 0.5*cos(M_PI*dr/rc_val);

				auto G2_val = exp_rij*fc_rij;
				if (dr > rc_val) G2_val = 0.0;

				int idx = n*feature_size + ii;
				G2_ptr[idx] += G2_val;

				auto G2_deriv = exp_rij*(-2.0*eta_val*fc_rij*rij_rs - M_PI/(2.0*rc_val)*sin(M_PI*dr/rc_val));
				if (dr > rc_val) G2_deriv= 0.0;

				int idx3 = 3*(n*feature_size + ii);
				//std::cout << n << " " << j << " " << dr << " " << rc_val << " " << G2_val << " " << G2_deriv << std::endl;
				dG2_ptr[idx3  ] += G2_deriv*dx;
				dG2_ptr[idx3+1] += G2_deriv*dy;
				dG2_ptr[idx3+2] += G2_deriv*dz;

				/*
				std::cout << n << " " << j << " " << dr << " " << rc_val << " " << 
						ii4 << " " << G2_deriv << " " << dx << " " << dy << " " << dz << " , " << 
						G2_ptr[idx] << " " << dG2_ptr[idx3] << " " << dG2_ptr[idx3+1] << " " << dG2_ptr[idx3+2] << std::endl;
				*/

			}
		}
	}

	return std::make_tuple(G2,dG2);
}

struct Dense
{
	int width, height, n; 
	int nodes; // nodes per layer. height == nodes

	std::vector<float> w, b; 

	Dense(int _width, int _height, bool randomize = true)
		: width(_width), height(_height), nodes(_height)
	{
		assert(height>0); assert(width>0); 

		n = height*width;

		if(randomize)
		{
			b.resize(height);
			w.resize(height*width);
			for(int i=0; i < b.size(); i++) b[i] = 2.0*((float)rand()/RAND_MAX-0.5);
			for(int i=0; i < w.size(); i++) w[i] = 2.0*((float)rand()/RAND_MAX-0.5);
		}
	}

	void set_w(std::vector<float> & _w)
	{
		w.resize(0);
		std::copy(_w.begin(), _w.end(), back_inserter(w));
	}

	void set_b(std::vector<float> & _b)
	{
		b.resize(0);
		std::copy(_b.begin(), _b.end(), back_inserter(b));
	}

};

float tanh(float x)
{
	float e2x = std::exp(2*x);
	return (e2x - 1.0)/(e2x + 1.0);
}

float dtanh(float x)
{
	return 1.0 - tanh(x);
}

float relu(float x)
{
	return x < 0.0 ? 0.0 : x;
};

float drelu(float x)
{
	return x < 0.0 ? 0.0 : 1.0;
};

float actfunc(float x)
{
	return tanh(x);	
}

float dactfunc(float x)
{
	return dtanh(x);
}

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

	void set_wb_aenet(const std::vector<float> & W, const std::vector<int> & nnodes)
	{
		//std::cout << nnodes.size() << " " << nn.size() << std::endl;
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

	std::tuple<std::vector<float>,std::vector<float>> 
		predict(const std::vector<float> & feature, std::vector<float> & dfeature, const int batch_size=1)
	{
		std::cout << "\n\n   Entering Net::predict()   \n\n";
		std::cout << "size(batch,feature,nn),nn[0].width : " << batch_size << " " 
			<< feature.size() << " " << nn.size() << " " << nn[0].width << std::endl; 

		//assert(feature.size() == batch_size*nn[0].width); 

		std::vector<float> vin,vout; 
		std::vector<float> dvin,dvout,dsigma; 
		std::vector<float> energies, forces; 

		std::copy(feature.begin(), feature.end(), back_inserter(vin));
		std::copy(dfeature.begin(), dfeature.end(), back_inserter(dvin));

		/*
		std::cout << "vin ===============================\n";
		for(int ih=0; ih<vin.size(); ih++) std::cout << vin[ih] << " "; std::cout << std::endl;
		std::cout << "===============================\n";
		std::cout << "dvin ===============================\n";
		for(int ih=0; ih<dvin.size(); ih++) std::cout << dvin[ih] << " "; std::cout << std::endl;
		std::cout << "===============================\n";
		*/

		for(int n=0; n<batch_size; n++)
		{
			const int fsize_per_particle = nn[0].width; 
			vin.resize(fsize_per_particle, 0.0);
			dvin.resize(3*fsize_per_particle, 0.0);

			for(int i=0; i<fsize_per_particle; i++) vin[i]=feature[n*fsize_per_particle+i];
			for(int i=0; i<3*fsize_per_particle; i++) dvin[i]=dfeature[n*3*fsize_per_particle+i];

			for(int i=0; i<nn.size(); i++)
			{
				const int height = nn[i].height;
				const int width = nn[i].width;
	
				//vout.resize(height);
				//for(int ih=0; ih<height; ih++) vout[ih]=nn[i].b[ih];
	
				vout.resize(0);
				std::copy(nn[i].b.begin(), nn[i].b.end(), back_inserter(vout));
	
//std::cout << "ilayer,vin,vout: " << i << " " << vin.size() << " " << vout.size() << std::endl;
	
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
				for(int ih=0; ih<height; ih++) dsigma[ih]=dactfunc(vout[ih]);
				for(int ih=0; ih<height; ih++) vout[ih]=actfunc(vout[ih]);
/*
std::cout << "=============================== vout after activation \n";
for(int ih=0; ih<height; ih++) std::cout << vout[ih] << " "; std::cout << std::endl;
std::cout << "===============================\n";
*/
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
			std::copy(vin.begin(), vin.end(), back_inserter(energies));
			std::copy(dvin.begin(), dvin.end(), back_inserter(forces));
		}

		return std::make_tuple(energies,forces);
	}
};

struct RXMDNN
{
	MDFrame md; 
	NeighborList nbr; 

	std::vector<AenetParams> aps;
	std::vector<Params> pms;

	std::vector<Net> nns; 

	int natoms; 

	// use XYZ file to construct MDFrame and NeighborList (for testing purpose)
	RXMDNN(std::string const & structfile, std::string const & paramfile) 
	{
		std::string filename(structfile);
		std::ifstream fin(filename);
	
		md = read_single_mdframe(fin);
		md.print();

		nbr = NeighborList(md);
	
		std::string paramfiles(paramfile);
		std::ifstream pin(paramfiles);

		std::string line; 
		while(std::getline(pin,line))
		{
			if(line.compare(0,5,"model") != 0) continue;

			auto ap = AenetParams(line);
			aps.push_back(ap);
	
			auto pm = Params(ap);
			pms.push_back(pm);

			std::vector<int> nodes = ap.get_nodes();

			auto nn = Net(nodes);
			nn.set_wb_aenet(ap.W, ap.nnodes);

			nns.push_back(nn);
		}

	}

	// Setup model parameters. no XYZ is used. 
	RXMDNN(int const _natoms, std::string const & paramfile) 
	{
		natoms = _natoms; 
		std::string paramfiles(paramfile);
		std::ifstream pin(paramfiles);

		std::string line; 
		while(std::getline(pin,line))
		{
			if(line.compare(0,5,"model") != 0) continue;

			auto ap = AenetParams(line);
			aps.push_back(ap);
	
			auto pm = Params(ap);
			pms.push_back(pm);

			std::vector<int> nodes = ap.get_nodes();

			auto nn = Net(nodes);
			nn.set_wb_aenet(ap.W, ap.nnodes);

			nns.push_back(nn);
		}

	}

	float get_maxrc_rxmdnn()
	{
		float maxrc=0.0; 
		for(auto pm : pms) for(auto rc : pm.rc) if(maxrc < rc) maxrc = rc;

		return maxrc; 
	}

	//void predict(int const natoms, int const maxnbrs, std::vector<float> const & nbrdist)
	void predict(int const natoms, int const atom_type, int const maxnbrs, void* nbrdist_ptr)
	{
		std::cout << "\n\nEntering RXMDNN::predict()\n\nnatoms,atom_type,maxnbrs: " << 
			natoms << " " << atom_type << " " << maxnbrs << std::endl;

		//atom_type is 1-indexed while pms&&nn are 0-indexed. 
		int atom_type_idx = atom_type - 1; 

		auto & pm = pms[atom_type_idx]; 

		auto feature = featurize_nbrlist(natoms, maxnbrs, nbrdist_ptr, pm);
		auto G2 = std::get<0>(feature);
		auto dG2 = std::get<1>(feature);

		/*
		std::cout << "\nG2 size " << G2.size() << std::endl;
		for(int j=0; j<natoms; j++) 
		{
			std::cout << std::setw(3) << j << ": ";
			//for(int i=0; i<pm.feature_size; i++) std::cout << std::scientific << 
			for(int i=0; i<10; i++) std::cout << std::scientific << 
				std::setprecision(5) << std::setw(10) << G2[j*pm.feature_size+i] << " "; 
			std::cout << std::endl;
		}
		std::cout << "\ndG2 size " << dG2.size() << std::endl;
		for(int j=0; j<natoms; j++) 
		{
			std::cout << std::setw(3) << j << ": ";
			//for(int i=0; i<pm.feature_size; i++) 
			for(int i=0; i<10; i++) 
				std::cout << std::scientific << std::setprecision(5) << std::setw(10) << dG2[3*(j*pm.feature_size+i)] << " "; 
			std::cout << std::endl;
		}
		std::cout << std::endl;
		*/

		auto result = nns[atom_type_idx].predict(G2, dG2, natoms); 

		auto energy = std::get<0>(result);
		auto force = std::get<1>(result);

		std::cout << "size(energy,force) : " << energy.size() << " " << force.size() << std::endl;
		for(int i=0; i<energy.size(); i++)
			std::cout << i << " " << energy[i] << " " << force[3*i] << " " << force[3*i+1] << " " << force[3*i+2] << std::endl;
	};
};

std::unique_ptr<RXMDNN> rxmdnn_ptr; 

extern "C" void init_rxmdnn_hybrid(int natoms)
{
	std::cout << "foo from init\n";
	rxmdnn_ptr = std::make_unique<RXMDNN>(natoms,"rxmdnn.in");
}

extern "C" void predict_rxmdnn_hybrid(int natoms, int atom_type, int maxnbrs, void *nbrdist_ptr)
{
	std::cout << "foo from predict_hybrid\n";
	rxmdnn_ptr->predict(natoms, atom_type, maxnbrs, nbrdist_ptr);
}

extern "C" void init_rxmdnn(void)
{
	std::cout << "foo from init\n";
	rxmdnn_ptr = std::make_unique<RXMDNN>("pto.xyz","rxmdnn.in");
}

extern "C" void predict_rxmdnn(void)
{
	std::cout << "foo from predict\n";
	int atom_type = 1;
	rxmdnn_ptr->predict(rxmdnn_ptr->md.natoms, atom_type, rxmdnn_ptr->md.natoms, (void *) rxmdnn_ptr->nbr.nbrdist.data()); 
}

extern "C" void get_maxrc_rxmdnn(double & maxrc)
{
	maxrc = rxmdnn_ptr->get_maxrc_rxmdnn(); 
	std::cout << "foo from maxrc " << maxrc << std::endl;
}
