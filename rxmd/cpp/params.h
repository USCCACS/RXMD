#include <stdint.h>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

struct Params
{
    public:

        Params(const std::string ffFileName)
        {
            std::ifstream ffFile(ffFileName.c_str());

			// read header of force field file
            std::getline(ffFile, ffHeader); 

			std::string line; 

			std::stringstream ss;

// atom-type independent params
            std::getline(ffFile, line); 
			ss.str(""); ss << line; 
			ss >> npar;

			vpar.resize(npar,-1.0); 
			for(int i=0; i<npar; ++i) 
			{
            	std::getline(ffFile, line); 
				ss.str(""); ss << line; 
				ss >> vpar[i];
			}

			// Constant parameters reset to actual use
			pvdW1 = vpar[28]; 
			pvdW1h = 0.5*pvdW1; 
			pvdW1inv = 1.0/pvdW1;

			// a small modification in sigma-bond prime
			vpar30 = vpar[29];

// single-atom depend params,  Nr of types of atoms
            std::getline(ffFile, line); 
			ss.str(""); ss << line; 
			ss >> nso;

			inxn2.resize(nso);
			inxn3.resize(nso);
			inxn3hb.resize(nso);
			inxn4.resize(nso);
			for(int i1=0; i1<nso; ++i1)
			{
				inxn2[i1].resize(nso,0);
				inxn3[i1].resize(nso);
				inxn3hb[i1].resize(nso);
				inxn4[i1].resize(nso);
				for(int i2=0; i2<nso; ++i2)
				{	
					inxn3[i1][i2].resize(nso,0);
					inxn3hb[i1][i2].resize(nso,0);
					inxn4[i1][i2].resize(nso);
					for(int i3=0; i3<nso; ++i3)
						inxn4[i1][i2][i3].resize(nso,0);
				}
			}

			atmname.resize(nso,""); rat.resize(nso,0); Val.resize(nso,0); mass.resize(nso,0);
			rvdw1.resize(nso,0); eps.resize(nso,0); gam.resize(nso,0); rapt.resize(nso,0);
			Vale.resize(nso,0); alf.resize(nso,0); vop.resize(nso,0); Valboc.resize(nso,0);
			povun5.resize(nso,0); chi.resize(nso,0); eta.resize(nso,0); vnq.resize(nso,0);
			plp2.resize(nso,0); bo131.resize(nso,0); bo132.resize(nso,0); bo133.resize(nso,0);
			povun2.resize(nso,0); pval3.resize(nso,0); Valval.resize(nso,0); pval5.resize(nso,0);

			// Skip 3 lines
			for(int iskip=0; iskip<3; ++iskip)
			{
	            std::getline(ffFile, line); 
			}

			double dnull; // dummy variable 
			for(int i1=0; i1<nso; ++i1)
			{
            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> atmname[i1] >> rat[i1] >> Val[i1] >> mass[i1] >> 
					rvdw1[i1] >> eps[i1] >> gam[i1] >> rapt[i1] >> Vale[i1]; 

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >>  alf[i1] >> vop[i1] >> Valboc[i1] >> povun5[i1] >> dnull >> chi[i1] >> eta[i1] >> dnull;

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> vnq[i1] >> plp2[i1] >> dnull >> bo131[i1] >> bo132[i1] >> bo133[i1] >> dnull >> dnull;

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> povun2[i1] >> pval3[i1] >> dnull >> Valval[i1] >> pval5[i1];
			}

			// allocate derived variables
			nlpopt.resize(nso,0); Valangle.resize(nso,0);

			r0s.resize(nso); r0p.resize(nso); r0pp.resize(nso);
			rvdW.resize(nso); Dij.resize(nso); alpij.resize(nso);
			gamW.resize(nso); gamij.resize(nso);

			for(int i1=0; i1<nso; ++i1)
			{
				r0s[i1].resize(nso,0); r0p[i1].resize(nso,0); r0pp[i1].resize(nso,0);
				rvdW[i1].resize(nso,0); Dij[i1].resize(nso,0); alpij[i1].resize(nso,0);
				gamW[i1].resize(nso,0); gamij[i1].resize(nso,0);

				nlpopt[i1] = 0.5*(Vale[i1] - Val[i1]);
				Valangle[i1] = Valboc[i1];

				// Calc default r0s, r0p, r0pp:
				for(int i2=0; i2<nso; ++i2)
				{
					// Terms for the Bond Order Calculation:
					r0s[i1][i2] = 0.5*(rat[i1]+rat[i2]);
					r0p[i1][i2] = 0.5*(rapt[i1]+rapt[i2]);
					r0pp[i1][i2] = 0.5*(vnq[i1]+vnq[i2]);

					// Terms used in van der Waals calc: 
					rvdW[i1][i2] = sqrt( 4.0*rvdw1[i1]*rvdw1[i2] );
					Dij[i1][i2] = sqrt( eps[i1]*eps[i2] );
					alpij[i1][i2] = sqrt( alf[i1]*alf[i2] );
					gamW[i1][i2] = sqrt( vop[i1]*vop[i2] );
					gamij[i1][i2] = pow( gam[i1]*gam[i2] , -1.5); // gamcco in reac.f
				}
			}

// 2-atom dependent parameters:
            std::getline(ffFile, line);
			ss.str(""); ss << line; 
			ss >> nboty; 

			pbo1.resize(nboty);pbo2.resize(nboty);pbo3.resize(nboty);
			pbo4.resize(nboty);pbo5.resize(nboty);pbo6.resize(nboty);bom.resize(nboty);
			pboc1.resize(nboty);pboc2.resize(nboty);pboc3.resize(nboty);pboc4.resize(nboty);pboc5.resize(nboty);
			Desig.resize(nboty); Depi.resize(nboty);Depipi.resize(nboty);pbe1.resize(nboty);pbe2.resize(nboty);
			povun1.resize(nboty);ovc.resize(nboty); v13cor.resize(nboty);

			// skip one line
            std::getline(ffFile, line);

			for(int ih=0; ih<nboty; ++ih)
			{
				int typea, typeb; 

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> typea >> typeb >> Desig[ih] >> Depi[ih] >> Depipi[ih] >> 
					pbe1[ih] >> pbo5[ih] >> v13cor[ih] >> pbo6[ih] >> povun1[ih]; 

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss  >>  pbe2[ih] >> pbo3[ih] >> pbo4[ih] >> bom[ih] >> pbo1[ih] >> pbo2[ih] >> ovc[ih] >> dnull; 

				// NOTE: the original ReaxFF uses 1-based index for atom type,
				// while this code follows the C++ convention, i.e. 0-based index. 
				typea--; typeb--;

				inxn2[typea][typeb] = ih; 
				inxn2[typeb][typea] = ih; 
			}

			// required by input file backwards setup
			for(int ih=0; ih<nboty; ++ih)
			{
				pboc1[ih] = vpar[1];
				pboc2[ih] = vpar[2];
			}

			for(int i1=0; i1<nso; ++i1)
			{
				for(int i2=0; i2<nso; ++i2)
				{
					int inxn=inxn2[i1][i2];
					pboc3[inxn] = sqrt(bo132[i1]*bo132[i2]);  //  be careful about variable name
					pboc4[inxn] = sqrt(bo131[i1]*bo131[i2]);  //  bo132 -> pboc2, bo131->pbo4  
					pboc5[inxn] = sqrt(bo133[i1]*bo133[i2]);
				}
			}

			// off-diagonal term corrections:
            std::getline(ffFile, line);
			ss.str(""); ss << line; 
			ss >> nodmty; 

			for(int i2=0; i2<nodmty; ++i2)
			{
				int nodm1,nodm2; 
				double deodmh,rodmh,godmh,rsig,rpi,rpi2; 

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> nodm1 >> nodm2 >> deodmh >> rodmh >> godmh >> rsig >> rpi >> rpi2; 

				// NOTE: the original ReaxFF uses 1-based index for atom type,
				// while this code follows the C++ convention, i.e. 0-based index. 
				nodm1--; nodm2--;

				if (rsig > 0.0) r0s[nodm1][nodm2]=rsig; 
				if (rsig > 0.0) r0s[nodm2][nodm1]=rsig; 
				if (rpi > 0.0)  r0p[nodm1][nodm2]=rpi; 
				if (rpi > 0.0)  r0p[nodm2][nodm1]=rpi; 
				if (rpi2 > 0.0) r0pp[nodm1][nodm2]=rpi2; 
				if (rpi2 > 0.0) r0pp[nodm2][nodm1]=rpi2; 
				if (rodmh > 0.0) rvdW[nodm1][nodm2]=2.0*rodmh; 
				if (rodmh > 0.0) rvdW[nodm2][nodm1]=2.0*rodmh; 
				if (deodmh > 0.0) Dij[nodm1][nodm2]=deodmh; 
				if (deodmh > 0.0) Dij[nodm2][nodm1]=deodmh; 
				if (godmh > 0.0) alpij[nodm1][nodm2]=godmh; 
				if (godmh > 0.0) alpij[nodm2][nodm1]=godmh; 
			}

// Valency term params 
            std::getline(ffFile, line);
			ss.str(""); ss << line; 
			ss >> nvaty; 

			pval1.resize(nvaty,0),pval2.resize(nvaty,0),pval4.resize(nvaty,0),pval6.resize(nvaty,0);
			pval7.resize(nvaty,0),pval8.resize(nvaty,0),pval9.resize(nvaty,0),pval10.resize(nvaty,0);
			theta00.resize(nvaty,0);
			ppen1.resize(nvaty,0),ppen2.resize(nvaty,0),ppen3.resize(nvaty,0),ppen4.resize(nvaty,0);
			pcoa1.resize(nvaty,0),pcoa2.resize(nvaty,0),pcoa3.resize(nvaty,0),pcoa4.resize(nvaty,0);

			for(int i=0; i<nvaty; ++i)
			{
				int i1,i2,i3;

            	std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> i1 >> i2 >> i3 >> theta00[i] >> 
					pval1[i] >> pval2[i] >> pcoa1[i] >> pval7[i] >> ppen1[i] >> pval4[i]; 

				// NOTE: the original ReaxFF uses 1-based index for atom type,
				// while this code follows the C++ convention, i.e. 0-based index. 
				i1--; i2--; i3--;

				inxn3[i1][i2][i3] = i;
				inxn3[i3][i2][i1] = i;
			}

			// Corrections in the Valence terms:
			
			inxn3[0][1][1]=0; //react.f, line 4933
			inxn3[1][1][0]=0; //react.f, line 4933

			for(int i=0; i<nvaty; ++i)
			{
			// Valency Terms which do not depend on inxn type:
				pval6[i] = vpar[14];
				pval8[i] = vpar[33];
				pval9[i] = vpar[16];
				pval10[i] = vpar[17];
			// Penalty Terms which do not depend on inxn type:
				ppen2[i] = vpar[19];
				ppen3[i] = vpar[20];
				ppen4[i] = vpar[21];
			// 3body Conjugation Terms which do not depend on type:
				pcoa2[i] = vpar[2];
				pcoa3[i] = vpar[38];
				pcoa4[i] = vpar[30];
			// theta00 given in degrees, but used in radians. Convert by:
				theta00[i] = (M_PI/180.0)*theta00[i];
			}

// 4-body term params
            std::getline(ffFile, line);
			ss.str(""); ss << line; 
			ss >> ntoty; 

			ptor1.resize(ntoty,0),ptor2.resize(ntoty,0),ptor3.resize(ntoty,0),ptor4.resize(ntoty,0);
			V1.resize(ntoty,0),V2.resize(ntoty,0),V3.resize(ntoty,0);
			pcot1.resize(ntoty,0),pcot2.resize(ntoty,0);

			for(int i=0; i<ntoty; ++i)
			{
				int i1,i2,i3,i4;
            	std::getline(ffFile, line);
				ss.str(""); ss << line; 

				ss >> i1 >> i2 >> i3 >> i4 >> V1[i] >> V2[i] >> V3[i] >> 
					ptor1[i] >> pcot1[i] >> dnull >> dnull; 

				// NOTE: the original ReaxFF uses 1-based index for atom type,
				// while this code follows the C++ convention, i.e. 0-based index. 
				i1--; i2--; i3--; i4--;

				// Set up inxn4 lookup table using a condensed input format.
				// If i1 and i4 are zero in the original format (-1 in this code), 
				// i1 and i4 represent all atom types (i.e. wildcard).
				if(i1==-1 && i4==-1) 
				{
					for(int i1 = 0; i1 < nso; ++i1)
					{
						for(int i4 = 0; i4 < nso; ++i4)
						{
							inxn4[i1][i2][i3][i4] = i;
							inxn4[i4][i2][i3][i1] = i;
							inxn4[i1][i3][i2][i4] = i;
							inxn4[i4][i3][i2][i1] = i;
						}
					}
				}
				else if(i1==-1 ||  i4==-1) 
				{
					std::cerr << "\nERROR: asymmetry was fund in force field file. \n" << line << std::endl;
				}
				else
				{
					inxn4[i1][i2][i3][i4] = i;
					inxn4[i4][i2][i3][i1] = i;
					inxn4[i1][i3][i2][i4] = i;
					inxn4[i4][i3][i2][i1] = i;
				}
			}

			// and a few params that don't depend on atom type 
			for(int i1=0; i1<ntoty; ++i1)
			{
				ptor2[i1] = vpar[23];
				ptor3[i1] = vpar[24];
				ptor4[i1] = vpar[25];
				pcot2[i1] = vpar[27];
			}

// hydrogen-bond params
            std::getline(ffFile, line);
			ss.str(""); ss << line; 
			ss >> nhbty; 

			phb1.resize(nhbty);phb2.resize(nhbty);phb3.resize(nhbty);r0hb.resize(nhbty);

			for(int i=0; i<nhbty; ++i)
			{
				int i1,i2,i3;
				std::getline(ffFile, line);
				ss.str(""); ss << line; 
				ss >> i1 >> i2 >> i3 >> r0hb[i] >> phb1[i] >> phb2[i] >> phb3[i];

				// NOTE: the original ReaxFF uses 1-based index for atom type,
				// while this code follows the C++ convention, i.e. 0-based index. 
				i1--; i2--; i3--; 

				inxn3hb[i1][i2][i3] = i;    // Note: inxn3hb(i,j,k) /= inxn3hb(k,j,i)
			}

// close force field file
			ffFile.close();
        };


	public:

		typedef std::vector<double> d1var; 
		typedef std::vector<d1var>  d2var; 

		typedef std::vector<int> i1var; 
		typedef std::vector<i1var>  i2var; 
		typedef std::vector<i2var>  i3var; 
		typedef std::vector<i3var>  i4var; 

		std::string ffHeader; 

		i2var inxn2;
        i3var inxn3, inxn3hb;
        i4var inxn4; 

// Parameters that count number of entries in each field:
		int32_t nso;		// # of atom types
		int32_t nboty;		// # of bond terms
		int32_t nodmty;		// # of off-diag terms
		int32_t npar, nvaty, ntoty, nhbty;

		d1var vpar; 

// Constant parameters reset to actual use
		double pvdW1, pvdW1h, pvdW1inv; 

// a small modification in sigma-bond prime
		double vpar30;

// Parameters with 1-atom Depend,  Nr of types of atoms
		d1var rat,rapt,vnq;
		d2var r0s,r0p,r0pp; 

		std::vector<std::string> atmname; 
		d1var Val,Valboc, mass, bo131,bo132,bo133, Vale,plp1,nlpopt,plp2; 
		d1var povun2,povun3,povun4, povun5,povun6,povun7,povun8; 

// Valency Terms (j-dependancy only):
		d1var pval3,pval5, Valangle,Valval; 

// Van der Waals Terms:
		d1var rvdw1, eps, alf, vop; 
		d2var Dij, rvdW, alpij, gamW; 

// Coulomb & Charge equilibration:
		d1var chi, eta, gam; 
		d2var gamij; 

// 2-atom dependent parameters:
		d1var pbo1,pbo2,pbo3,pbo4,pbo5,pbo6,bom;
		d1var pboc1,pboc2,pboc3,pboc4,pboc5;
		d1var povun1,ovc, v13cor; 

// Bond Energy parameters (eq. 6)
		d1var Desig,Depi,Depipi, pbe1,pbe2;

// 3-atom dependent parameters
		d1var pval1,pval2,pval4,pval6,pval7,pval8,pval9,pval10;
		d1var ppen1,ppen2,ppen3,ppen4,pcoa1,pcoa2,pcoa3,pcoa4;
		d1var theta00;

// 4-atom dependent parameters
		d1var ptor1,ptor2,ptor3,ptor4, V1,V2,V3, pcot1,pcot2;

// hydrogen-bond params
		d1var r0hb,phb1,phb2,phb3;
};
