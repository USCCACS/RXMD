#include <iostream>

using namespace std; 

extern "C" void CppInterface(void *param_ptr, void *pos_ptr, void *atype_ptr, void *nbrlist_ptr, void *bo_ptr); 

void CppInterface(void *param_ptr, void *pos_ptr, void *atype_ptr, void *nbrlist_ptr, void *bo_ptr) 
{
	int *param = (int *)param_ptr;
	double *pos = (double *)pos_ptr;
	double *atype = (double *)atype_ptr;
	int *nbrlist = (int *)nbrlist_ptr;
	double *bo= (double *)bo_ptr;

	int myid=param[0]; 
	int NATOMS=param[1]; 
	int NBUFFER_P=param[2]; 
	int NBUFFER_N=param[3]; 
	int MAXNEIGHBS=param[4]; 
	int STRIDE=NBUFFER_P+NBUFFER_N+1;

    for (int i=0; i<10; ++i)
	{
		int NumNeighbors=nbrlist[i];
		cout << "myid = " << myid << " i = " << i << " nnbr = " << NumNeighbors << ": " ; 
		for (int j = 1; j <= NumNeighbors; ++j) 
		{
			int nbr = nbrlist[i+j*STRIDE];
			cout << " " << nbr <<
			" (" << pos[3*nbr] << "," << pos[3*nbr+1] << "," << pos[3*nbr+2] << ")," ; 
		}
		cout << endl; 

/*
		cout << "myid = " << myid << " i = " << i << " i1 = " << i1 << " atype = " << atype[i1] << 
		" (x,y,z) = (" << pos[i1*3-2] << ", " << pos[i1*3-1] << ", " << pos[i1*3] << ")" << endl; 
*/
	}

	return; 
}
