#include <iostream>
#include "Basic.h"

using namespace std; 

extern "C" void CppInterface(void *param_ptr, void *pos_ptr, void *atype_ptr, void *nbrlist_ptr, void *bo_ptr); 

void CppInterface(void *param_ptr, void *pos_ptr, void *atype_ptr, void *nbrlist_ptr, void *bo_ptr) 
{
	int *param = (int *)param_ptr;
	double *pos = (double *)pos_ptr;
	double *atype = (double *)atype_ptr;
	int *nbrlist = (int *)nbrlist_ptr;
	double *bo= (double *)bo_ptr;

	Basic b(param, pos, atype, nbrlist);

	cout << b << endl; 
	//b.PrintNeighbors();

	return; 
}
