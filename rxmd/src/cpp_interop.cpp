#include <iostream>

using namespace std; 

extern "C" void cpp_interop(void *param_ptr, void *pos_ptr, void *atype_ptr); 

void cpp_interop(void *param_ptr, void *pos_ptr, void *atype_ptr)
{
	int *param = (int *)param_ptr;
	double *pos = (double *)pos_ptr;
	double *atype = (double *)atype_ptr;

	int NATOMS=param[0]; 
	int NBUFFER_P=param[1]; 
	int NBUFFER_N=param[2]; 

    for (int i=0; i<NATOMS; i++)
	{
		cout << "i = " << i << " atype = " << atype[i] << endl; 
	}

	return; 
}
