#include <iostream>
#include <iomanip>

class Basic
{
	public:
		Basic(int *param_ptr, double *pos_ptr, double *atom_type_ptr, int *nbrlist_ptr) : 
			param(param_ptr), pos(pos_ptr), atom_type(atom_type_ptr), nbrlist(nbrlist_ptr) {}; 

		int PrintNeighbors();

		const int& get_param(int i) const { return param[i]; };
		const double& get_pos(int i) const { return pos[i]; };
		const int& get_nbrlist(int i) const { return nbrlist[i]; };

	private:
		int *param;
		double *pos;
		double *atom_type;
		int *nbrlist;
};

std::ostream& operator<<(std::ostream& os, const Basic& obj); 
