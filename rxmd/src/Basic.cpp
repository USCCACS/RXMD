#include <iomanip>
#include "Basic.h"

using namespace std; 

int Basic::PrintNeighbors()
{
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
	}

	return 0; 
}

std::ostream& operator<<(std::ostream& os, const Basic& obj)
{
    const int myid=obj.get_param(0);
    const int NATOMS=obj.get_param(1);
    const int NBUFFER_P=obj.get_param(2);
    const int NBUFFER_N=obj.get_param(3);
    const int MAXNEIGHBS=obj.get_param(4);
    const int STRIDE=NBUFFER_P+NBUFFER_N+1;

    for (int i=0; i<10; ++i)
    {
        int NumNeighbors=obj.get_nbrlist(i);
        os << "myid = " << myid << " i = " << i << " nnbr = " << NumNeighbors << ": " ;
        for (int j = 1; j <= NumNeighbors; ++j)
        {
            os << std::setprecision(4);
            int nbr = obj.get_nbrlist(i+j*STRIDE);
            os << " " << nbr <<
            " (" << obj.get_pos(3*nbr) << "," << obj.get_pos(3*nbr+1) << "," << obj.get_pos(3*nbr+2) << ")," ;
            os << std::setprecision(7);
        }
        os << std::endl;
    }

    return os;
};
