extern "C" void init_rxmdtorch(int myrank);
extern "C" void get_nn_force_torch(int nlocal, int ntotal, int nbuffer, void *pos_ptr, void *type_ptr, void *force_ptr, void *nbr_ptr, double &energy);
extern "C" void get_maxrc_rxmdnn(double & maxrc);
