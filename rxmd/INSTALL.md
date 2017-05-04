To build rxmd, make.inc must be updated.
Users could copy paste the config file of a specific machines from the config folder.
Or customize it for the target machine.

Choose the build target for make
  make <OPTIONS>

 OPTIONS:
   all
   init
   rxmd
   clean

parallel building is supported.
make all -j24
