# HTheta
Particle-in-cell code for the simulation of the theta-angle dynamics inside a Hall thruster

# Compilation and running instructions (not on the CINECA cluster)
1) Compile the latest code version, by typing in the /src folder the following command:
make clean –f makefile.txt
make HTHETA –f makefile.txt
2) Copy the simulation input file «sim_name.inp» (this can be prepared from the template provided in /inp/readme.txt) as the standard input file «params.inp», located in the /inp folder:
cp /inp/sim_name.inp /inp/params.inp
3) Run the executable located in the /bin folder:
/bin/htheta.exe &
4) Resuls of the simulations are stored in the /out folder, with names «sim_name_history.out» and «sim_name_parameters.out»

# Compilation and running instructions (on the CINECA cluster on a single node)
