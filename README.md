# HTheta
Particle-in-cell code for the simulation of the theta-angle dynamics inside a Hall thruster

# Compilation and running instructions (not on the CINECA cluster)
1) Create the environmental variable HTHETA pointing to the location of the repository (this can be added to .bashrc to avoid having to do it every time the shell is entered)
export HTHETA=local_path_2_main_folder_repository  
2) Compile the latest code version, by typing in the /src folder the following command:
make clean –f makefile.txt
make HTHETA –f makefile.txt
3) Copy the simulation input file «sim_name.inp» (this can be prepared from the template provided in /inp/readme.txt) as the standard input file «params.inp», located in the /inp folder:
cp /inp/sim_name.inp /inp/params.inp
4) Run the executable located in the /bin folder:
/bin/htheta.exe &
5) Resuls of the simulations are stored in the /out folder, with names «sim_name_history.out» and «sim_name_parameters.out»

# Compilation and running instructions (on the CINECA cluster on a single node)
