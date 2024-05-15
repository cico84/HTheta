# HTheta
Particle-in-cell code for the simulation of the theta-angle dynamics inside a Hall thruster

Before trying to compile and run the code, create the environmental variable HTHETA pointing to the location of the repository (this can be added to .bashrc to avoid having to do it every time the shell is entered) by typing:
export HTHETA=local_path_2_main_folder_repository  

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
1) Copy the simulation input file «sim_name.inp» (this can be prepared from the template provided in /inp/readme.txt) as the standard input file «params.inp», located in the /inp folder:
cp /inp/sim_name.inp /inp/params.inp
2) Go to the /bin folder and execute the program by typing:
./run_code_g100.sh
This will compile the code and run it by loading the required modules (for both compilation on front end and execution on compute nodes)
3) Resuls of the simulations are stored in the /out folder, with names «sim_name_history.out» and «sim_name_parameters.out»