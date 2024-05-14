This folder contains the input file, containing the input parameters definition

INPUT TEMPLATE
&sim_settings
out_name = 'sim01'   ! Output file name
nthreads = 20        ! Number of Open MP threads
w        = 5.d8      ! Macroparticle's statistical weight
dt       = 1.d-12    ! Timestep [s]
yd       = 0.05      ! Azimuthal domain length [m]
zch      = 0.025     ! Channel length [m]
zacc     = 0.015     ! Axial domain length [m]
Ez0      = 20000     ! Axial electric field [V/m]
Br0      = 0.018     ! Radial magnetic field [T]
n0       = 5.e17     ! Average plasma density [/m3]
Te0      = 116000    ! Initial electron temperature [K]
Ti0      = 2320      ! Initial ion temperature [k]
npmax    = 500000000 ! Maximum number of macro-particles in the domain
npic     = 10000000  ! Total number of PIC cycles
ny       = 5000      ! Number of cells along the azimuthal direction
/	