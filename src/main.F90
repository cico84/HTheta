!******************************************************************************          
!*                                  HTHETA                                    *
!*                1D{theta}-3V FULLY KINETIC PIC SIMULATION                   *
!*                    OF ACCELERATION REGION IN SPT-100                       *               
!*                       (LOCAL FIELD APPROXIMATION)                          *
!*                                                                            *
!******************************************************************************

!  GOAL

!  Simulate the dynamics of Hall thruster along the azimuthal direction
!  Study of the anomalous axial drift

!  RECORD OF THE REVISIONS

!    Date          Programmer            Description of the changes
!  ========    ==================    ========================================================
!  01/03/18    Francesco Taccogna    Asadi-Bogopolsky first version
!  01/05/24    Filippo Cichocki      Parallelization with Open MP
           
!  CHARACTERISTICS OF MODEL:

!  1. 1D Cartesian model in the azimuthal direction: theta -> y
!  2. Local field approximation (densXe,Ez,Br -> fixed value)
!  3. Each macro-particle represents a charge per unit surface
!  4. Uniform spatial grid
!  5. PIC for electrons and ions ("linear weighting" function, CIC)
!  6. Electrostatic Poisson's equation solution
!  7. Electric field with centered 2-points gradient scheme
!  8. Leap-frog method for particle push
!  9. No collisions

! NOMENCLATURE:
! ny             = N. of cells along the azimuthal direction (along z there is only one cell) 
! y(j)           = Grid point azimuthal coordinate
! ype (i)        = Azimuthal coordinate of the macro-particle "i"
! vxpe(i)        = x-component (radial) of the macro-particle velocity
! vype(i)        = y-component (azimuthal) of the macro-particle velocity
! vzpe(i)        = z-component (axial) of the macro-particle velocity
! rhoi(j)        = Ion charge density at the grid node "j"
! rhoe(j)        = Electron charge density at the grid node "j"
! phi (j)        = Electric potential at grid node "j"
! apoi,bpoi,cpoi -> Poisson's equation discretization coefficients in Cartesian topology  
! Ey  (j)        = Azimuthal electric field at grid point "j"

!******************************************************************************
!*                             DEFINITION     PHASE                           *
!******************************************************************************

      ! Modules
      
      ! ! NVTX profiling module (unused in Open MP)
      ! module nvtx

      !       use iso_c_binding
      !       implicit none
            
      !       integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00', Z'00ff00ff', Z'0000ffff', Z'00ff0000', Z'00ffffff']
      !       character(len=256),private :: tempName
            
      !       type, bind(C):: nvtxEventAttributes
      !             integer(C_INT16_T):: version=1
      !             integer(C_INT16_T):: size=48 !
      !             integer(C_INT):: category=0
      !             integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
      !             integer(C_INT):: color
      !             integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
      !             integer(C_INT):: reserved0
      !             integer(C_INT64_T):: payload   ! union uint,int,double
      !             integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
      !             type(C_PTR):: message  ! ascii char
      !       end type
            
      !       interface nvtxRangePush
      !             ! push range with custom label and standard color
      !             subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
      !                   use iso_c_binding
      !                   character(kind=C_CHAR,len=*) :: name
      !             end subroutine
                  
      !             ! push range with custom label and custom color
      !             subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
      !                   use iso_c_binding
      !                   import:: nvtxEventAttributes
      !                   type(nvtxEventAttributes):: event
      !             end subroutine
      !       end interface
            
      !       interface nvtxRangePop
      !             subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
      !             end subroutine
      !       end interface
            
      !       contains
            
      !             subroutine nvtxStartRange(name,id)
      !                   character(kind=c_char,len=*) :: name
      !                   integer, optional:: id
      !                   type(nvtxEventAttributes):: event
                        
      !                   tempName=trim(name)//c_null_char
                        
      !                   if ( .not. present(id)) then
      !                         call nvtxRangePush(tempName)
      !                   else
      !                   event%color=col(mod(id,7)+1)
      !                   event%message=c_loc(tempName)
      !                   call nvtxRangePushEx(event)
      !                   end if
      !             end subroutine
                  
      !             subroutine nvtxEndRange
      !                   call nvtxRangePop
      !             end subroutine
      
      ! end module nvtx

      module const_phys
            ! Physical and mathematical constants
            implicit none
            save
            double precision, parameter :: pi	  = acos(-1.d0)    ! Pi Greek
            double precision, parameter :: eps0	  = 8.854188d-12   ! Vacuum Dielectric Constant
            double precision, parameter :: q	  = 1.602189d-19   ! Electric charge [C]
            double precision, parameter :: me	  = 9.1093897d-31  ! Electron mass [Kg]
            double precision, parameter :: Mi	  = 2.1801714d-25	 ! Xenon atom Mass (Xe) [Kg]
            double precision, parameter :: kB	  = 1.380662d-23   ! Boltzmann Constant [J/K]
            double precision, parameter :: JtoeV  = 6.24146d18 	 ! Conversion factor from J to eV
            double precision, parameter :: eVtoK  = 11600          ! Conversion factor from eV to K
            double precision, parameter :: duepi  = 2.*pi
            double precision            :: w                       ! Macroparticle's statistical weight (INPUT)
            double precision            :: dt                      ! Timestep [s] (INPUT)
            double precision            :: conste
            double precision            :: consti
            double precision            :: wq
      end module const_phys

      module system
            ! Physical external parameters
            implicit none
            save
            double precision :: yd                  ! Azimuthal domain length [m] (INPUT)
            double precision :: zch                 ! Channel length [m] (INPUT)
            double precision :: zacc                ! Axial domain length [m] (INPUT)
            double precision :: Ez0                 ! Axial electric field [V/m] (INPUT)
            double precision :: Br0                 ! Radial magnetic field [T] (INPUT)
            double precision :: n0                  ! Average plasma density [/m3] (INPUT)
            double precision :: Te0                 ! Initial electron temperature [K] (INPUT)
            double precision :: Ti0                 ! Initial ion temperature [k] (INPUT)
      end module system

      module grid
            use const_phys
            use system
            implicit none
            save
            integer                                     :: ny         ! Number of cells along azimuthal direction (INPUT)
            double precision                            :: dy, duedy
            ! Datasets from 0 to ny
            double precision, allocatable, dimension(:) :: y, vol
      end module grid

      module poi
            use grid
            implicit none
            save
            double precision                            :: apoi, bpoi, cpoi
            ! Datasets from 0 to ny
            double precision, allocatable, dimension(:) :: rhoe,rhoi,phi,Ey,dpoi
      end module poi

      module part
            use grid
            implicit none
            save
            integer            :: npmax,npic
            integer 		 :: npe,npi
            ! Datasets from 1 to npmax
            double precision, allocatable, dimension(:) :: ype, zpe, vxpe, vype, vzpe
            double precision, allocatable, dimension(:) :: ypi, zpi, vxpi, vypi, vzpi
            double precision                            :: vxpeprox, vxpiprox
      end module part

      module diagn
            implicit none
            save
            character(len=100) :: out_name  ! Name of output files related to this simulation (INPUT)
            character(len=100) :: out_path
            double precision   :: mob, Ee_ave, Ei_ave, Debye, omegae, omegace, vthetae, vzi, vthe, CFLe1, CFLe2
            double precision   :: time_ini, time_push, time_poisson, time_scatter, time_othr, time_pic_cycle

            integer(8)         :: nticks_sec
            integer(8)         :: nticks_max
            integer(8)         :: t00, t0, t1, t2, t3, t4, t5

            integer            :: npinit, nthreads
      end module diagn

      module rand
            implicit none
            save
            double precision                     :: rs, rs1, rs2
            integer*8			             :: iseed
            !$omp threadprivate(iseed)
            integer*8, dimension(:), allocatable :: seedNum ! Random seeds of the parallel threads
      end module rand

      program HTheta
            
            use poi
            use part
            use diagn
            use rand
            use const_phys
            use system
            use grid 
            use poi
            use part
            use omp_lib
            !use nvtx

            implicit none
            integer          :: ipic, j, i 
            integer*8        :: ini_seed

      !******************************************************************************
      !******************************************************************************


      !******************************************************************************
      !*                                                                            *
      !*                       INITIALIZATION       PHASE                           *
      !*                                                                            *  
      !******************************************************************************

            call system_clock(COUNT_RATE = nticks_sec, COUNT_MAX = nticks_max)

            call system_clock(t00)

            ! Initialize computation times
            time_ini     = 0.0d0
            time_push    = 0.0d0 
            time_poisson = 0.0d0 
            time_scatter = 0.0d0 
            time_othr    = 0.0d0

            !call nvtxStartRange('FULL PROGRAM')

            !call nvtxStartRange('INPUT READING')
            
            ! Read input parameters
            call read_input_parameters(out_name, out_path, nthreads)

            ! Allocate data vectors
            call allocate_fields()

            ! Set number of active Open MP threads
            call omp_set_num_threads(nthreads)

            write(*,*) " Number of requested Open MP threads:", nthreads

            ! Plasma parameters estimations
            ! Debye length
            Debye=dsqrt(eps0*kB*Te0/(n0*q*q))   
            ! Plasma/cyclotron frequency
            omegae=dsqrt(n0*q*q/(me*eps0))
            omegace=q*Br0/me
            ! Azimuthal drift velocity of electrons
            vthetae=Ez0/Br0
            ! Axial ion velocity increment
            vzi = dsqrt(2.*q*Ez0*zacc/Mi)
            ! Number of initial electrons/ions (every particle is a charge per unit surface C/m^2)      
            npinit = nint(yd*n0/w)         

            !call nvtxEndRange

            !******************************************************************************

            !call nvtxStartRange('MESH GENERATION')
            ! Mesh generation
            call mesh
            !call nvtxEndRange

            ! Constraint CFL      
            vthe  = dsqrt(8.*kB*Te0/(pi*me))
            CFLe1 = vthe*dt/dy
            CFLe2 = vthetae*dt/dy     
            
            open  (3, file = trim(out_path)//'/'//trim(out_name)//'_parameters.out',status='unknown')
            write (3,*) 'lD,dy=       ',Debye,dy
            write (3,*) 'dtp,dtB,dt=  ',1./omegae,1./omegace,dt
            write (3,*) 'CFLe=  	 ',CFLe1,CFLe2
            write (3,*) 'vthe,vde,vzi=',vthe,vthetae,vzi
            write (3,*) '# part,ppc=  ',npinit,float(npinit)/float(ny-1)
            write (3,*) 'Total particle size [MB] treated by the PUSH function', &
                        1d-6*(10d0 * float(npinit) * 8d0 + 2d0 * float(ny+1) * 8d0 )
            close (3)      
            
            !******************************************************************************

            ! Poisson's equation coefficients: cartesian coordinate {y}
            ! apoi*phi(j-1)+bpoi*phi(j)+cpoi*phi(j+1)=dpoi(j)=-(rhoi(j)-rhoe(j))*dy**2/eps0
            ! PERIODICITY CONDITIONS ALONG Y: phi(0)=phi(ny)

            apoi=1.
            cpoi=apoi
            bpoi=-2.

            ! Reference potential at y = y(0) = y(ny)
            phi(ny) = 0.
            phi(0 ) = 0.

            !******************************************************************************   

            ! Initialization of random numbers generator
            !iseed = -1
            !gen   = ran2(iseed)
            ini_seed = 1000000
            allocate( seedNum(nthreads) )
            call ini_seeds (nthreads, ini_seed, seedNum )

            !******************************************************************************

            ! Initial particle distribution
            call init
            
            call system_clock(t1)
            time_ini = dble(t1 - t00) / dble(nticks_sec)

            !******************************************************************************      
            open (11,file=trim(out_path)//'/'//trim(out_name)//'_history.out',status='unknown',position='append')
            open (12,file=trim(out_path)//'/'//trim(out_name)//'_cpu_times.out',status='unknown',position='append')

            write(12,'(1A25,1E10.3)') " Initialization time [s]:", time_ini
            write(12,'(1A60)') "  SCATTER [s]    POISSON [s]     PUSH [s]        OTHR [s]   "
            ! *****************************************************************************
            ! ******************************** PIC cycle **********************************
            ! *****************************************************************************
            do ipic = 1, npic

                  call system_clock(t0)
                  
                  ! Initialize computation times within the PIC cycle
                  time_scatter = 0.0d0 
                  time_poisson = 0.0d0 
                  time_push    = 0.0d0 
                  time_othr    = 0.0d0

                  !call nvtxStartRange('SCATTER PHASE')
                  ! Weight particles to the nodes of the mesh to obtain the charge density
                  call scatter
                  !call nvtxEndRange

                  call system_clock(t1)
                  time_scatter = dble(t1 - t0) / dble(nticks_sec)

                  ! Compute Poisson's equation source term
                  do j = 0, ny
                        dpoi(j) = -(rhoi(j)-rhoe(j)) * dy**2 / eps0
                  end do
            
                  !call nvtxStartRange('FIELD SOLVE')
                  ! Solve for the self-consistent azimuthal electric field
                  call fieldsolve
                  !call nvtxEndRange
                  call system_clock(t2)
                  time_poisson = dble(t2 - t1) / dble(nticks_sec)

                  !call nvtxStartRange('PUSH')
                  ! Update macro-particles positions and velocities
                  call push
                  !call nvtxEndRange
                  call system_clock(t3)
                  time_push = dble(t3 - t2) / dble(nticks_sec)

                  ! Diagnostics: averaged energy and mobility
                  Ee_ave = 0.
                  Ei_ave = 0.
                  mob    = 0.

                  !$omp parallel default(none) firstprivate(npe, npi) private(i) shared(vxpe,vype,vzpe,vxpi,vypi,vzpi)  &
                  !$omp                        reduction(+: Ee_ave, Ei_ave, mob)
                  !$omp do
                  do i = 1, npe
                        Ee_ave = Ee_ave + (vxpe(i)**2+vype(i)**2+vzpe(i)**2)
                        mob    = mob + vzpe(i)
                  end do
                  !$omp end do
                  !$omp do
                  do i = 1, npi
                        Ei_ave = Ei_ave + (vxpi(i)**2+vypi(i)**2+vzpi(i)**2)
                  end do
                  !$omp end do
                  !$omp end parallel 

                  Ee_ave = Ee_ave*JtoeV*0.5*me/npe
                  Ei_ave = Ei_ave*JtoeV*0.5*Mi/npi
                  mob    = mob/(npe*Ez0)

                  call system_clock(t4)
                  time_othr = dble(t4 - t3) / dble(nticks_sec)

                  if ( mod(ipic,1) .eq. 0 ) then
                        write(11,101) ipic*dt, Ee_ave, Ei_ave, mob, phi(ny/2), Ey(ny/2), rhoe(ny/2)/q, rhoi(ny/2)/q
                        write(12,102) time_scatter, time_poisson, time_push, time_othr
                  end if

                  call system_clock(t5)
                  time_pic_cycle = time_pic_cycle + dble(t5 - t0) / dble(nticks_sec)

            ! end of PIC cycle
            end do

101         format (8(2x,1pg13.5e3))
102         format (4(2x,1e12.3,1x))
            !******************************************************************************
            !******************************************************************************

            call deallocate_fields()

            !call nvtxEndRange

            time_pic_cycle = time_pic_cycle / npic
            call system_clock(t5)
            write(12,'(1A30,1E10.3)') " Total computational time [s]:",  dble(t5 - t00) / dble(nticks_sec)
            write(12,'(1A30,1E10.3)') " Average PIC cycle time   [s]:", time_pic_cycle

            close (11)
            close (12)
                        
      end program

      !******************************************************************************
      !******************************************************************************

      !******************************************************************************        
      !******************************************************************************
      !*                                                                            *
      !*                              LIST OF SUBROUTINES                           *
      !*                                                                            *
      !******************************************************************************        
      !******************************************************************************

      !******************************************************************************
      !******************************************************************************
      subroutine read_input_parameters(out_name, out_path, nthreads)
            
            use const_phys
            use system
            use part
            use grid
            implicit none
            character(len=100), intent(out) :: out_name, out_path
            character(len=100)              :: home, inp_file
            integer,            intent(out) :: nthreads 
            integer                         :: istatus

            namelist /sim_settings/ out_name, nthreads, w, dt, yd, zch, zacc, Ez0, Br0, &
                                    n0, Te0, Ti0, npmax, npic, ny

            ! Read the standard input file to read the simulation name (ID 1)
            call getenv('HTHETA', home)
            inp_file = trim(home)//'/inp/params.inp'
            out_path = trim(home)//'/out/'

            ! Read the input parameters file (ID 10)
            open ( unit = 10, file = inp_file, iostat = istatus, form = 'formatted',  &
                  access = 'sequential', status = 'old' )
            if (istatus > 0) then
                  write(*,*) '  ERROR: Input parameters file is invalid'
                  write(*,*) " Input parameters file:", inp_file
                  stop
            end if

            read (10, sim_settings)
            conste = q*dt/(2.*me)
            consti = q*dt/Mi
            wq     = w*q 

            write(*,*) " Read parameters:"
            write(*,*) " Output file name:", out_name
            write(*,*) " Macro-particle statistical weight:", w
            write(*,*) " PIC time step [s]:", dt
            write(*,*) " Extension along the azimuthal direction [m]:", yd
            write(*,*) " Extension of channel [m]:", zch
            write(*,*) " Acceleration region extension [m]:", zacc
            write(*,*) " Axial electric field [V/m]:", Ez0
            write(*,*) " Radial magnetic field [T]:", Br0
            write(*,*) " Plasma density [/m3]:", n0
            write(*,*) " Electron and ion temperatures [K]:", Te0, Ti0
            write(*,*) " Max. number of macro-particles per species:", npmax
            write(*,*) " Number of PIC time steps:", npic
            write(*,*) " Number of cells along the azimuthal direction:", ny
            write(*,*) " HTHETA: running simulation..."

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************
      subroutine allocate_fields()
            
            ! This subroutine allocates the arrays required along the azimuthal direction
            ! and those related to macro-particles

            use const_phys
            use system
            use grid
            use part
            use poi
            implicit none

            ! Datasets from 0 to ny
            allocate( y   (0:ny) )
            allocate( vol (0:ny) )
            
            allocate( rhoe(0:ny) )
            allocate( rhoi(0:ny) )
            allocate( phi (0:ny) )
            allocate( Ey  (0:ny) )
            allocate( dpoi(0:ny) )
            
            ! Datasets from 1 to npmax
            allocate(  ype(1:npmax) ) 
            allocate(  zpe(1:npmax) ) 
            allocate( vxpe(1:npmax) ) 
            allocate( vype(1:npmax) ) 
            allocate( vzpe(1:npmax) ) 
            allocate(  ypi(1:npmax) ) 
            allocate(  zpi(1:npmax) ) 
            allocate( vxpi(1:npmax) ) 
            allocate( vypi(1:npmax) ) 
            allocate( vzpi(1:npmax) )

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************

      subroutine deallocate_fields()
            
            ! This subroutine allocates the arrays required along the azimuthal direction
            ! and those related to macro-particles

            use grid
            use part
            use poi
            implicit none

            ! Datasets from 0 to ny
            deallocate( y, vol, rhoe, rhoi, phi, Ey, dpoi )
            
            ! Datasets from 1 to npmax
            deallocate(  ype, zpe, vxpe, vype, vzpe, ypi, zpi, vxpi, vypi, vzpi )

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************
            
      subroutine mesh 
          
            use grid
            implicit none
            integer :: j
      
            y(0)  =  0.
            dy    = yd/ny
            duedy = 2.*dy
            do j = 1, ny
                  y(j) = y(j-1) + dy
            end do
            
            ! Node volumes
            do j = 0, ny
                  vol(j) = dy
            end do

            return

      end subroutine
                     
      !******************************************************************************
      !******************************************************************************

      subroutine init 
          
            use grid
            use part
            use diagn
            use rand
            implicit none
            integer                    :: i
            double precision           :: duekteme, duektimi
            double precision, external :: ran2
      
            npe = 0
            npi = 0

            iseed    = seedNum(1)
      
            duekteme = -2.*kB*Te0 / me
            duektimi = -2.*kB*Ti0 / Mi
            
            ! Uniform distribution of initial electron and ion macro-particles
            do i=1,npinit
            
                  ! Electrons macro-particles:
                  npe=npe+1
                  rs=ran2(iseed)
                  ype(npe)=rs*y(ny)
                  rs=ran2(iseed)
                  zpe(npe)=(zch-zacc)+rs*zacc

                  ! Full Maxwellian distribution (Box-Muller transformation)         
                  rs1=ran2(iseed)
                  rs2=ran2(iseed)
                  if ((MOD(i+1,2)).eq.(MOD(2,2))) then         
                  vxpe(npe)=DSQRT(duekteme*LOG(rs1))*DCOS(duepi*rs2)
                  vxpeprox=DSQRT(duekteme*LOG(rs1))*DSIN(duepi*rs2)
                  else
                  vxpe(npe)=vxpeprox
                  end if
                  rs1=ran2(iseed)
                  rs2=ran2(iseed)
                  vype(npe)=DSQRT(duekteme*LOG(rs1))*DCOS(duepi*rs2) 
                  vzpe(npe)=DSQRT(duekteme*LOG(rs1))*DSIN(duepi*rs2)
                  
                  ! Ion macro-particles:
                  npi=npi+1
                  ypi(npi)=ype(npe)
                  rs=ran2(iseed)
                  zpi(npi)=(zch-zacc)+rs*zacc

                  ! Full Maxwellian distribution (Box-Muller transformation): vzi        
                  rs1=ran2(iseed)
                  rs2=ran2(iseed) 
                  if ((MOD(i+1,2)).eq.(MOD(2,2))) then         
                  vxpi(npi)=DSQRT(duektimi*LOG(rs1))*DCOS(duepi*rs2)
                  vxpiprox=DSQRT(duektimi*LOG(rs1))*DSIN(duepi*rs2)
                  else
                  vxpi(npi)=vxpiprox
                  end if
                  rs1=ran2(iseed)
                  rs2=ran2(iseed)
                  vypi(npi)=DSQRT(duektimi*LOG(rs1))*DCOS(duepi*rs2)   
                  vzpi(npi)=DSQRT(duektimi*LOG(rs1))*DSIN(duepi*rs2)

            end do

      return
      end subroutine

      !******************************************************************************
      !******************************************************************************                              

      subroutine scatter 
     
            use poi
            use part
            implicit none
            integer          :: i,j, jp
            double precision :: wy
      
            ! Charge density initialization to 0
            do j = 0, ny
                  rhoe(j)=0.
                  rhoi(j)=0.
            end do
            
            !$omp parallel default(none) firstprivate(wq) private(i, j, wy, jp) shared(ny, vol, npe, npi,     &
            !$omp                        y, dy, ype, ypi) reduction(+: rhoe, rhoi)  
            !$omp do 
            ! Electron charge deposition on the mesh points 
            do i = 1, npe     
                  ! Charge density weighting (linear weighting, CIC)
                  jp         =int(ype(i) / dy) + 1
                  wy         =   (y(jp)-ype(i)) / dy     
                  rhoe(jp-1) =    wy*wq + rhoe(jp-1)
                  rhoe(jp  ) =   (1d0-wy)*wq + rhoe(jp)
            end do
            !$omp end do
            !$omp do 
            ! Ion charge deposition on the mesh points 
            do i = 1, npi       
                  ! Charge density weighting (linear weighting, CIC) 
                  jp         =int(ypi(i) / dy) + 1
                  wy         =   (y(jp)-ypi(i)) / dy     
                  rhoi(jp-1) =    wy*wq + rhoi(jp-1)
                  rhoi(jp  ) =   (1d0-wy)*wq + rhoi(jp)
            end do
            !$omp end do
            !$omp end parallel

            ! Periodic boundary conditions for both ion and electron charge density
            rhoe(0)  = rhoe(0) + rhoe(ny)   
            rhoe(ny) = rhoe(0)
            rhoi(0)  = rhoi(0) + rhoi(ny)   
            rhoi(ny) = rhoi(0)      
            do j = 0, ny
                  rhoe(j)=rhoe(j)/vol(j)
                  rhoi(j)=rhoi(j)/vol(j)
            end do

            return

      end subroutine
            
       
      !******************************************************************************
      !******************************************************************************

      subroutine fieldsolve 
      
            use poi
            use part
            implicit none
            integer          :: j
            double precision :: rnap(ny-1),snap(ny-1),tnap(ny-1),den

            !******************************************************************************       
      
            ! Poisson's equation in cartesian coordinates {theta}->{y}
            ! apoi*phi(j-1) + bpoi*phi(j) + cpoi*phi(j+1) = dpoi(j)
            ! Thomas tridiagonal algorithm for periodic boundary conditions: phi(ny) = phi(0) = 0.
            ! [M. Napolitano, Comm. Appl. Num. Meth. 1, 11-15 (1985)]
            rnap(1)=-cpoi/bpoi
            snap(1)=-apoi/bpoi
            tnap(1)=dpoi(1)/bpoi
            do j=2,ny-1
                  den=(apoi*rnap(j-1)+bpoi)
                  rnap(j)=-cpoi/den
                  snap(j)=-apoi*snap(j-1)/den
                  tnap(j)=(dpoi(j)-apoi*tnap(j-1))/den
            end do
            do j=ny-1,1,-1
                  phi(j)=rnap(j)*phi(j+1)+tnap(j)
            end do

            !******************************************************************************

            ! Electric field equation solution
            do j=1,ny-1
                  Ey(j)=-(phi(j+1)-phi(j-1))/duedy
            end do
            ! Periodic boundary conditions
            Ey(0) = -(phi(1)-phi(ny-1))/duedy
            Ey(ny)= Ey(0)
            
            !******************************************************************************
      
            return
      end subroutine
      

      !******************************************************************************
      !******************************************************************************

      subroutine push
          
            use poi
            use part
            use rand
            use omp_lib
            implicit none
            integer                    :: i,ie,ii,thread_num,jp
            double precision           :: tB,tBB,duekteme,duektimi,vmod,ang
            double precision           :: vyea,vyeb,vyec,vzea,vzeb,vzec,wy,Eyp
            double precision           :: y_min, y_max
            double precision, external :: ran2

            ! Leapfrog method (Boris algorithm)
            tB  = conste*Br0
            tBB = 2.*tB/(1.+tB**2)
            duekteme = -2.*kB*Te0/me   ! Electrons leap frog constant
            duektimi = -2.*kB*Ti0/Mi   ! Ions leap frog constant
            y_min    = y( 0)
            y_max    = y(ny)

            !$omp parallel default(none) firstprivate(tB, tBB, duekteme, duektimi, conste, consti, Ez0, dt, zch, zacc, &
            !$omp                                     y_min, y_max)                                                    &
            !$omp                             private(thread_num, ie, ii, i, jp, wy, Eyp, vyea, vzea, vyeb, vzeb,      &
            !$omp                                     vyec, vzec, rs1, rs2, ang, vmod, vxpeprox, vxpiprox)             &
            !$omp                              shared(ny, y, Ey, npe, npi, vxpe, vype, vzpe, vxpi, vypi, vzpi, ype,    &
            !$omp                                     zpe, ypi, zpi, seedNum, dy)
            thread_num  = omp_get_thread_num()
            iseed       = seedNum (thread_num+1)
            ie          = 0
            ii          = 0

            ! Electrons loop
            !$omp do
            do i = 1, npe

                  ! Compute index of occupied cell, weighting factor and electric field at particle
                  jp         =int(ype(i) / dy) + 1
                  wy         =   (y(jp)-ype(i)) / dy     
                  Eyp        = wy*Ey(jp-1) + (1.-wy)*Ey(jp)
                  
                  ! First half acceleration by electric field
                  vyea    = vype(i) - conste*Eyp
                  vzea    = vzpe(i) - conste*Ez0
                  ! Full rotation around magnetic field
                  vyeb    = vyea - vzea*tB
                  vzeb    = vzea + vyea*tB 
                  vyec    = vyea - vzeb*tBB
                  vzec    = vzea + vyeb*tBB
                  ! Second half acceleration by electric field  
                  vype(i) = vyec - conste*Eyp               
                  vzpe(i) = vzec - conste*Ez0
                  ! Coordinates updates (cartesian approximation)
                  ype (i) = ype(i) + vype(i)*dt
                  zpe (i) = zpe(i) + vzpe(i)*dt
                  ! Periodic boundary conditions
                  if      ( ype(i) .le. y_min  ) then
                        ype(i) = y_max  + ype(i)
                  else if (ype(i)  .ge. y_max ) then
                        ype(i) = ype(i) -  y_max
                  end if
                  ! ! Refresh particles
                  ! if ((zch-zpe(i)).ge.zacc) then
                  !       ie=ie+1
                  !       zpe(i)=zch
                  !       ! Full-Maxwellian distribution (Box-Muller transformation)         
                  !       rs1=ran2(iseed)
                  !       vmod=DSQRT(duekteme*DLOG(rs1))
                  !       rs2=ran2(iseed)
                  !       ang=duepi*rs2
                  !       if ((MOD(ie+1,2)).eq.(MOD(2,2))) then         
                  !             vxpe(i)=vmod*DCOS(ang)
                  !             vxpeprox=vmod*DSIN(ang)
                  !       else
                  !             vxpe(i)=vxpeprox
                  !       end if
                  !       rs1=ran2(iseed)
                  !       vmod=DSQRT(duekteme*DLOG(rs1))
                  !       rs2=ran2(iseed)
                  !       ang=duepi*rs2
                  !       vype(i)=vmod*DCOS(ang) 
                  !       vzpe(i)=vmod*DSIN(ang)
                  ! end if
            end do
            !$omp end do

            !$omp do
            ! Ions loop   
            do i = 1, npi

                  ! Compute index of occupied cell, weighting factor and electric field at particle
                  jp         =int(ypi(i) / dy) + 1
                  wy         =   (y(jp)-ypi(i)) / dy     
                  Eyp        = wy*Ey(jp-1) + (1.-wy)*Ey(jp)

                  vypi(i)    = vypi(i) + consti*Eyp
                  vzpi(i)    = vzpi(i) + consti*Ez0
                  ypi(i)     = ypi (i) + vypi(i)*dt
                  zpi(i)     = zpi (i) + vzpi(i)*dt

                  ! Periodic boundary conditions
                  if      ( ypi(i) .lt. y_min ) then
                        ypi(i) = y_max  + ypi(i)
                  else if ( ypi(i) .ge. y_max ) then
                        ypi(i) = ypi(i) -  y_max
                  end if

                  ! ! Refresh particles
                  ! if (zpi(i).ge.zch) then
                  !       ii=ii+1
                  !       zpi(i)=zch-zacc
                  !       ! Full-Maxwellian distribution (Box-Muller transformation)         
                  !       rs1=ran2(iseed)
                  !       vmod=DSQRT(duektimi*DLOG(rs1))
                  !       rs2=ran2(iseed)
                  !       ang=duepi*rs2
                  !       if ((MOD(ii+1,2)).eq.(MOD(2,2))) then         
                  !             vxpi(i)=vmod*DCOS(ang)
                  !             vxpiprox=vmod*DSIN(ang)
                  !       else
                  !             vxpi(i)=vxpiprox
                  !       end if
                  !       rs1=ran2(iseed)
                  !       vmod=DSQRT(duektimi*DLOG(rs1))
                  !       rs2=ran2(iseed)
                  !       ang=duepi*rs2       
                  !       vypi(i)=vmod*DCOS(ang) 
                  !       vzpi(i)=vmod*DSIN(ang)
                  ! end if
            end do
            !$omp end do
            !$omp end parallel

            !******************************************************************************
            !******************************************************************************        

            return
      end subroutine

      !******************************************************************************
      !******************************************************************************

  
      !******************************************************************************
      !*                                                                            *
      !*                               LIST OF FUNCTIONS                            *
      !*                                                                            *
      !******************************************************************************


      !*****************************************************************************       
      function ran2(i)
  
            !Long period (> 2.E18) random number generator of l'Ecuyer with Bays-Durham shuffle and added safeguards (from Numerical Recipes).
            !Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
            !Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
            !RNMX should approximate the largest floating value that is less than 1.
            implicit none
            
            integer*8, parameter :: im1=2147483563, im2=2147483399, imm1=im1-1
            real*8, parameter    :: am=1./im1
            integer*8, parameter :: ia1=40014, ia2=40692, iq1=53668,iq2=52774, ir1=12211
            integer*8, parameter :: ir2=3791, ntab=32, ndiv=1+imm1/ntab
            real*8, parameter    :: epsr = 1.2e-7, rnmx = 1.-epsr
            
            integer*8 :: naux
            integer*8, save :: iv(ntab)=[(0, naux = 1, ntab)]
            integer*8, save :: iy=0, idum2 = 123456789      !!!data idum2/123456789/, iv/ntab*0/, iy/0/
            integer*8, save  :: idum
            !$omp threadprivate(iv, iy, idum2, idum)
            integer*8  :: j, k
            
            integer*8, intent(inout) :: i
            real*8 :: ran2
          
            idum=i
        
            if (idum.le.0) then              
              if (idum.eq.0) idum=-1          
              idum=-idum                
              idum2=idum
              do j=ntab+8,1,-1             
                k=idum/iq1
                if (idum.lt.0) idum=idum+im1
                if (j.le.ntab) iv(j)=idum
              end do
              iy=iv(1)
            end if
        
            k=idum/iq1                     
            idum=ia1*(idum-k*iq1)-k*ir1
        
            if (idum.lt.0) idum=idum+im1
        
            k=idum2/iq2
            idum2=ia2*(idum2-k*iq2)-k*ir2
        
            if (idum2.lt.0) idum2=idum2+im2 
        
            j=1+iy/ndiv 
            iy=iv(j)-idum2                               
            iv(j)=idum
                                 
            if (iy.lt.1) iy=iy+imm1
        
            i=idum
        
            ran2=min(am*iy,rnmx)          
        
            return
          end function
        
          subroutine ini_seeds (nthreads, seedNum0, seeds)
            ! This subroutine initializes the seeds of a ran2 pseudo random number generator
            ! ----------------------------------------------- INPUTS/OUTPUTS  ----------------------------------------------------
            integer*4, intent(in)                        :: nthreads    ! Number of parallel threads
            integer*8, intent(in)                        :: seedNum0    ! Global simulation seed number
            ! --------------------------------------------------- OUTPUTS --------------------------------------------------------
            integer*8, dimension(nthreads), intent(out)  :: seeds       ! Initial seeds for the parallel threads
            ! ---------------------------------------------- INTERNAL VARIABLES --------------------------------------------------
            integer*8                                    :: i, seed
            real*8                                       :: ran_num
            real*8, external                             :: ran2
            ! -------------------------------------- EXECUTABLE PART OF THE SUBROUTINE -------------------------------------------
            seed      = - seedNum0      ! Must be negative
            ! First thread uses the global simulation seed
            seeds(1)  = - seedNum0      ! Must be negative
            ! The other threads use a random generated seed number (negative)
            do i = 2, nthreads
              ran_num  = ran2( seed )
              seeds(i) = - int(ran_num * 2147483647.0, kind=8)   ! Must be negative
            end do
            write(*,*) " Initial seeds of the used threads:", seeds(1:nthreads)
        
          end subroutine ini_seeds

!*****************************************************************************
!*****************************************************************************
