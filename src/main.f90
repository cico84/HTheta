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
! Eype(i)        = Azimuthal electric field at macro-particle position

!******************************************************************************
!*                             DEFINITION     PHASE                           *
!******************************************************************************

      ! Modules

      MODULE const_phys
            ! Physical and mathematical constants
            IMPLICIT NONE
            SAVE
            DOUBLE PRECISION, PARAMETER :: pi	  = acos(-1.d0)    ! Pi Greek
            DOUBLE PRECISION, PARAMETER :: eps0	  = 8.854188d-12   ! Vacuum Dielectric Constant
            DOUBLE PRECISION, PARAMETER :: q	  = 1.602189d-19   ! Electric charge [C]
            DOUBLE PRECISION, PARAMETER :: me	  = 9.1093897d-31  ! Electron mass [Kg]
            DOUBLE PRECISION, PARAMETER :: Mi	  = 2.1801714d-25	 ! Xenon atom Mass (Xe) [Kg]
            DOUBLE PRECISION, PARAMETER :: kB	  = 1.380662d-23   ! Boltzmann Constant [J/K]
            DOUBLE PRECISION, PARAMETER :: JtoeV  = 6.24146d18 	 ! Conversion factor from J to eV
            DOUBLE PRECISION, PARAMETER :: eVtoK  = 11600          ! Conversion factor from eV to K
            DOUBLE PRECISION, PARAMETER :: duepi  = 2.*pi
            DOUBLE PRECISION            :: w                       ! Macroparticle's statistical weight (INPUT)
            DOUBLE PRECISION            :: dt                      ! Timestep [s] (INPUT)
            DOUBLE PRECISION            :: conste
            DOUBLE PRECISION            :: consti
            DOUBLE PRECISION            :: wq
      END MODULE const_phys

      MODULE system
            ! Physical external parameters
            IMPLICIT NONE
            SAVE
            DOUBLE PRECISION :: yd                  ! Azimuthal domain length [m] (INPUT)
            DOUBLE PRECISION :: zch                 ! Channel length [m] (INPUT)
            DOUBLE PRECISION :: zacc                ! Axial domain length [m] (INPUT)
            DOUBLE PRECISION :: Ez0                 ! Axial electric field [V/m] (INPUT)
            DOUBLE PRECISION :: Br0                 ! Radial magnetic field [T] (INPUT)
            DOUBLE PRECISION :: n0                  ! Average plasma density [/m3] (INPUT)
            DOUBLE PRECISION :: Te0                 ! Initial electron temperature [K] (INPUT)
            DOUBLE PRECISION :: Ti0                 ! Initial ion temperature [k] (INPUT)
      END MODULE system

      MODULE grid
            USE const_phys
            USE system
            IMPLICIT NONE
            SAVE
            INTEGER                                     :: ny         ! Number of cells along azimuthal direction (INPUT)
            DOUBLE PRECISION                            :: dy, duedy
            ! Datasets from 0 to ny
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: y, vol
      END MODULE grid

      MODULE poi
            USE grid
            IMPLICIT NONE
            SAVE
            DOUBLE PRECISION                            :: apoi, bpoi, cpoi
            ! Datasets from 0 to ny
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rhoe,rhoi,phi,Ey,dpoi
      END MODULE poi

      MODULE part
            USE grid
            IMPLICIT NONE
            SAVE
            INTEGER            :: npmax,npic
            INTEGER 		 :: npe,npi
            ! Datasets from 1 to npmax
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ype, zpe, vxpe, vype, vzpe
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ypi, zpi, vxpi, vypi, vzpi
            DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: wye, Eype, wyi, Eypi
            DOUBLE PRECISION                            :: vxpeprox, vxpiprox
            INTEGER         , ALLOCATABLE, DIMENSION(:) :: jpe, jpi
      END MODULE part

      MODULE diagn
            IMPLICIT NONE
            SAVE
            CHARACTER(LEN=100) :: out_name  ! Name of output files related to this simulation (INPUT)
            CHARACTER(LEN=100) :: out_path
            DOUBLE PRECISION   :: mob, Ee_ave, Ei_ave, Debye, omegae, omegace, vthetae, vzi, vthe, CFLe1, CFLe2
            INTEGER            :: npinit
      END MODULE diagn

      MODULE rand
            IMPLICIT NONE
            SAVE
            INTEGER          :: iseed
            DOUBLE PRECISION :: rs, rs1, rs2
      END MODULE rand

      program HTheta
            
            USE poi
            USE part
            USE diagn
            USE rand
            USE const_phys
            USE system
            USE grid 
            USE poi
            USE part

      !******************************************************************************
      !******************************************************************************


      !******************************************************************************
      !*                                                                            *
      !*                       INITIALIZATION       PHASE                           *
      !*                                                                            *  
      !******************************************************************************

            ! Read input parameters
            call read_input_parameters(out_name, out_path)

            ! Allocate data vectors
            call allocate_fields()

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

            !******************************************************************************

            ! Mesh generation
            call mesh

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
            iseed = -1
            gen   = ran2(iseed)

            !******************************************************************************

            ! Initial particle distribution
            call init
            
            !******************************************************************************      

            ! *****************************************************************************
            ! ******************************** PIC cycle **********************************
            ! *****************************************************************************
            do ipic = 1, npic

                  ! Weight particles to the nodes of the mesh to obtain the charge density
                  call scatter

                  ! Compute Poisson's equation source term
                  do j = 0, ny
                        dpoi(j) = -(rhoi(j)-rhoe(j)) * dy**2 / eps0
                  end do
            
                  ! Solve for the self-consistent azimuthal electric field
                  call fieldsolve

                  ! Update macro-particles positions and velocities
                  call push

                  ! Diagnostics: averaged energy and mobility
                  Ee_ave = 0.
                  Ei_ave = 0.
                  mob    = 0.
                  do i = 1, npe
                        Ee_ave = Ee_ave + (vxpe(i)**2+vype(i)**2+vzpe(i)**2)
                        mob    = mob + vzpe(i)
                  end do
                  Ee_ave = Ee_ave*JtoeV*0.5*me/npe
                  mob    = mob/(npe*Ez0)
            
                  do i = 1, npi
                        Ei_ave = Ei_ave + (vxpi(i)**2+vypi(i)**2+vzpi(i)**2)
                  end do
                  Ei_ave = Ei_ave*JtoeV*0.5*Mi/npi

                  if ( mod(ipic,1) .eq. 0 ) then
                        open (11,file=trim(out_path)//'/'//trim(out_name)//'_history.out',status='unknown',position='append')
                        write(11,101) ipic*dt, Ee_ave, Ei_ave, mob, phi(ny/2), Ey(ny/2), rhoe(ny/2)/q, rhoi(ny/2)/q
                        close (11)
                  end if
            
      101  format (8(2x,1pg13.5e3))
            
            ! end of PIC cycle
            end do

            !******************************************************************************
            !******************************************************************************

            call deallocate_fields()
                        
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
      subroutine read_input_parameters(out_name, out_path)
            
            character(LEN=100), intent(out) :: out_name, out_path
            USE const_phys
            USE system

            namelist /sim_settings/ out_name, w, dt, yd, zch, zacc, Ez0, Br0, n0, &
                                    Te0, Ti0, npmax, npic, ny

            ! Read the standard input file to read the simulation name (ID 1)
            call getenv('HTHETA', home)
            inp_file = trim(home)//'/inp/params.inp'
            out_path = trim(home)//'/out/'

            ! Read the input parameters file (ID 10)
            open ( unit = 10, file = inp_file, iostat = istatus, form = 'formatted',  &
                  access = 'sequential', status = 'old' )
            if (istatus > 0) then
                  write(*,*) '  ERROR: Input parameters file is invalid'
                  write(*,*) " Input parameters file:", paths%inp_file_path
                  stop
            end if

            read (10, sim_settings)
            conste = q*dt/(2.*me)
            consti = q*dt/Mi
            wq     = w*q 

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************
      subroutine allocate_fields()
            
            ! This subroutine allocates the arrays required along the azimuthal direction
            ! and those related to macro-particles

            USE const_phys
            USE system
            USE grid
            USE part
            USE poi

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
            allocate(  wye(1:npmax) )
            allocate( Eype(1:npmax) )
            allocate(  wyi(1:npmax) )
            allocate( Eypi(1:npmax) )
            allocate(  jpe(1:npmax) )
            allocate(  jpi(1:npmax) )

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************

      subroutine deallocate_fields()
            
            ! This subroutine allocates the arrays required along the azimuthal direction
            ! and those related to macro-particles

            USE grid
            USE part
            USE poi

            ! Datasets from 0 to ny
            deallocate( y, vol rhoe, rhoi, phi, Ey, dpoi)
            
            ! Datasets from 1 to npmax
            deallocate(  ype, zpe, vxpe, vype, vzpe, ypi, zpi, vxpi, vypi, vzpi, &
                         wye, Eype, wyi, Eypi, jpe, jpi)

            return

      end subroutine

      !******************************************************************************
      !******************************************************************************
            
      subroutine mesh 
          
            USE grid
            IMPLICIT NONE
            INTEGER :: j
      
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
          
            USE grid
            USE part
            USE diagn
            USE rand
            IMPLICIT NONE
            INTEGER                    :: i
            DOUBLE PRECISION           :: duekteme,duektimi,vmod,ang
            DOUBLE PRECISION, EXTERNAL :: ran2
      
            npe = 0
            npi = 0
      
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
     
            USE poi
            USE part
            IMPLICIT NONE
            INTEGER          :: i,j
      
            ! Charge density initialization to 0
            do j = 0, ny
                  rhoe(j)=0.
                  rhoi(j)=0.
            end do

            ! Electron charge deposition on the mesh points 
            do i = 1, npe     
                  ! Charge density weighting (linear weighting, CIC)
                  jpe(i)=int(ype(i)/dy)+1         
                  wye(i)=(y(jpe(i))-ype(i))/dy     
                  rhoe(jpe(i)-1)=wye(i)*wq+rhoe(jpe(i)-1)
                  rhoe(jpe(i))=(1.-wye(i))*wq+rhoe(jpe(i))
            end do

            ! Periodic boundary conditions
            rhoe(0)  = rhoe(0) + rhoe(ny)   
            rhoe(ny) = rhoe(0)
            do j = 0, ny
                  rhoe(j)=rhoe(j)/vol(j)
            end do

            ! Ion charge deposition on the mesh points 
            do i=1,npi       
                  ! Charge density weighting (linear weighting, CIC) 
                  jpi(i)=int(ypi(i)/dy)+1         
                  wyi(i)=(y(jpi(i))-ypi(i))/dy     
                  rhoi(jpi(i)-1)=wyi(i)*wq+rhoi(jpi(i)-1)
                  rhoi(jpi(i))=(1.-wyi(i))*wq+rhoi(jpi(i))                   
            end do

            ! Periodic boundary condition 
            rhoi(0)=rhoi(0)+rhoi(ny)   
            rhoi(ny)=rhoi(0)      
            do j=0,ny
                  rhoi(j)=rhoi(j)/vol(j)
            end do
            
            return

      end subroutine
            
       
      !******************************************************************************
      !******************************************************************************

      subroutine fieldsolve 
      
            USE poi
            USE part
            IMPLICIT NONE
            INTEGER          :: j
            DOUBLE PRECISION :: rnap(ny-1),snap(ny-1),tnap(ny-1),den

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
          
            USE poi
            USE part
            USE rand
            IMPLICIT NONE
            INTEGER                    :: i,ie,ii
            DOUBLE PRECISION           :: tB,tBB,duekteme,duektimi,vmod,ang
            DOUBLE PRECISION           :: vyea,vyeb,vyec,vzea,vzeb,vzec
            DOUBLE PRECISION, EXTERNAL :: ran2

            ! Electrons      

            ! Leapfrog method (Boris algorithm)
            tB  = conste*Br0
            tBB = 2.*tB/(1.+tB**2)
            duekteme = -2.*kB*Te0/me

            ie = 0
            do i = 1, npe
                  Eype(i)=wye(i)*Ey(jpe(i)-1)+(1.-wye(i))*Ey(jpe(i))
                  ! first half acceleration by electric field
                  vyea=vype(i)-conste*Eype(i)         
                  vzea=vzpe(i)-conste*Ez0
                  ! Full rotation around magnetic field
                  vyeb=vyea-vzea*tB
                  vzeb=vzea+vyea*tB 
                  vyec=vyea-vzeb*tBB
                  vzec=vzea+vyeb*tBB
                  ! Second half acceleration by electric field  
                  vype(i)=vyec-conste*Eype(i)                
                  vzpe(i)=vzec-conste*Ez0
                  ! Coordinates updates (cartesian approximation)
                  ype(i)=ype(i)+vype(i)*dt
                  zpe(i)=zpe(i)+vzpe(i)*dt
                  ! Periodic boundary conditions
                  if (ype(i).le.y(0)) then
                  ype(i)=y(ny)+ype(i)
                  else if (ype(i).ge.y(ny)) then
                  ype(i)=ype(i)-y(ny)
                  end if
                  ! Refresh particles
                  if ((zch-zpe(i)).ge.zacc) then
                        ie=ie+1
                        zpe(i)=zch
                        ! Full-Maxwellian distribution (Box-Muller transformation)         
                        rs1=ran2(iseed)
                        vmod=DSQRT(duekteme*DLOG(rs1))
                        rs2=ran2(iseed)
                        ang=duepi*rs2
                        if ((MOD(ie+1,2)).eq.(MOD(2,2))) then         
                              vxpe(i)=vmod*DCOS(ang)
                              vxpeprox=vmod*DSIN(ang)
                        else
                              vxpe(i)=vxpeprox
                        end if
                        rs1=ran2(iseed)
                        vmod=DSQRT(duekteme*DLOG(rs1))
                        rs2=ran2(iseed)
                        ang=duepi*rs2
                        vype(i)=vmod*DCOS(ang) 
                        vzpe(i)=vmod*DSIN(ang)
                  end if
            end do

            ! Ions    
            duektimi = -2.*kB*Ti0/Mi
            ii=0  
            do i=1,npi
                  Eypi(i)=wyi(i)*Ey(jpi(i)-1)+(1.-wyi(i))*Ey(jpi(i))
                  vypi(i)=vypi(i)+consti*Eypi(i)
                  vzpi(i)=vzpi(i)+consti*Ez0
                  ypi(i)=ypi(i)+vypi(i)*dt
                  zpi(i)=zpi(i)+vzpi(i)*dt
                  ! Periodic boundary conditions
                  if (ypi(i).lt.y(0)) then
                        ypi(i)=y(ny)+ypi(i)
                  else if (ypi(i).ge.y(ny)) then
                        ypi(i)=ypi(i)-y(ny)
                  end if
                  ! Refresh particles
                  if (zpi(i).ge.zch) then
                        ii=ii+1
                        zpi(i)=zch-zacc
                        ! Full-Maxwellian distribution (Box-Muller transformation)         
                        rs1=ran2(iseed)
                        vmod=DSQRT(duektimi*DLOG(rs1))
                        rs2=ran2(iseed)
                        ang=duepi*rs2
                        if ((MOD(ii+1,2)).eq.(MOD(2,2))) then         
                              vxpi(i)=vmod*DCOS(ang)
                              vxpiprox=vmod*DSIN(ang)
                        else
                              vxpi(i)=vxpiprox
                        end if
                        rs1=ran2(iseed)
                        vmod=DSQRT(duektimi*DLOG(rs1))
                        rs2=ran2(iseed)
                        ang=duepi*rs2       
                        vypi(i)=vmod*DCOS(ang) 
                        vzpi(i)=vmod*DSIN(ang)
                  end if
            end do


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

      FUNCTION ran2(i)
      
            ! Random number generator
            ! Long period (> 2.E18) random number generator of l'Ecuyer with Bays-Durham shuffle and added safeguards.
            ! Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
            ! Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
            ! RNMX should approximate the largest floating value that is less than 1.
            
            IMPLICIT REAL*8 (a-h,o-z)
            parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1)
            parameter (ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211)
            parameter (ir2=3791,ntab=32,ndiv=1+imm1/ntab)
            parameter (epsr=1.2e-7,rnmx=1.-epsr)
            integer iv(ntab)
            save iv,iy,idum2
            data idum2/123456789/, iv/ntab*0/, iy/0/
            
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

!*****************************************************************************
!*****************************************************************************
