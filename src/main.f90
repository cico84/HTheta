!******************************************************************************          
!*                                                                            *
!*                1D{theta}-3V FULLY KINETIC PIC SIMULATION                   *
!*                    OF ACCELERATION REGION IN SPT-100                       *               
!*                       (LOCAL FIELD APPROXIMATION)                          *
!*                                                                            *
!******************************************************************************



!  GOAL

!  Simulare la dinamica del plasma in un settore del dominio azimutale
!  Studio dell'instabilità di deriva ExB

 
!  RECORD DELLE REVISIONI

!    Data        Programmatore               Descrizione delle modifiche
!  ========    ==================    ========================================================
!  **/03/18    Francesco Taccogna    Codice originale per periodo di ricerca Asadi-Bogopolsky
           
  
!  CHARACTERISTICS OF MODEL:

!  1. Modello monodimensionale cartesiano nella coordinata azimutale: theta -> y
!  2. Local field approximation (densXe,Ez,Br -> fixed value)
!  3. Ciascuna macroparticella e' una carica per unita' di superficie
!  4. Griglia spaziale uniforme
!  5. PIC per gli elettroni e gli ioni (funzione peso "linear weighting")
!  6. Risoluzione eq. di Poisson (SOR+Chebyschev)
!  7. Risoluzione eq. di Maxwell (gradiente a due punti) 
!  8. Integrazione numerica delle equazioni del moto (metodo leapfrog) 
!  9. e/Xe collisions: Monte Carlo (Xsect from Szabo)

! NOMENCLATURE:

! ny             = # punti di griglia 
! y(j)           = coordinata del punto di griglia
! ype/i(i)       = coordinata della macroparticella
! vxpe/i(i)      = componente x (radiale) della velocita' della macroparticella
! vype/i(i)      = componente y (azimutale) della velocita' della macroparticella
! vzpe/i(i)      = componente z (assiale) della velocita' della macroparticella
! rhoi(j)        = densita' di carica ionica sui punti di griglia
! rhoe(j)        = densita' di carica elettronica sui punti di griglia
! phi(j)         = potenziale elettrico sui punti di griglia
! apoi,bpoi,cpoi -> coeff. dell'eq. di Poisson in coordinate cartesiane {y} 
! Ey(j)          = componente y del campo elettrico sui punti di griglia
! Eype/i(i)      = comp. y del campo elettrico nella posizione della macropart.


!******************************************************************************
!*                             DEFINITION     PHASE                           *
!******************************************************************************

!     elenco moduli

MODULE const_phys
! Physical and mathematical constants
IMPLICIT NONE
SAVE
DOUBLE PRECISION, PARAMETER :: pi	  = acos(-1.d0)		! Pi Greek
DOUBLE PRECISION, PARAMETER :: eps0	  = 8.854188d-12  	! Vacuum Dielectric Constant
DOUBLE PRECISION, PARAMETER :: q	  = 1.602189d-19  	! Electric charge [C]
DOUBLE PRECISION, PARAMETER :: me	  = 9.1093897d-31 	! Electron mass [Kg]
DOUBLE PRECISION, PARAMETER :: Mi	  = 2.1801714d-25	! Xenon atom Mass (Xe) [Kg]
DOUBLE PRECISION, PARAMETER :: kB	  = 1.380662d-23  	! Boltzmann Constant [J/K]
DOUBLE PRECISION, PARAMETER :: JtoeV  = 6.24146d18 		! Conversion factor from J to eV
DOUBLE PRECISION, PARAMETER :: eVtoK  = 11600        	! Conversion factor from eV to K
DOUBLE PRECISION, PARAMETER :: w 	  = 5.d8			! Macroparticle's statistical weight
DOUBLE PRECISION, PARAMETER :: dt	  = 1.d-12			! Timestep
DOUBLE PRECISION, PARAMETER :: conste = q*dt/(2.*me)
DOUBLE PRECISION, PARAMETER :: consti = q*dt/Mi
DOUBLE PRECISION, PARAMETER :: duepi  = 2.*pi
DOUBLE PRECISION, PARAMETER :: wq     = w*q
END MODULE const_phys

MODULE system
! Physical external parameters
IMPLICIT NONE
SAVE
DOUBLE PRECISION, PARAMETER :: yd	= 0.05		! lunghezza dominio assiale
DOUBLE PRECISION, PARAMETER :: zch  = 0.025		! posizione assiale interfaccia regione di ionizzazione-accelerazione
DOUBLE PRECISION, PARAMETER :: zacc	= 0.015		! lunghezza dominio radiale
DOUBLE PRECISION, PARAMETER :: Ez0  = 20000		! corrente di scarica
DOUBLE PRECISION, PARAMETER :: Br0  = 0.018		! valore massimo del modulo del campo magnetico
DOUBLE PRECISION, PARAMETER :: n0   = 5.e17		! densità massima di plasma stimata
DOUBLE PRECISION, PARAMETER :: Te0  = 116000	! temperatura con la quale gli elettroni entrano nel canale acceleratore
DOUBLE PRECISION, PARAMETER :: Ti0  = 2320		! temperatura propellente neutro	
END MODULE system

MODULE grid
USE const_phys
USE system
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: ny=5000
DOUBLE PRECISION   :: dy,duedy,y(0:ny),vol(0:ny)
END MODULE grid

MODULE poi
USE grid
IMPLICIT NONE
SAVE
DOUBLE PRECISION   :: rhoe(0:ny),rhoi(0:ny),phi(0:ny),Ey(0:ny),apoi,bpoi,cpoi,dpoi(0:ny)
END MODULE poi

MODULE part
USE grid
IMPLICIT NONE
SAVE
INTEGER, PARAMETER :: npmax=500000000,npic=10000000
INTEGER 		   :: npe,npi
DOUBLE PRECISION   :: ype(npmax),zpe(npmax),vxpe(npmax),vype(npmax),vzpe(npmax)
DOUBLE PRECISION   :: ypi(npmax),zpi(npmax),vxpi(npmax),vypi(npmax),vzpi(npmax)
DOUBLE PRECISION   :: wye(npmax),Eype(npmax),wyi(npmax),Eypi(npmax)
DOUBLE PRECISION   :: vxpeprox,vxpiprox                                    
INTEGER            :: jpe(npmax),jpi(npmax)
END MODULE part

MODULE diagn
IMPLICIT NONE
SAVE
DOUBLE PRECISION :: mob,Ee_ave,Ei_ave,Debye,omegae,omegace,vthetae,vzi,vthe,CFLe1,CFLe2
INTEGER			 :: npinit
END MODULE diagn

MODULE rand
IMPLICIT NONE
SAVE
INTEGER			 :: iseed
DOUBLE PRECISION :: rs,rs1,rs2
END MODULE rand


program SPT1D_theta
      
USE poi
USE part
USE diagn
USE rand
            

!******************************************************************************
!******************************************************************************


!******************************************************************************
!*                                                                            *
!*                       INITIALIZATION       PHASE                           *
!*                                                                            *  
!******************************************************************************

!     plasma parameters estimations
!     Debye length
      Debye=dsqrt(eps0*kB*Te0/(n0*q*q))   
!     plasma/cyclotron frequency
      omegae=dsqrt(n0*q*q/(me*eps0))
      omegace=q*Br0/me
!     velocità elettronica di deriva azimutale
      vthetae=Ez0/Br0
!     velocità ionica assiale
      vzi=dsqrt(2.*q*Ez0*zacc/Mi)
!     number of initial electrons/ions (every particle is a charge per unit surface C/m^2)      
      npinit=nint(yd*n0/w)         


!******************************************************************************

!     mesh generation
      call mesh

!     constraint CFL      
      vthe=dsqrt(8.*kB*Te0/(pi*me))
      CFLe1=vthe*dt/dy
      CFLe2=vthetae*dt/dy     
      
      open (3,file='parameter.res',status='unknown')
           write (3,*) 'lD,dy=       ',Debye,dy
           write (3,*) 'dtp,dtB,dt=  ',1./omegae,1./omegace,dt
           write (3,*) 'CFLe=  	     ',CFLe1,CFLe2
           write (3,*) 'vthe,vde,vzi=',vthe,vthetae,vzi
           write (3,*) '# part,ppc=  ',npinit,float(npinit)/float(ny-1)
      close (3)      
      
!******************************************************************************

!     Poisson coefficients: cartesian coordinate {y}
!     apoi*phi(j-1)+bpoi*phi(j)+cpoi*phi(j+1)=dpoi(j)=-(rhoi(j)-rhoe(j))*dy**2/eps0
!     condizioni al contorno: periodicity: phi(0)=phi(ny)

      apoi=1.
      cpoi=apoi
      bpoi=-2.

!     Reference potential at y=y(0)=y(ny)
      phi(ny)=0.
      phi(0)=0.


!******************************************************************************   

!     initialization of random numbers generator
      iseed=-1
      gen=ran2(iseed)


!******************************************************************************

!     initial particle distribution
      call init
      

!******************************************************************************      

!     PIC cycle

      do ipic=1,npic

         call scatter

!        Poisson equation source term
         do j=0,ny
            dpoi(j)=-(rhoi(j)-rhoe(j))*dy**2/eps0
         end do
         
         call fieldsolve

         call push

!        diagnostics: averaged energy; mobility
         Ee_ave=0.
         mob=0.
         do i=1,npe
            Ee_ave=Ee_ave+(vxpe(i)**2+vype(i)**2+vzpe(i)**2)
            mob=mob+vzpe(i)
         end do
         Ee_ave=Ee_ave*JtoeV*0.5*me/npe
         mob=mob/(npe*Ez0)
         
         Ei_ave=0.
         do i=1,npi
            Ei_ave=Ei_ave+(vxpi(i)**2+vypi(i)**2+vzpi(i)**2)
         end do
         Ei_ave=Ei_ave*JtoeV*0.5*Mi/npi

         if (mod(ipic,1).eq.0) then
         open (10,file='history.res',status='unknown',position='append')
         write(10,101) ipic*dt,Ee_ave,Ei_ave,mob,phi(ny/2),Ey(ny/2),rhoe(ny/2)/q,rhoi(ny/2)/q
         close (10)
         end if
         
 101  format (8(2x,1pg13.5e3))
      
!     end PIC cycle
      end do

!******************************************************************************
!******************************************************************************

                             
      end program


!******************************************************************************
!******************************************************************************


!******************************************************************************        
!******************************************************************************
!*                                                                            *
!*                              ELENCO SUBROUTINES                            *
!*                                                                            *
!******************************************************************************        
!******************************************************************************


!******************************************************************************
!******************************************************************************
        
      subroutine mesh 
          
      USE grid
      IMPLICIT NONE
      INTEGER :: j
     
      y(0)=0.
      dy=yd/ny
      duedy=2.*dy
      do j=1,ny
         y(j)=y(j-1)+dy
      end do
      
!     cell volumes
      do j=0,ny
         vol(j)=dy
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
     
      npe=0
      npi=0
    
      duekteme = -2.*kB*Te0/me
      duektimi = -2.*kB*Ti0/Mi
      
!     numero iniziale di elettroni ed ioni: uniform distribution
      do i=1,npinit
       
!      electrons
       npe=npe+1
       rs=ran2(iseed)
       ype(npe)=rs*y(ny)
       rs=ran2(iseed)
       zpe(npe)=(zch-zacc)+rs*zacc
!      full Maxwellian distribution (Box-Muller transformation)         
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
      
!      ions:
       npi=npi+1
       ypi(npi)=ype(npe)
       rs=ran2(iseed)
       zpi(npi)=(zch-zacc)+rs*zacc
!      drifted Maxwellian distribution (Box-Muller transformation): vzi        
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
    
!     charge density reset
      do j=0,ny
         rhoe(j)=0.
         rhoi(j)=0.
      end do

!     electron charge deposite on the mesh points 
      do i=1,npe     
!        charge density weighting (linear weighting)
         jpe(i)=int(ype(i)/dy)+1         
         wye(i)=(y(jpe(i))-ype(i))/dy     
         rhoe(jpe(i)-1)=wye(i)*wq+rhoe(jpe(i)-1)
         rhoe(jpe(i))=(1.-wye(i))*wq+rhoe(jpe(i))        
      end do      
!     periodic boundary condition 
      rhoe(0)=rhoe(0)+rhoe(ny)   
      rhoe(ny)=rhoe(0)
      do j=0,ny
         rhoe(j)=rhoe(j)/vol(j)
      end do

!     ion charge deposite on the mesh points 
      do i=1,npi       
!        charge density weighting (linear weighting) 
         jpi(i)=int(ypi(i)/dy)+1         
         wyi(i)=(y(jpi(i))-ypi(i))/dy     
         rhoi(jpi(i)-1)=wyi(i)*wq+rhoi(jpi(i)-1)
         rhoi(jpi(i))=(1.-wyi(i))*wq+rhoi(jpi(i))                   
      end do       
!     periodic boundary condition 
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
 
!     Poisson equation in cartesian coordinate {theta}->{y}
!     apoi*phi(j-1)+bpoi*phi(j)+cpoi*phi(j+1)=dpoi(j)
!     Thomas tridiagonal algorithm for periodic boundary conditions: phi(ny)=phi(0)=0.
!     [M. Napolitano, Comm. Appl. Num. Meth. 1, 11-15 (1985)]
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

!     Electric field equation solution
      do j=1,ny-1
         Ey(j)=-(phi(j+1)-phi(j-1))/duedy
      end do
!     Periodic boundary conditions
      Ey(0)=-(phi(1)-phi(ny-1))/duedy
      Ey(ny)=Ey(0)
      

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

!     electrons      

!     leapfrog method (Boris algorithm)
      tB=conste*Br0
      tBB=2.*tB/(1.+tB**2)
      duekteme = -2.*kB*Te0/me

      ie=0
      do i=1,npe
       Eype(i)=wye(i)*Ey(jpe(i)-1)+(1.-wye(i))*Ey(jpe(i))
!      first half acceleration by electric field
       vyea=vype(i)-conste*Eype(i)         
       vzea=vzpe(i)-conste*Ez0
!      full rotation around magnetic field
       vyeb=vyea-vzea*tB
       vzeb=vzea+vyea*tB 
       vyec=vyea-vzeb*tBB
       vzec=vzea+vyeb*tBB
!      second half acceleration by electric field  
       vype(i)=vyec-conste*Eype(i)                
       vzpe(i)=vzec-conste*Ez0
!      coordinates updates (cartesian approximation)
       ype(i)=ype(i)+vype(i)*dt
       zpe(i)=zpe(i)+vzpe(i)*dt
!      periodic boundary conditions
       if (ype(i).le.y(0)) then
          ype(i)=y(ny)+ype(i)
       else if (ype(i).ge.y(ny)) then
          ype(i)=ype(i)-y(ny)
       end if
!      refresh particles
       if ((zch-zpe(i)).ge.zacc) then
          ie=ie+1
          zpe(i)=zch
!         drifting-Maxwellian distribution (Box-Muller transformation)         
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

!     ions    
      duektimi = -2.*kB*Ti0/Mi
      ii=0  
      do i=1,npi
       Eypi(i)=wyi(i)*Ey(jpi(i)-1)+(1.-wyi(i))*Ey(jpi(i))
       vypi(i)=vypi(i)+consti*Eypi(i)
       vzpi(i)=vzpi(i)+consti*Ez0
       ypi(i)=ypi(i)+vypi(i)*dt
       zpi(i)=zpi(i)+vzpi(i)*dt
!      periodic boundary conditions
       if (ypi(i).lt.y(0)) then
          ypi(i)=y(ny)+ypi(i)
       else if (ypi(i).ge.y(ny)) then
          ypi(i)=ypi(i)-y(ny)
       end if
!      refresh particles
       if (zpi(i).ge.zch) then
          ii=ii+1
          zpi(i)=zch-zacc
!         drifting-Maxwellian distribution (Box-Muller transformation)         
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
!*                                  ELENCO FUNCTIONS                          *
!*                                                                            *
!******************************************************************************


!*****************************************************************************       

      FUNCTION ran2(i)
      
!     Random number generator
!     Long period (> 2.E18) random number generator of l'Ecuyer with Bays-Durham shuffle and added safeguards.
!     Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values).
!     Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
!     RNMX should approximate the largest floating value that is less than 1.
      
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
