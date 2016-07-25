!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Introductory code: just a neutral brush!
!
! Serial version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call read ! read input from input.txt
call allocation ! allocate arrays in memory
call kai ! generate poor solvent
call solve ! solve the system
end

subroutine solve
use brush
use partfunc
use layer
use volume
use bulk
use seed1
use longs
use kai
implicit none

integer ncha
integer *4 ier ! Kinsol error flag
real*8 pi ! pi
real*8 Na ! avogadros' number              
parameter (Na=6.02d23)
integer cc,ccc

real*8 x1(ntot),xg1(ntot)   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

integer i,j,k,m,ii,flag,c, jj ! dummy indices

INTEGER temp
real*8 tempr
real*8 tmp

REAL*8 xfile(ntot)                        
real*8 algo, algo2                  

real*8 chains(3,long,ncha_max) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
real*8 chainsw(ncha_max), sumweight_tosend
real*8 zp(long)

real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*26 denssolfilename  ! contains the denisty of the solvent
character*28 denspolfilename

integer countfile         ! enumerates the outputfiles 
integer countfileuno     ! enumerates the outputfiles para una corrida
integer conf              ! counts number of conformations

integer readsalt          !integer to read salt concentrations

seed=435               ! seed for random number generator

print*, 'Program Simple Brush'
print*, 'GIT Version: ', _VERSION


!     initializations of variables 
pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice
conf=0                    ! counter for conformations

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol= ((4.0/3.0)*pi*(0.3)**3)/vsol  ! volume polymer segment in units of vsol

! bulk

xsolbulk = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     init guess all 1.0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,n
xg1(i)=1.0
x1(i)=1.0
zc(i)= (i-0.5) * delta
enddo

!     init guess from files fort.100 (solvent) and fort.200 (potential)                      

if (infile.ge.1) then
do i=1,n
read(100,*)tempr,x1(i) 
xg1(i) = x1(i)  ! solvent
enddo   
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, 'Calling RIS chain generator'
in1n = 0

call initcha              ! init matrices for chain generation
conf=0                    ! counter of number of conformations

newcuantas = 0

do while (conf.lt.cuantas)
call cadenas1(chains,ncha) ! generate only chains with first segment at z > 0

do j=1,ncha
  if(conf.lt.cuantas) then
   conf=conf+1

   flag = 0
      do k=1,long
        temp=int(chains(1,k,j)/delta)+1  ! put segments into the correct layer
           if(temp.gt.ntot)flag = 1
      enddo ! k

      if(flag.eq.0) then
      newcuantas=newcuantas+1
       do k=1,long
        temp=int(chains(1,k,j)/delta)+1  ! put segments into the correct layer
        in1n(newcuantas,temp) =  in1n(newcuantas,temp) + 1
       enddo ! k
      endif
   endif
enddo ! j
enddo ! while

print*,"Chains ready"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     computation starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initializations of input depended variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do cc = 1, nst
do ccc = 1, nsigma

st=sts(cc)
sigma=sigmas(ccc)

countfileuno=cc
countfile=ccc

iter=0                    ! iteration counter

do i=1,n             ! initial gues for x1
xg1(i)=x1(i)
enddo

! Call solver 

   iter = 0
   print*, 'solve: Enter solver ', ntot, ' eqs'
   call call_kinsol(x1, xg1, ier)

do i=1,n
avsol(i)=x1(i) ! retrive xsol from solution
enddo

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', countfileuno,'.',countfile,'.dat'
write(denspolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitypolymer.',countfileuno,'.',countfile,'.dat'
write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', countfileuno,'.',countfile,'.dat'

open(unit=310,file=sysfilename)
open(unit=321,file=denspolfilename)
open(unit=330,file=denssolfilename)

do i=1,n
write(321,*)zc(i),avpol(i)
write(330,*)zc(i),avsol(i)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     additional system information
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(310,*)'system      = neutral polymer'
write(310,*)'fnorm       = ', norma ! residual size of iteration vector
write(310,*)'error       = ',error
write(310,*)'q           = ',q
write(310,*)'length seg  = ',0.50 ! value see subroutine cadenas
write(310,*)'delta       = ',delta
write(310,*)'vsol        = ',vsol
write(310,*)'vpol        = ',vpol*vsol

write(310,*)'cuantas     = ',cuantas
write(310,*)'newcuantas     = ',newcuantas
write(310,*)'iterations  = ',iter

close(310)
CLOSE(321)
close(330)

call free_energy

enddo
enddo

stop

end

