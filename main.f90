!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Introductory code: just a neutral brush!
!
! Serial version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call read ! read input from input.txt
call allocation ! allocate arrays in memory
call kai ! generate poor solvent
call creador
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
use string
implicit none

integer ncha
integer *4 ier ! Kinsol error flag
real*8 pi ! pi
real*8 Na ! avogadros' number              
parameter (Na=6.02d23)
integer cc,ccc

real*8 x1((ntot+1)*(NS-2)),xg1((ntot+1)*(NS-2))   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 

integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

integer i,j,k,m,ii,flag,c, jj ! dummy indices

INTEGER temp
real*8 tempr
real*8 tmp

real*8 algo, algo2                  

real*8 zp(long)

real*8 sum,sumel          ! auxiliary variable used in free energy computation  
real*8 sumpi,sumrho,sumrhopol, sumrho2, sumrho2mol !suma de la fraccion de polimero

! single layer files
character*18 sysfilename      ! contains value of free energy, input parameter etc
character*26 denssolfilename  ! contains the denisty of the solvent
character*28 denspolfilename

integer countfile         ! enumerates the outputfiles 
integer countfileuno     ! enumerates the outputfiles para una corrida

integer readsalt          !integer to read salt concentrations

seed=435               ! seed for random number generator

print*, 'Program Simple Brush'
print*, 'GIT Version: ', _VERSION


!     initializations of variables 
pi=dacos(-1.0d0)          ! pi = arccos(-1) 
itmax=200                 ! maximum number of iterations       
n=ntot                    ! size of lattice

vsol=0.030                ! volume solvent molecule in (nm)^3
vpol= ((4.0/3.0)*pi*(0.3)**3)/vsol  ! volume polymer segment in units of vsol

! bulk

xsolbulk = 1.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read first and last
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(unit=20, file='first.dat')
open(unit=21, file='last.dat')

do i=1,n
read(20,*)tempr,xfirst(i) 
read(21,*)tempr,xlast(i) 
enddo   

close(20)
close(21)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  make initial guess by interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do ii = 1, NS-2
do i=1,n
xg1(i+(ii-1)*n) = (xlast(i)-xfirst(i))*float(ii)/float(NS-1) + xfirst(i)
enddo
enddo

!!! Initial guess for LM = 0
do ii = 1, NS-2
xg1(ntot*(NS-2)+ii) = 0 !-0.27213976464426065
enddo

x1 = xg1

do i = 1,n
zc(i)= (i-0.5) * delta
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     computation starts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     initializations of input depended variables 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do cc = 1, nst
do ccc = 1, nsigma

st=sts(cc)
sigma=sigmas(ccc)

countfile=cc+ccc-1

iter=0                    ! iteration counter

open(unit=10, file='out.out')

! Call solver 

   iter = 0
   print*, 'solve: Enter solver ', (NS-2)*(ntot+1), ' eqs'
   call call_kinsol(x1, xg1, ier)

do ii = 1, NS-2
do i=1,n
avsol(i,ii+1)=xg1(i+(ii-1)*ntot)  ! solvent density=volume fraction
enddo
enddo
do ii =1, NS-2
print*, 'LM',ii, xg1(ntot*(NS-2)+ii)
enddo

avsol(:,1)=xfirst(:)
avsol(:,NS)=xlast(:)

do ii = 1, NS

countfileuno = ii

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', countfileuno,'.',countfile,'.dat'
write(denspolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitypolymer.',countfileuno,'.',countfile,'.dat'
write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', countfileuno,'.',countfile,'.dat'

open(unit=310,file=sysfilename)
open(unit=321,file=denspolfilename)
open(unit=330,file=denssolfilename)

do i=1,n
write(321,*)zc(i),avpol(i, ii)
write(330,*)zc(i),avsol(i, ii)
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
write(310,*)'index  = ',ii

close(310)
CLOSE(321)
close(330)

call free_energy(ii)

enddo ! ii

enddo
enddo

stop

end

