!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Introductory code: just a neutral brush!
!
! Serial version
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use seed1
call initmpi
call read ! read input from input.txt
call allocation ! allocate arrays in memory
call MPIdist ! distribute grafting points among processors
call kai ! generate poor solvent
seed = 435
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
use MPI
implicit none

real*8 suma
real*8, external :: LINTERPOL
real*8 xpos
integer ix, iy
integer, external :: imap, mapx, mapy

integer ncha
integer *4 ier ! Kinsol error flag
real*8 pi ! pi
real*8 Na ! avogadros' number              
parameter (Na=6.02d23)
integer cc,ccc
real*8 vin(NS0), vout(NS-2)
real*8 xinput(ntot,NS0)
real*8 xoutput(ntot,NS)
real*8 xxin(NS0), xxout(NS-2)
integer*4 NI, NO
real*8 x1((ntot+1)*(NS-2)),xg1((ntot+1)*(NS-2))   ! density solvent iteration vector
real*8 zc(ntot)           ! z-coordinate layer 
integer mid
integer n                 ! number of lattice sites
integer itmax             ! maximum number of iteration allowed for 
real*8 fnorm              ! L2 norm of residual vector function fcn

integer i,j,k,m,ii,flag,c, jj ! dummy indices

INTEGER temp, tempx, tempy
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


integer tag, source
parameter(tag = 0)
integer err
integer ierror


seed=435               ! seed for random number generator

if(rank.eq.0)print*, 'Program Simple Brush'
if(rank.eq.0)print*, 'GIT Version: ', _VERSION


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

pro = 0.0

open(unit=30, file='firstp.dat')
open(unit=31, file='lastp.dat')

do ix = 1, dimx
do i=1,newcuantas
read(30,*)pro(i,ix,1)
read(31,*)pro(i,ix,NS)
enddo
enddo

close(21)
close(31)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read probs from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 2, NS0-1
do i=1,newcuantas
do ix = 1, dimx
read(800+ii,*, IOSTAT=ierror)pro(i,ix,ii)
enddo
enddo ! i
enddo ! ix
close(800+cc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  make initial guess by interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)print*, 'Interpolation'

do i = 1, NS0
xxin(i) = float(i-1)/float(NS0-1)
enddo
do i = 2, NS-1
xxout(i-1) = float(i-1)/float(NS-1)
enddo


NI = NS0
NO = NS-2

do i = 1, newcuantas
do ix = 1, dimx
 vin(1) = pro(i,ix,1)
 vin(NS0) = pro(i,ix,NS)
 do j = 2, NS0-1
  vin(j) = pro(i,ix,j)
 enddo
do j = 1,NS-2
  xpos = xxout(j)
  vout(j) = LINTERPOL (NI, xxin, vin, xpos , IERR)
enddo

do j = 2, NS-1
pro(i,ix,j) = vout(j-1)
enddo

enddo
enddo

!!!!!!!!!!!!!! Auxiliary fields... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!do i = 1, newcuantas
!print*, i, pro(i,1,1),pro(i,1,3),pro(i,1,5)
!enddo

! OJO
!do ii = 1, NS
!print*, ii, pro(1,1,ii)
!enddo 
!do ii = 1, NS
!print*, ii, pro(10,1,ii)
!enddo 
!stop


call integration


! OK


do ii = 1, NS

call free_energy(ii)
countfileuno = ii

write(sysfilename,'(A7,BZ,I3.3,A1,I3.3,A4)')'system.', countfileuno,'.',countfile,'.dat'
write(denspolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitypolymer.',countfileuno,'.',countfile,'.dat'
write(denssolfilename,'(A15,BZ,I3.3,A1,I3.3,A4)')'densitysolvent.', countfileuno,'.',countfile,'.dat'

open(unit=310,file=sysfilename)
open(unit=321,file=denspolfilename)
open(unit=330,file=denssolfilename)

do i=1,n
ix = mapx(i)
iy = mapy(i)

write(321,*)ix,iy,avpol(i, ii)
write(330,*)ix,iy,avsol(i, ii)
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


enddo ! ii

enddo
enddo
call MPI_FINALIZE(ierr) ! finaliza MPI
stop

end


subroutine MPIdist
use MPI
use brush
implicit none
integer ppc
integer i
ppc = int(dimx/size)
if(rank.eq.0)print*,'ppc', ppc
do i = 1, size
startx(i) = 1 + ppc*(i-1)
endx(i) = startx(i)+ppc-1
enddo
endx(size) = dimx
do i = 1, size
if(rank.eq.0)print*,'Proc', i,'start', startx(i), 'end',endx(i)
enddo
end 

