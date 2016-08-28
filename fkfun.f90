subroutine fkfun(x,f,ier2)
use brush
use partfunc
use layer
use volume
use bulk
use longs
use kai
use string
use segregated
implicit none
integer*4 ier2
real*8 x(ntot+1),f(ntot+1)
real*8 xh(ntot,NS)
real*8 xpot(ntot, NS)
integer i,j,k1,k2,ii,kk, jj,iz       ! dummy indices
integer err
integer n
real*8 algo, algo1, algo3
real*8 algo2
real*8 xtotal(ntot, NS)
integer ix,iy,jx,jy,kx,ky,k
integer, external :: imap, mapx, mapy
real*8 avpol_tmp(ntot)
real*8 LM0(NS-2)
real*8 aa
integer fl(2)
real*8 arc1,arc2
real*8 pro0(cuantas,NS)
integer iter2
real*8 norma2
real*8 error2
integer xx

error2 = 1.0d-5
shift = 1.0
n = ntot

! Retrive value from xoutput
xh = xoutput
LM0 = LMoutput

! Retrive values from input vector
! Overwrites xoutput for NS_current
ii = NS_current
do i=1,ntot                 
xh(i,ii)=x(i)  ! solvent density=volume fraction   
enddo
LM0(ii-1) = x(ntot+1)

! calculate xtotal and avpol for all

xtotal = 0.0
do ii = 1, NS
do i = 1,n
xtotal(i, ii) = 1.0-xh(i, ii)
avpol(i, ii) = 1.0-xh(i, ii)
enddo
enddo

! init variables
ii = NS_current
avpol(:,ii) = 0.0 
q(:,ii)=0.0d0                   ! init q to zero

! calculation of xpot
do ii = 1, NS
do i = 1, ntot

xpot(i,ii) = xh(i,ii)**(vpol)
 ix=mapx(i)
 iy=mapy(i)
do jx = -Xulimit, Xulimit 
  do jy = -Xulimit, Xulimit 
  kx = ix+jx
  kx= mod(kx-1+50*dimx, dimx) + 1
  ky = iy+jy
  if((ky.le.dimy).and.(ky.ge.1)) then
    k = imap(kx,ky)
    xpot(i,ii) = xpot(i,ii)*dexp(Xu(jx,jy)*xtotal(k,ii)*st/(vpol*vsol))
  endif
  enddo
 enddo
enddo
enddo

ii = NS_current ! add LM
do j = 1, ntot
xpot(j,ii) = xpot(j,ii)*exp(LM0(ii-1)*(xtotal(j,ii+1)-xtotal(j,ii-1)))
enddo

jj = ii - 1 ! Lagrange multiplier index for ii

do xx = 1, dimx
avpol_tmp = 0.0
do i=1,newcuantas ! loop over cuantas

pro(i,xx,ii) = shift

    do j=1, long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     pro(i,xx,ii)= pro(i,xx,ii) * xpot(k,ii)
    enddo

    q(xx,ii)=q(xx,ii)+pro(i,xx,ii)

    do j=1,long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,xx,ii)*sigma*vsol/delta*vpol
     ky = dimy-ky+1
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,xx,ii)*sigma*vsol/delta*vpol ! opposing wall
    enddo
 
enddo ! i

avpol(:,ii) = avpol(:,ii) +avpol_tmp(:)/q(xx,ii)
pro(:,xx,ii) = pro(:,xx,ii)/q(xx,ii)
enddo ! xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of arclenght vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

arc2 = 0.0
arc1 = 0.0 
do j = 1, ntot
arc1 = arc1 + (avpol(j,NS_current)-avpol(j,NS_current-1))**2
arc2 = arc2 + (avpol(j,NS_current+1)-avpol(j,NS_current))**2
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contruction of f and the volume fractions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ii = NS_current
do i=1,ntot
 f(i)=xh(i,ii)+avpol(i,ii)-1.0d0
enddo

f(ntot+1) = -(arc1-arc2)/(arc1+arc2)

iter=iter+1

algo1 = 0.0
do i = 1, ntot+1
 algo1 = algo1 + f(i)**2
end do
if(iter.eq.1)normaini = algo1

norma=algo1
print*, 'Outer Loop:', iter, norma, 'NS_current', NS_current
!print*, LM0(1), arc0, arc(1)
flush(10)

!do i = 1, (NS-2)*(ntot+1)
!print*, 'out', i, f(i)
!enddo
ier2 = 0
return
end

