subroutine fkfun(x,f,ier2)
use brush
use partfunc
use layer
use volume
use bulk
use longs
use kai
implicit none
integer*4 ier2
real*8 x(ntot),f(ntot)
real*8 xh(ntot)
real*8 xpot(ntot)
integer i,j,k1,k2,ii, jj,iz       ! dummy indices
integer err
integer n
real*8 algo
real*8 algo2
real*8 xtotal(ntot)
integer ix,iy,jx,jy,kx,ky,k
integer, external :: imap, mapx, mapy
real*8 avpol_tmp(ntot)

shift = 1.0d100

n = ntot

do i=1,n                 
xh(i)=x(i)  ! solvent density=volume fraction   
enddo

xtotal = 0.0
do i = 1,n
xtotal(i) = 1.0-xh(i)
enddo

avpol = 0.0 
avpol_tmp = 0.0 

! calculation of xpot
do i = 1, ntot
xpot(i) = xh(i)**(vpol)
 ix=mapx(i)
 iy=mapy(i)
do jx = -Xulimit, Xulimit 
  do jy = -Xulimit, Xulimit 
  kx = ix+jx
  kx= mod(kx-1+50*dimx, dimx) + 1
  ky = iy+jy
  if((ky.le.dimy).and.(ky.ge.1)) then
    k = imap(kx,ky)
    xpot(i) = xpot(i)*dexp(Xu(jx,jy)*xtotal(k)*st/(vpol*vsol))
  endif
enddo
enddo
enddo

!    probability distribution

do ii = 1, dimx ! loop over grafting pos
q=0.0d0                   ! init q to zero
avpol_tmp=0.0
do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, long
     k = in1n(i,j)
     kx=mapx(k)+(ii-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     pro(i,ii)= pro(i,ii) * xpot(k)
    enddo

    q(ii)=q(ii)+pro(i,ii)

    do j=1,long
     k = in1n(i,j)
     kx=mapx(k)+(ii-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,ii)*sigma*vsol/delta*vpol
     ky = dimy-ky+1
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,ii)*sigma*vsol/delta*vpol ! opposing wall
    enddo
 
enddo ! i

! norm by q

avpol_tmp = avpol_tmp/q(ii)
avpol = avpol + avpol_tmp

enddo ! ii
! contruction of f and the volume fractions

do i=1,n
 f(i)=xh(i)+avpol(i)-1.0d0
enddo

iter=iter+1

algo = 0.0
do i = 1, n
 algo = algo + f(i)**2
end do

algo2 = 0.0  
do i = 1, n
 algo2 = algo2+avpol(i)
enddo
algo2 = algo2*delta/(vpol*vsol)/2/sigma


PRINT*, iter, algo, algo2
norma=algo

return
end
