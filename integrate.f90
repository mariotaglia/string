subroutine integration
use brush
use partfunc
use layer
use volume
use bulk
use longs
use kai
use string
use MPI


implicit none
real*8, external :: LINTERPOL
integer*4 ier2
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
real*8 avpol2(ntot,2:NS-1)
real*8 aa
integer fl(2)
real*8 arc(NS), sumarc(NS), arc_tmp(NS)
real*8 pro0(cuantas,dimx,NS)
real*8 pro_old(cuantas,dimx,NS)
integer iter2
real*8 norma_tmp
real*8 error2
integer xx
real*8 xxout(2:NS-1)
real*8 vin(1:NS)
real*8 xpos
integer fs, ls

do i = 2, NS-1
xxout(i) = float(i-1)/float(NS-1)
enddo

norma = 1.0
q = 0.0

do while (norma.gt.error)

iter=iter+1
pro_old = pro

!print*,'!1!'

! 1. Calculate avpol from the probabilities


if(iter.eq.1) then
 fs = 1
 ls = NS ! calculate avpol for extremes,  only for first iter
else
 fs = 2
 ls = NS-1
endif


do ii = fs, ls ! loop over nodes
avpol_tmp = 0.0
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
    do j=1,long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,xx,ii)*sigma*vsol/delta*vpol ! pro is normed
     ky = dimy-ky+1
     k = imap(kx,ky)
     avpol_tmp(k)=avpol_tmp(k)+pro(i,xx,ii)*sigma*vsol/delta*vpol ! opposing wall
    enddo
enddo ! i
enddo ! xx

if(iter.eq.1) then
  call MPI_REDUCE(avpol_tmp(:), avpol(:,ii), ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
  if((ii.gt.1).and.(ii.lt.NS))avpol2(:,ii) = avpol(:,ii)
else
  call MPI_REDUCE(avpol_tmp(:), avpol2(:,ii), ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
enddo ! ii

if(iter.eq.1) then
  call MPI_BCAST(avpol, ntot*NS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err) ! broadcast to everyone
else
  call MPI_BCAST(avpol2, ntot*(NS-2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err) ! broadcast to everyone
  avpol(:,2:NS-1) = avpol2(:,2:NS-1)
endif

xh = 1.0-avpol

!OJO
!do ii = 2, NS-1
!do i = 1, ntot
!print*, i, ii, avpol2(i,ii), xh2(i,ii)
!enddo
!enddo
!stop


! 2. Calculation of xpot
do ii = fs, ls
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
    xpot(i,ii) = xpot(i,ii)*dexp(Xu(jx,jy)*avpol(k,ii)*st/(vpol*vsol))
  endif
  enddo !jx
 enddo ! jy
enddo ! i
enddo ! ii

!OJO
!do ii = 2, NS-1
!do i = 1, ntot
!print*, i, ii, xpot(i,ii)
!enddo
!enddo
!stop



! 3. Calculation of pro0

q(:,2:NS-1)=0.0d0                   ! init q to zero

do ii = fs, ls ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
pro0(i,xx,ii) = 1.0

    do j=1, long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     pro0(i,xx,ii)= pro(i,xx,ii) * xpot(k,ii)
    enddo

    q(xx,ii)=q(xx,ii)+pro0(i,xx,ii)
enddo ! i
pro0(:,xx,ii) = pro0(:,xx,ii)/q(xx,ii)
enddo ! xx
enddo ! ii

!do i = 1, newcuantas
!print*, i, pro0(i,1,2),pro0(i,1,3),pro0(i,1,4)
!enddo
!stop

! 4. Evolve system
! dP/dt = k(P0-P)

do ii = 2, NS-1 ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
pro(i,xx,ii)=STEP*(pro0(i,xx,ii)-pro(i,xx,ii)) + pro(i,xx,ii) ! normalization is conserved
enddo
enddo
enddo

!do i = 1, newcuantas
!print*, i, pro(i,1,1),pro(i,1,2),pro(i,1,3),pro(i,1,4),pro(i,1,5)
!enddo
!stop

! 5. Calculate distance between beads
arc = 0.0

arc(1) = 0.0 ! position of the fist point

do ii = 2, NS ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
arc(ii) = arc(ii)+ abs(pro(i,xx,ii-1)-pro(i,xx,ii))
enddo
enddo
enddo

arc_tmp = arc

call MPI_REDUCE(arc_tmp, arc, NS-1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_BCAST(arc, NS-1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

sumarc = 0.0

! calculate sum of arc
do ii = 2, NS
do jj = 2,ii
sumarc(ii) = sumarc(ii)+arc(jj) 
enddo
enddo

sumarc = sumarc/sumarc(NS)

!6. Redistribute prob

do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas

 vin(:) = pro(i,xx,:)
do j = 2,NS-1
  xpos = xxout(j)
  pro(i,xx,j) = LINTERPOL (NS, sumarc, vin, xpos , IERR)
enddo

enddo
enddo

! 7. Renorm q 

!q = 0.0
!do ii = 2, NS-1
!do xx = startx(rank+1), endx(rank+1)
!do i = 1, newcuantas
!q(xx,ii) = q(xx,ii) + pro(i,xx,ii)
!enddo
!pro(:,xx,ii) = pro(:,xx,ii)/q(xx,ii)
!enddo ! xx
!print*,'q', ii, q(1,ii)
!enddo ! ii

!do i = 1, newcuantas
!print*,i, pro(i,1,1),pro(i,1,2),pro(i,1,3),pro(i,1,4),pro(i,1,5)
!enddo


! calculate displacement

norma_tmp = 0.0
do ii = 2, NS-1
do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas
norma_tmp=norma_tmp + abs(pro(i,xx,ii)-pro_old(i,xx,ii))
enddo
enddo ! xx
enddo ! ii

!print*, norma_tmp

call MPI_REDUCE(norma_tmp, norma, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_BCAST(norma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

!print*, norma

norma = norma/STEP
if(rank.eq.0) print*, iter, norma
!stop
enddo ! norma

!print*, q

avsol = xh
end

