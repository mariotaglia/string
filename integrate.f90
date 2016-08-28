subroutine integrate
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
integer*4 ier2
real*8 xh(ntot,NS)
real*8 xpot(ntot, 2:NS-1)
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
real*8 xh2(ntot,2:NS-1)
real*8 aa
integer fl(2)
real*8 arc(NS), sumarc(NS), arc_tmp(NS)
real*8 pro0(cuantas,dimx,NS)
real*8 pro_old(cuantas,dimx,NS)
integer iter2
real*8 norma2
real*8 error2
integer xx


do i = 2, NS-1
xxout(i) = float(i-1)/float(NS-1)
enddo

norma = 1.0

do while (norma.lt.error)

iter=iter+1
pro_old = pro

! 1. Calculate avpol from the probabilities
do ii = 2, NS-1 ! loop over nodes

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

call MPI_REDUCE(avpol_tmp(:), avpol2(:,ii), ntot, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
enddo ! ii

call MPI_BCAST(avpol2, ntot*(NS-2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err) ! broadcast to everyone

xh2 = 1.0-avpol2

! 2. Calculation of xpot
do ii = 2, NS-1
do i = 1, ntot
xpot(i,ii) = xh2(i,ii)**(vpol)
 ix=mapx(i)
 iy=mapy(i)
do jx = -Xulimit, Xulimit 
  do jy = -Xulimit, Xulimit 
  kx = ix+jx
  kx= mod(kx-1+50*dimx, dimx) + 1
  ky = iy+jy
  if((ky.le.dimy).and.(ky.ge.1)) then
    k = imap(kx,ky)
    xpot(i,ii) = xpot(i,ii)*dexp(Xu(jx,jy)*avpol2(k,ii)*st/(vpol*vsol))
  endif
  enddo
 enddo
enddo
enddo

! 3. Calculation of pro0

q=0.0d0                   ! init q to zero

do ii = 2, NS-1 ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
pro0(i,xx,ii) = shift

    do j=1, long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     pro(i,xx,ii)= pro(i,xx,ii) * xpot(k,ii)
    enddo

    q(xx,ii)=q(xx,ii)+pro0(i,xx,ii)
enddo ! i
pro0(:,xx,ii) = pro0(:,xx,ii)/q(xx,ii)
enddo ! xx

! 4. Evolve system
! dP/dt = k(P0-P)

do ii = 2, NS-1 ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
pro(i,xx,ii)=STEP*(pro0(i,xx,ii)-pro(i,xx,ii)) ! normalization is conserved
enddo
enddo
enddo


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

arc_temp = arc

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

q = 0.0
do ii = 2, NS-1
do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas
q(xx,ii) = q(xx,ii) + pro(i,xx,ii)
enddo
pro(:,xx,ii) = pro(:,xx,ii)/q(xx,ii)
pro(:,xx,ii) = pro(:,xx,ii)/q(xx,ii)
enddo ! xx
enddo ! ii

! calculate displacement

norma_tmp = 0.0
do ii = 2, NS-1
do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas
norma_tmp=norma_tmp + abs(pro(i,xx,ii)-pro_old(i,xx,ii))
enddo
enddo ! xx
enddo ! ii


call MPI_REDUCE(norma_tmp, norma, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_BCAST(norma, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

if(rank.eq.0) print*, iter, norma

enddo ! norma

end

