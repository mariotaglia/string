module intg
use string
integer fs,ls
real*8, allocatable :: xpot(:,:)
real*8, allocatable :: pro0(:,:,:)
real*8, allocatable :: xxout(:)
endmodule


subroutine integration
use maps
use brush
use partfunc
use MPI
use intg

implicit none
integer i, j, xx, ii
real*8 pro_old(cuantas,dimx,NS)
real*8 norma_tmp
integer err

! Allocate commons in intg

allocate(xpot(ntot,NS))
allocate (pro0(cuantas,dimx,NS))
allocate (xxout(2:NS-1))

! Set initial 
do i = 2, NS-1
xxout(i) = float(i-1)/float(NS-1)
enddo

norma = 1.0e10
q = 0.0

do while (norma.gt.error)

iter=iter+1
pro_old = pro

!print*,'!1!'


! fs and ls have the first and last node to consider
if(iter.eq.1) then
 fs = 1
 ls = NS ! calculate avpol for extremes,  only for first iter
else
 fs = 2
 ls = NS-1
endif


! 1. Calculate avpol from the probabilities
! NEED TO OPTIMIZE
if(usecsr.eq.0)call calc_avpol
if(usecsr.eq.1)call calc_avpol_csr

!print*, avpol

! 2. Calculation of xpot
call calc_xpot

! 3. Calculation of pro0
! NEED TO OPTIMIZE
call calc_pro0

! 4. Evolve system
call evolve

! 5. Calculate distance between beads
call redistrib

! 6. Renorm q 
call renormq

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

norma = norma/STEP
if(rank.eq.0) print*, iter, norma
if(mod(iter,100).eq.0)write(10,*)iter, norma
flush(10)

! save probs
if(mod(iter,500).eq.0) then
  call saveprobs
endif ! iter/500

enddo ! norma
end



subroutine calc_avpol ! calculate avpol from probabilities
use maps
use brush
use partfunc
use layer
use volume
use bulk
use longs
use string
use MPI
use intg
implicit none
real*8 avpol_tmp(ntot)
real*8 avpol_tmp2(ntot)
real*8 avpol2(ntot,2:NS-1)
integer ii,xx,k,kx,ky,i, j, kk
integer err

do ii = fs, ls ! loop over nodes
avpol_tmp = 0.0


do xx = startx(rank+1), endx(rank+1)
!print*, ii,xx
avpol_tmp2 = 0.0
do i=1,newcuantas ! loop over cuantas
    do j=1,long
     k = in1n(i,j)
     kx=mapx(k)+(xx-1)
     kx= mod(kx-1+50*dimx, dimx) + 1
     ky=mapy(k)
     k = imap(kx,ky)
     avpol_tmp2(k)=avpol_tmp2(k)+pro(i,xx,ii) ! pro is normed
    enddo
enddo ! i

! opposing wall
do k = 1, ntot
avpol_tmp(k) = avpol_tmp(k) + avpol_tmp2(k)
kx = mapx(k)
ky = mapy(k)
ky = dimy-ky+1
kk = imap(kx,ky)
avpol_tmp(kk) = avpol_tmp(kk) + avpol_tmp2(k)
enddo
enddo ! xx

avpol_tmp = avpol_tmp*sigma*vsol*vpol/delta


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

end subroutine

subroutine calc_xpot
use maps
use brush
use layer
use volume
use bulk
use longs
use kai
use string
use MPI
use intg

implicit none
integer ii,i ,ix,iy,jx,jy,kx,ky,k

avsol = 1.0-avpol
do ii = fs, ls
do i = 1, ntot
xpot(i,ii) = avsol(i,ii)**(vpol)
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
end subroutine


subroutine calc_pro0
use maps
use brush
use partfunc
use string
use MPI
use intg
use longs

implicit none
integer ii, i, xx, j,k, kx,ky

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
     pro0(i,xx,ii)= pro0(i,xx,ii) * xpot(k,ii)
    enddo

    q(xx,ii)=q(xx,ii)+pro0(i,xx,ii)
enddo ! i
pro0(:,xx,ii) = pro0(:,xx,ii)/q(xx,ii)
enddo ! xx
enddo ! ii
end subroutine

subroutine evolve
use maps
use brush
use string
use MPI
use intg

implicit none
integer ii,xx,i

do ii = 2, NS-1 ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
pro(i,xx,ii)=STEP*(pro0(i,xx,ii)-pro(i,xx,ii)) + pro(i,xx,ii) ! normalization is conserved
enddo
enddo
enddo
end subroutine

subroutine redistrib
use maps
use brush
use string
use MPI
use intg


implicit none
real*8 arc(NS), sumarc(NS), arc_tmp(NS)
integer ii,xx,i, jj, j
real*8, external :: LINTERPOL
real*8 vin(1:NS)
real*8 xpos
integer err

arc = 0.0
arc(1) = 0.0 ! position of the fist point

do ii = 2, NS ! loop over beads
do xx = startx(rank+1), endx(rank+1)
do i=1,newcuantas ! loop over cuantas
arc(ii) = arc(ii)+ abs(pro(i,xx,ii-1)-pro(i,xx,ii))
enddo
enddo ! xx
enddo ! ii

arc_tmp = arc

call MPI_REDUCE(arc_tmp, arc, NS, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
call MPI_BCAST(arc, NS, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, err)

sumarc = 0.0

! calculate sum of arc
do ii = 2, NS
do jj = 2,ii
sumarc(ii) = sumarc(ii)+arc(jj) 
enddo
enddo

sumarc = sumarc/sumarc(NS)

! Redistribute prob

do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas

 vin(:) = pro(i,xx,:)
do j = 2,NS-1
  xpos = xxout(j)
  pro(i,xx,j) = LINTERPOL (NS, sumarc, vin, xpos , IERR)
enddo

enddo
enddo

end subroutine



subroutine renormq
use maps
use brush
use partfunc
use string
use MPI
use intg


implicit none
integer ii, xx,i
real*8 qq(dimx,NS)

qq = 0.0
do ii = 2, NS-1
do xx = startx(rank+1), endx(rank+1)
do i = 1, newcuantas
qq(xx,ii) = qq(xx,ii) + pro(i,xx,ii)
enddo
pro(:,xx,ii) = pro(:,xx,ii)/qq(xx,ii)
enddo ! xx
enddo ! ii
end subroutine


subroutine saveprobs
use maps
use brush
use string
use MPI
use intg
implicit none
real*8 pro_tosend(cuantas,dimx)
integer ii, ix, i

pro_tosend = 0.0
do ii = 2,NS-1
do ix = startx(rank+1), endx(rank+1)
pro_tosend(:,ix) = pro(:,ix,ii)
enddo

! recover pro()
call MPI_REDUCE(pro_tosend(:,:), pro(:,:,ii), dimx*cuantas, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
enddo ! ii

if(rank.eq.0) then
 do ii = 1, NS
 do i=1,newcuantas
 do ix = 1, dimx
  write(1900+ii,*)pro(i,ix,ii)
 enddo
 enddo
 close(1900+ii)
 enddo
endif ! rank

end subroutine
