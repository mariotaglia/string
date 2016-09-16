module csr
real*8, allocatable :: inc_values(:,:)
integer, allocatable :: inc_columns(:,:)
integer, allocatable :: pntrb(:,:)
integer, allocatable :: pntre(:,:)
integer, allocatable :: gidx(:)
endmodule


subroutine in2csr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine transforms the in1n matrix into a csr format 
! to use the intel mkl libraries
!
! note that there is a csr matrix per grafting point
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use longs
use brush
use maps
use mpi
use csr
implicit none

real*8 sparse(dimx*dimy) ! space matrix for a given conformation
real*8 nse(dimx*dimy) ! non-sparse element matrix for a given conformation
integer nsp(dimx*dimy) ! non-sparse position matrix for a given conformation
integer pos, posx, posy
integer nonzero ! number of non-zeros in sparse matrix
integer i, j, k, idx, xx, xxx
integer flag

if(rank.eq.0) print*, 'Create CSR matrixes for mkl'

nonzero = 0


! DEBUG1. Print in1n
!do i = 1, newcuantas
!do j = 1, long
!print*, i,j,in1n(i,j)
!enddo
!enddo
!stop

! 1. Determine non-zeros
! it is the same for all processors and grafting points

do i = 1, newcuantas
 sparse = 0.0
 do j = 1, long
  if(sparse(in1n(i,j)).eq.0.0) then
    nonzero = nonzero + 1
  endif
    sparse(in1n(i,j)) = sparse(in1n(i,j)) + 1.0
 enddo ! j
! DBUG : PRINT SPARSE
!print*, i, sparse
enddo ! i
!stop

!2. Allocate arrays for csr format

ALLOCATE (inc_values(nonzero, ppc(rank+1)))
ALLOCATE (inc_columns(nonzero, ppc(rank+1)))
ALLOCATE (pntrb(newcuantas, ppc(rank+1)))
ALLOCATE (pntre(newcuantas, ppc(rank+1)))
ALLOCATE (gidx(ppc(rank+1)))

!3. Dump in1n into csr arrays

do xx = startx(rank+1), endx(rank+1) ! loop over grafting positions of processor rank
xxx = xx - startx(rank+1) + 1 ! index of csr matrix

gidx(xxx) = 0 ! global index for matrix xxx

do i = 1, newcuantas
 nse = 0.0 
 nsp = 0
 idx = 0
 do j = 1, long
   flag = 0

! shift to grafting point
   pos = in1n(i,j)
   posx = mapx(pos)+(xx-1)
   posx = mod(posx-1+50*dimx,dimx) + 1
   posy = mapy(pos)
   pos = imap(posx,posy) 
!

   do k = 1, idx 
     if(nsp(k).eq.pos) then ! position already saved
       nse(k) = nse(k) + 1.0 ! increment by one
       flag = 1 ! will not increase index
       exit ! exit k loop
     endif
   enddo ! k
   if(flag.eq.0) then ! increase index
    idx = idx + 1
    nse(idx) = 1.0
    nsp(idx) = pos
   endif 
 enddo ! j

! print*, i
! print*, nse
! print*, nsp

! nse and nsp now contains the element and position of conformation i

 pntrb(i,xxx) = gidx(xxx) + 1 ! position of first element in the row in global list  

do j = 1, idx ! values to add
 gidx(xxx) = gidx(xxx) + 1 ! increase total counter by one
 inc_values(gidx(xxx),xxx) = nse(j) ! store value
 inc_columns(gidx(xxx),xxx) = nsp(j) ! store position (1...dimx*dimy)
enddo ! j




enddo ! i

pntre(1:newcuantas-1,xxx) = pntrb(2:newcuantas,xxx)
pntre(newcuantas,xxx) = gidx(xxx) + 1

enddo ! xx

!print*, nonzero
!print*, gidx(1)
!print*, inc_values
!print*, inc_columns
!print*, pntrb
!print*, pntre
!stop


end subroutine


subroutine calc_avpol_csr ! calculate avpol from probabilities
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
use csr


implicit none
real*8 avpol_tmp(ntot)
real*8 avpol_tmp2(ntot)
real*8 avpol2(ntot,2:NS-1)
integer ii,xx,k,kx,ky,i, j, kk,xxx
integer err
CHARACTER matdescra(6)

matdescra(1) = "G"
matdescra(4) = "F"


do ii = fs, ls ! loop over nodes
avpol_tmp = 0.0

do xx = startx(rank+1), endx(rank+1)
xxx = xx - startx(rank+1) + 1 ! index of csr matrix

avpol_tmp2 = 0.0

 call mkl_dcsrmv('T', newcuantas, dimx*dimy, 1.0d+0, matdescra, inc_values(:,xxx), inc_columns(:,xxx), pntrb(:,xxx), pntre(:,xxx), & 
 pro(:,xx,ii), 0.0d+0, avpol_tmp2)

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

subroutine calc_pro0_csr
use maps
use brush
use partfunc
use string
use MPI
use intg
use longs
use csr

implicit none
integer ii, i, xx, j,k, kx,ky,xxx
real*8 lnpro(newcuantas)
real*8 lnxpot(ntot)
CHARACTER matdescra(6)


matdescra(1) = "G"
matdescra(4) = "F"

q(:,2:NS-1)=0.0d0                   ! init q to zero

do ii = fs, ls ! loop over beads

lnxpot(:) = log(xpot(:,ii))

do xx = startx(rank+1), endx(rank+1)
xxx = xx - startx(rank+1) + 1 ! index of csr matrix

call mkl_dcsrmv('N', newcuantas, dimx*dimy, 1.0d+0, matdescra, &
 inc_values(:,xxx), inc_columns(:,xxx), pntrb(:,xxx), pntre(:,xxx), &
 lnxpot, 0.0d+0, lnpro)
   
pro0(:,xx,ii) = exp(lnpro(:))

q(xx,ii) = sum(pro0(:,xx,ii))
pro0(:,xx,ii) = pro0(:,xx,ii)/q(xx,ii)

enddo ! xx
enddo ! ii
end subroutine





