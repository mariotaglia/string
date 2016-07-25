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
real*8 x((ntot+1)*(NS-2)),f((ntot+1)*(NS-2))
real*8 xh(ntot,NS)
real*8 xpot(ntot, NS)
integer i,j,k1,k2,ii,kk, jj,iz       ! dummy indices
integer err
integer n
real*8 algo
real*8 algo2
real*8 xtotal(1-Xulimit:ntot+Xulimit)
real*8 LM(NS-2) ! Lagrange multipliers
real*8 aa
integer fl(2)
real*8 arc(NS-1), sumarc, arc0

shift = 1.0d100
n = ntot

! Retrive values from input vector
! solvent from 1 to NS-2
do ii = 1,NS-2
do i=1,ntot                 
xh(i,ii+1)=x(i+(ii-1)*ntot)  ! solvent density=volume fraction   
enddo
enddo
! LM from ntot*(NS-2) + 1 to ntot*(NS-2) + NS - 2
do ii = 1, NS-2
LM(ii) = x(ntot*(NS-2)+ii)
enddo

! Retrive solvent for first and last
xh(:,1) = xfirst(:)
xh(:,NS) = xlast (:)

xtotal = 0.0
do ii = 1, NS
do i = 1,n
xtotal(i, ii) = 1.0-xh(i, ii)
enddo
enddo

! init variables
avpol = 0.0 
avpol(:,1)=1-xh(:,1)
avpol(:,NS)=1-xh(:,NS)

! calculation of xpot
do ii = 1, NS
do i = 1, ntot
xpot(i,ii) = xh(i,ii)**(vpol)
do j = -Xulimit, Xulimit 
 xpot(i,ii) = xpot(i,ii)*dexp(Xu(j)*xtotal(i+j,ii)*st/(vpol*vsol))
enddo
enddo
enddo

!    probability distribution
q=0.0d0                   ! init q to zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALCULATE FIRST AND LAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fl(1) = 1 ! first
fl(2) = NS ! last

do kk = 1, 2

ii = fl(ii)
jj = ii - 1 ! Lagrange multiplier index for ii

do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, ntot
     pro(i,ii)= pro(i,ii) * xpot(j)**in1n(i,j)
    enddo

    q(ii)=q(ii)+pro(i, ii)

    do j=1, ntot
       avpol(j,ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,j)
       avpol(j,ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,ntot-j+1) ! opposing wall
    enddo
 
enddo ! i

enddo ! kk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LOOP OVER STRING BEADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 2, NS-1
jj = ii - 1 ! Lagrange multiplier index for ii
!SEGUIR!

do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, ntot
     pro(i,ii)= pro(i,ii) * xpot(j)**in1n(i,j)
    enddo

    call Wfunction(pro(i,ii),pro(i,ii-1), LM(jj)) ! solves for pro(i,ii) using W function, see notes

    q(ii)=q(ii)+pro(i,ii)

    do j=1, ntot
       avpol(j, ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,j)
       avpol(j, ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,ntot-j+1) ! opposing wall
    enddo
 
enddo ! i

enddo ! ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! norm by q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 1, NS
avpol(:,ii) = avpol(:,ii)/q(ii)
pro(:,ii) = pro(:, ii)/q(ii)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of arclenght vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sumarc = 0.0
do ii = 1, NS-1
arc(ii) = 0.0
do i = 1, newcuantas
arc(ii) = arc(ii) + (pro(i,ii)-pro(i,ii+1))**2
enddo
arc(ii) =sqrt(arc(ii))
sumarc = sumarc + arc(ii)
enddo ! ii
arc0 = sumarc / (float(NS)-1.0) ! target arclenght

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contruction of f and the volume fractions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do ii = 1,NS-2
do i=1,ntot                 
xh(i,ii+1)=x(i+(ii-1)*ntot)  ! solvent density=volume fraction   
enddo
enddo
! LM from ntot*(NS-2) + 1 to ntot*(NS-2) + NS - 2
do ii = 1, NS-2
LM(ii) = x(ntot*(NS-2)+ii)
enddo

do ii = 1, NS-2 ! Packing constraint
do i=1,ntot
 f(i+(ii-1)*ntot)=xh(i,ii+1)+avpol(i,ii+1)-1.0d0
enddo
enddo

do ii = 1, NS-2
 f(ntot*(NS-2)+ii) = arc(ii)-arc0
enddo

iter=iter+1

algo = 0.0
do i = 1, (NS-2)*(ntot+1)
 algo = algo + f(i)**2
end do

PRINT*, iter, algo
norma=algo

return
end
