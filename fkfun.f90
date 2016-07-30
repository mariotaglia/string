subroutine fkfun(x,f,ier2)
use brush
use partfunc
use layer
use volume
use bulk
use longs
use kai
use string

implicit none
integer*4 ier2
real*8 x((ntot+1)*(NS-2)),f((ntot+1)*(NS-2))
real*8 xh(ntot,NS)
real*8 xpot(ntot, NS)
integer i,j,k1,k2,ii,kk, jj,iz       ! dummy indices
integer err
integer n
real*8 algo, algo1, algo3
real*8 algo2
real*8 xtotal(1-Xulimit:ntot+Xulimit, NS)
real*8 LM0(NS-2)
real*8 aa
integer fl(2)
real*8 arc(NS-1), sumarc, arc0
real*8 pro0(cuantas,NS)
integer iter2
real*8 norma2
real*8 error2

error2 = 1.0d-5
shift = 1.0
n = ntot

! Retrive values from input vector
! solvent from 1 to NS-2
do ii = 1,NS-2
do i=1,ntot                 
xh(i,ii+1)=x(i+(ii-1)*ntot)  ! solvent density=volume fraction   
!print*, i,x(i+(ii-1)*ntot)
enddo
enddo

do ii = 1, NS-2
if (FIX.ne.1) then   
   LM0(ii) = (x(ntot*(NS-2)+ii))
else
   LM0(ii) = fixLM(ii)
endif
enddo

!LM0(ii) = exp(x(ntot*(NS-2)+ii))
!print*, ii,LM0(ii)

! Retrive solvent for first and last
do i =1, n
xh(i,1) = xfirst(i)
xh(i,NS) = xlast(i)
enddo

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

ii = fl(kk)
jj = ii - 1 ! Lagrange multiplier index for ii

do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, ntot
     pro(i,ii)= pro(i,ii) * xpot(j,ii)**in1n(i,j)
    enddo

    q(ii)=q(ii)+pro(i, ii)

    do j=1, ntot
       avpol(j,ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,j)
       avpol(j,ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,ntot-j+1) ! opposing wall
    enddo
 
enddo ! i

avpol(:,ii) = avpol(:,ii)/q(ii)
pro(:,ii) = pro(:, ii)/q(ii)

enddo ! kk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LOOP OVER STRING BEADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 2, NS-1

! update xpot

do j = 1, ntot
xpot(j,ii) = xpot(j,ii)*exp(LM0(ii-1)*(xtotal(j,ii)-xtotal(j,ii-1)))
enddo

jj = ii - 1 ! Lagrange multiplier index for ii

q(ii) = 0.0
avpol(:,ii) = 0.0

do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, ntot
     pro(i,ii)= pro(i,ii) * xpot(j, ii)**in1n(i,j)
    enddo

    q(ii)=q(ii)+pro(i,ii)

    do j=1, ntot
       avpol(j, ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,j)
       avpol(j, ii)=avpol(j, ii)+pro(i, ii)*sigma*vsol/delta*vpol*in1n(i,ntot-j+1) ! opposing wall
    enddo
 
enddo ! i

avpol(:,ii) = avpol(:,ii)/q(ii)
pro(:,ii) = pro(:,ii)/q(ii)
enddo ! ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of arclenght vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sumarc = 0.0
do ii = 1, NS-1
arc(ii) = 0.0
do j = 1, ntot
arc(ii) = arc(ii) + (avpol(j,ii+1)-avpol(j,ii))**2
enddo
sumarc = sumarc + arc(ii)
enddo ! ii
arc0 = sumarc / (float(NS)-1.0) ! target arclenght

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contruction of f and the volume fractions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 1, NS-2 ! Packing constraint
do i=1,ntot
 f(i+(ii-1)*ntot)=xh(i,ii+1)+avpol(i,ii+1)-1.0d0
enddo
enddo

do ii = 1, NS-2
if(FIX.ne.1) then
   f(ntot*(NS-2)+ii) = -(arc(ii)-arc0)/arc0
else
   f(ntot*(NS-2)+ii) = 0.0
endif
enddo

iter=iter+1

algo1 = 0.0
do i = 1, (NS-2)*(ntot)
 algo1 = algo1 + f(i)**2
end do
algo2 = 0.0
do i = (NS-2)*(ntot)+1, (NS-2)*(ntot+1)
 algo2 = algo2 + f(i)**2
enddo

norma=algo1
print*, 'Outer Loop:', iter, algo1, algo2, norma
!print*, LM0(1), arc0, arc(1)
flush(10)

!do i = 1, (NS-2)*(ntot+1)
!print*, 'out', i, f(i)
!enddo
ier2 = 0
return
end

