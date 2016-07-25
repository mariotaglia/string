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
integer i,j,k1,k2,ii, jj,iz       ! dummy indices
integer err
integer n
real*8 algo
real*8 algo2
real*8 xtotal(1-Xulimit:ntot+Xulimit)
real*8 LM(NS-2) ! Lagrange multipliers
real*8 aa

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
!  LOOP OVER STRING BEADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 2, NS-1

!SEGUIR!

do i=1,newcuantas ! loop over cuantas
pro(i) = shift

    do j=1, ntot
     pro(i)= pro(i) * xpot(j)**in1n(i,j)
    enddo

    q=q+pro(i)

    do j=1, ntot
       avpol(j)=avpol(j)+pro(i)*sigma*vsol/delta*vpol*in1n(i,j)
       avpol(j)=avpol(j)+pro(i)*sigma*vsol/delta*vpol*in1n(i,ntot-j+1) ! opposing wall
    enddo
 
enddo ! i

! norm by q

avpol = avpol/q

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
