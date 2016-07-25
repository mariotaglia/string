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
real*8 xtotal(1-Xulimit:ntot+Xulimit)

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

! calculation of xpot
do i = 1, ntot
xpot(i) = xh(i)**(vpol)
do j = -Xulimit, Xulimit 
 xpot(i) = xpot(i)*dexp(Xu(j)*xtotal(i+j)*st/(vpol*vsol))
enddo
enddo

!    probability distribution
q=0.0d0                   ! init q to zero

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
