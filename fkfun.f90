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
real*8 x(ntot*(NS-2)),f(ntot*(NS-2))
real*8 xh(ntot,NS)
integer i,j,k1,k2,ii,kk, jj,iz       ! dummy indices
integer err
integer n
real*8 algo, algo1, algo3
real*8 algo2
real*8 xtotal(1-Xulimit:ntot+Xulimit, NS)
real*8 aa
integer fl(2)
real*8 error2

! zspow
real*8 xgpow(2*(NS-2)), x1pow(2*(NS-2))
integer*4 ier

error2 = 1.0d-6
shift = 1.0
n = ntot

! Retrive values from input vector
! solvent from 1 to NS-2
do ii = 1,NS-2
do i=1,ntot                 
xh(i,ii+1)=x(i+(ii-1)*ntot)  ! solvent density=volume fraction   
enddo
enddo

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for LM and beta using zspow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii=1, NS-2
xgpow(ii) = LM(ii)
xgpow(ii+NS-2) = beta(ii)
enddo

x1pow = xgpow

iter2 = 1
call call_zspow(x1pow, x1pow, ier, 2*(NS-2))
stop

do ii=1, NS-2
LM(ii) = x1pow(ii) 
beta(ii) = x1pow(ii+NS-2)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! contruction of f and the volume fractions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 1, NS-2 ! Packing constraint
do i=1,ntot
 f(i+(ii-1)*ntot)=xh(i,ii+1)+avpol(i,ii+1)-1.0d0
enddo
enddo

iter=iter+1

algo1 = 0.0
do i = 1, (NS-2)*(ntot)
 algo1 = algo1 + f(i)**2
end do

!PRINT*, iter, algo1+algo2+algo3, LM(1), beta(1)
!PRINT*, iter, algo1+algo2+algo3,algo1,algo2,algo3, LM(1)
norma=algo1
print*, 'Outer Loop:', iter, norma



return
end

subroutine Wfunction(pro,prop,LM,beta) ! solves for pro(i,ii) using W function, see notes
implicit none
real*8 pro,prop,LM,beta
real*8 arg
integer*4 nb, l, nerror
real*8, external :: wapr
real*8 cutoff
real*8 pro0

real*8 qpro, qprop, qLM, qbeta, qpro0
real*8 L1, L2
real*8 lnarg
real*8 tmp

qpro0 = pro
cutoff = 20.0 ! minimum z to use asymptotic expression for W = exp(z), for z < w use WAPR 

qprop = prop
qbeta = beta
qLM = exp(LM) ! (log(LM)+1)

nb = 0 ! upper branch, LM<0
l = 0 ! not offset

lnarg=log(qLM*qpro0)+(qLM*qbeta*qprop)

if(lnarg.gt.cutoff) then
 L1 = lnarg
 L2 = log(lnarg)
 tmp = L1-L2+L2/L1
 tmp = tmp + L2*(-2.0+L2)/2/L1**2
 tmp = tmp + L2*(6.0-9.0*L2+2*L2**2)/6/L1**3
 tmp = tmp + L2*(-12.0+36.0*L2-22.0*L2**2+3.0*L2**3)/12.0/L1**4
 tmp = tmp + L2*(60.0-300.0*L2+350*L2**2-125*L2**3+12.0*L2**4)/60.0/L1**5
else ! use WAPR
 arg = exp(lnarg)
 tmp = wapr(arg,nb,nerror,l)
endif



if(nerror.eq.1) then
  print*,'Wfunction: Error in WAPR'
  print*, arg, LM, pro, prop, beta
  stop
endif

!print*, 'tmp', tmp
!print*, 'pro0', pro0
qpro = tmp/qLM
pro=qpro
end subroutine

subroutine FCN(xin, fout, nin, par)
use string, only : NS
use brush
use partfunc
use layer
use volume
use longs

implicit none
integer nin
real*8 xin(nin), fout(nin)
real*8 par(1)
integer *4 ier
integer ii,i
real*8 LM0(NS-2), LM1(NS-2)
real*8 beta0(NS-2), beta1(NS-2) ! q/q'
real*8 pro0(cuantas,NS)
real*8 arc(NS-1), sumarc, arc0
real*8 norma2
integer j,jj

! retrive LM and beta from input
do ii=1, NS-2
LM0(ii) = xin(ii)
beta0(ii) = xin(ii+NS-2)
print*, ii, LM0(ii), beta0(ii)
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  LOOP OVER STRING BEADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ii = 2, NS-1
jj = ii - 1 ! Lagrange multiplier index for ii

q(ii) = 0.0
avpol(:,ii) = 0.0

do i=1,newcuantas ! loop over cuantas

pro(i,ii) = shift

    do j=1, ntot
     pro(i,ii)= pro(i,ii) * xpot(j, ii)**in1n(i,j)
    enddo

!    if(pro(i,ii).eq.0.0) then
!      print*, i, ii
!      stop
!    endif
    pro0(i,ii) = pro(i,ii)
    call Wfunction(pro(i,ii),pro(i,ii-1), LM0(jj), beta0(jj)) ! solves for pro(i,ii) using W function, see notes

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
! norm by q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construction of arclenght vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sumarc = 0.0
do ii = 1, NS-1
arc(ii) = 0.0
do i = 1, newcuantas
arc(ii) = arc(ii) + abs(pro(i,ii)-pro(i,ii+1))
enddo
sumarc = sumarc + arc(ii)
enddo ! ii
arc0 = sumarc / (float(NS)-1.0) ! target arclenght

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do ii = 2, NS-1
jj=ii-1
LM1(jj) = 0.0
!print*, 'arc', jj, arc(jj), arc0
do i = 1, newcuantas
LM1(jj) = LM1(jj) + abs(log(pro(i,ii)*q(ii)/pro0(i,ii)))
enddo
LM1(jj) = log(LM1(jj)/arc0/beta0(jj))
beta1(jj)=q(ii)
enddo

norma2 = 0.0

do ii=1,NS-2
fout(ii) = 1-LM1(ii)/LM0(ii)
norma2 = norma2 + abs((LM1(ii)-LM0(ii))/LM1(ii))
enddo
do ii=1,NS-2
fout(ii+NS-2) = 1-beta1(ii)/beta0(ii)
norma2 = norma2 + abs((beta1(ii)-beta0(ii))/beta1(ii))
enddo

print*,LM0(1), LM1(1)
print*,beta0(1), beta1(1)
print*,'      Inner Loop:', iter2, norma2
iter2 = iter2 + 1

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subrutina que llama a zspow
!

subroutine call_zspow(x1, xg1, ier,ntot)
implicit none
      
integer ntot
integer NSIG
integer ITEMAX
real*8 par(1)
real*8 WK(((ntot)*(3*(ntot)+15))/2)
real*8 x1(ntot), xg1(ntot)
real*8 fnorm
integer*4 ier
integer ier_zspow
integer i

! IMSL

NSIG = 6
ITEMAX = 200

x1 = xg1

call ZSPOW(NSIG,ntot,ITEMAX,PAR, x1, fnorm, WK, ier_zspow)
print*, ier_zspow
stop
ier = -ier_zspow ! errores invertidos
return
end




