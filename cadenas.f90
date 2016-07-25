subroutine creador

use brush
implicit none
integer flag
integer ncha,j,k
integer temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHAIN GENERATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, 'Calling RIS chain generator'
in1n = 0

call initcha              ! init matrices for chain generation
conf=0                    ! counter of number of conformations

newcuantas = 0

do while (conf.lt.cuantas)
call cadenas1(chains,ncha) ! generate only chains with first segment at z > 0

do j=1,ncha
  if(conf.lt.cuantas) then
   conf=conf+1

   flag = 0
      do k=1,long
        temp=int(chains(1,k,j)/delta)+1  ! put segments into the correct layer
           if(temp.gt.ntot)flag = 1
      enddo ! k

      if(flag.eq.0) then
      newcuantas=newcuantas+1
       do k=1,long
        temp=int(chains(1,k,j)/delta)+1  ! put segments into the correct layer
        in1n(newcuantas,temp) =  in1n(newcuantas,temp) + 1
       enddo ! k
      endif
   endif
enddo ! j
enddo ! while

print*,"Chains ready"
end subroutine


subroutine mrrrr(a,b,c)
real*8 a(3,3),b(3,3),c(3,3)

do i=1,3
 do j=1,3
  c(i,j)=0
 enddo
enddo

do i=1,3
 do j=1,3
  do k=1,3
    c(i,j)=c(i,j)+a(i,k)*b(k,j)
  enddo
 enddo
enddo

return
end 

subroutine initcha
use pis
use matrices
use sines
implicit none

pi=acos(-1.0000000e0)
sitheta=sin(68.0*pi/180.0)
cotheta=cos(68.0*pi/180.0)
siphip=sin(120.0*pi/180.0)
cophip=cos(120.0*pi/180.0)

tt(1,1)=cotheta
tt(1,2)=sitheta
tt(1,3)=0.0
tt(2,1)=sitheta
tt(2,2)=-cotheta
tt(2,3)=0.0
tt(3,1)=0.0
tt(3,2)=0.0
tt(3,3)=-1.0

tp(1,1)=cotheta
tp(1,2)=sitheta
tp(1,3)=0.0
tp(2,1)=sitheta*cophip
tp(2,2)=-cotheta*cophip
tp(2,3)=siphip
tp(3,1)=sitheta*siphip
tp(3,2)=-cotheta*siphip
tp(3,3)=-cophip

tm(1,1)=cotheta
tm(1,2)=sitheta
tm(1,3)=0.0
tm(2,1)=sitheta*cophip
tm(2,2)=-cotheta*cophip
tm(2,3)=-siphip
tm(3,1)=-sitheta*siphip
tm(3,2)=cotheta*siphip
tm(3,3)=-cophip
return
end


!**************************************************************
!* rotates a given chains conformation                        *  
!* pre: xend = input chain                                    *
!* post: xendr= rotated chain                                 *
!*        test= if equal N the rotation took the chain to z<0 *
!*              and has to be rejected as allowd conf         * 
!**************************************************************     

subroutine rota(xend,xendr,n,test)
use seed1
use layer
implicit none

integer n
real*8 xend(3,n+5),rands,xendr(3,n+5)
real*8 theta,theta1,x(120),y(120),r,xp(120),yp(120),pi
integer test
integer k,m1,i
real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
real*8 alfa,gama,cga,a,b,c

pi=acos(-1.0)

fac=rands(seed)
fac1=rands(seed)
fac2=rands(seed)
alfa=fac*2*pi
cbe=fac1*2.0-1.0
gama=fac2*2*pi

sbe=(1-cbe**2)**0.5
cal=cos(alfa)
sal=sin(alfa)
cga=cos(gama)
sga=sin(gama)

test=1

do i=1,n              ! rotation segmentos

a=xend(1,i)
b=xend(2,i)
c=xend(3,i)

xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal

if(xendr(1,i).lt.0.0) then
!print*, i, xendr(1,i)
test=0  ! segment is at z < 0  
endif

xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe

enddo
return
end


!*************************************************************
!* Generates a chains conformations bases                    *
!* on a three state RIS-model see Flory book                 *
!* GENERA CADENAS DE PAH-Os                                  *
!*************************************************************
subroutine cadenas1(chains,ncha)
use longs
use seed1
use pis
use matrices
use sines
use brush
implicit none

real*8, external :: interaction11
integer i,state,ii,j,k1,k2,ncha
real*8 rn,dista
real*8 rands,angle
real*8 m(3,3), mm(3,3)
real*8 x(3),xend(3,long+5),xendr(3,long+5)
real*8 tmp
REAL*8 chains(3,long,ncha_max), chainsw(ncha_max)
integer LT
integer test
REAL*8 tolerancia    !tolerancia en el calculo de selfavoiding


tolerancia = 1.0e-5

223  xend(1,1)=0.0      ! first position 
xend(2,1)=0.0
xend(3,1)=0.0
rn=rands(seed)
angle=0.0

! transition probalibity
m(1,1)=cotheta          
m(1,2)=sitheta
m(1,3)=0.0
m(2,1)=cos(angle)*sitheta
m(2,2)=-cos(angle)*cotheta
m(2,3)=sin(angle)
m(3,1)=sin(angle)*sitheta
m(3,2)=-sin(angle)*cotheta
m(3,3)=-cos(angle)

x(1)=m(1,1)*lseg     
x(2)=m(2,1)*lseg
x(3)=m(3,1)*lseg

xend(1,2)=xend(1,1)+x(1)  ! second postion
xend(2,2)=xend(2,1)+x(2)
xend(3,2)=xend(3,1)+x(3)

do i=3,long          ! loop over remaining positions!

123     rn=rands(seed)
state=int(rn*3)        ! random select the state= {trans,gauch+,gauch-}

if (state.eq.3) then 
state=2
endif

if (state.eq.0) then
!*********************************** TRANS     
  call mrrrr(m,tt,mm)
  do ii=1,3
  do j=1,3
  m(ii,j)=mm(ii,j)
  enddo
  enddo
elseif (state.eq.1) then
!********************************** GAUCHE +
  call mrrrr(m,tp,mm)
  do ii=1,3
  do j=1,3
  m(ii,j)=mm(ii,j)
  enddo
  enddo
elseif (state.eq.2) then
!********************************** GAUCHE -
  call mrrrr(m,tm,mm)
  do ii=1,3
  do j=1,3
  m(ii,j)=mm(ii,j)
  enddo
  enddo
endif

x(1)=m(1,1)*lseg
x(2)=m(2,1)*lseg
x(3)=m(3,1)*lseg

xend(1,i)=xend(1,i-1)+x(1)   ! ith postion chain
xend(2,i)=xend(2,i-1)+x(2)
xend(3,i)=xend(3,i-1)+x(3)

enddo

dista=0.0                       ! check self avoiding constraint (segmentos)
do k1=1,long
do k2=k1+1,long
dista=(xend(1,k2)-xend(1,k1))**(2.0)
dista=dista+(xend(2,k2)-xend(2,k1))**(2.0)
dista=dista+(xend(3,k2)-xend(3,k1))**(2.0)
dista=sqrt(dista)+tolerancia
if (dista.lt.lseg) then
goto 223
endif
enddo
enddo

ncha=0

do i=1,12
call rota(xend,xendr,long,test)   ! rotate chain conformation ncha time

if (test.eq.1) then ! only those with all segments at z > 0
ncha=ncha+1
do j=1,long
 chains(1,j,ncha)=xendr(1,j)       ! output 
 chains(2,j,ncha)=xendr(2,j)
 chains(3,j,ncha)=xendr(3,j)
!print*, chains(1,j,ncha), j, ncha
enddo
endif

enddo

if (ncha.eq.0) goto 223

return
end



