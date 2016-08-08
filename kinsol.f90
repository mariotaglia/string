!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
use brush
use mkinsol
use string
implicit none

integer neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vv(*), ftem(*)
integer *4 ier ! Kinsol error flag

common /psize/ neq

neq = (ntot)*(NS-2)

do  i = 1, neq
   vv(i) = vv(i) * pp(i)
enddo
ier = 0

return
end

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * *
!c     The routine kpreco is the preconditioner setup routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
use brush
use mkinsol
use string
implicit none
integer *4 ier ! Kinsol error flag
integer neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vtemp1(*), vtemp2(*)

common /psize/ neq

!   pp(:) = 1.0

do i = 1, (NS-2)*(ntot)
   pp(i) = 1.0  / (1.0+exp(1.0-udata(i)))
enddo
!do i = 1, (NS-2)*ntot+1,(NS-2)*(ntot+1)
!   pp(i) = 1.0  / (1.0+exp(1.0-udata(i)))
!enddo


   ier = 0
return
end

subroutine call_kinsol(x1_old, xg1_old, ier)
use brush
use string

implicit none
integer *4 ier ! Kinsol error flag
integer i
real*8 x1((ntot)*(NS-2)), xg1((ntot)*(NS-2))
real*8 x1_old((ntot)*(NS-2)), xg1_old((ntot)*(NS-2))
integer*8 iout(15) ! Kinsol additional output information
real*8 rout(2) ! Kinsol additional out information
integer*8 msbpre
real*8 fnormtol, scsteptol
real*8 scale((ntot)*(NS-2))
real*8 uscale((ntot)*(NS-2))
real*8 fscale((ntot)*(NS-2))
real*8 constr((ntot)*(NS-2))
integer*4  globalstrat, maxl, maxlrst
integer neq ! Kinsol number of equations
integer*4 max_niter
common /psize/ neq ! Kinsol
integer ierr

neq=(ntot)*(NS-2)

! INICIA KINSOL

msbpre  = 1 ! maximum number of iterations without prec. setup (?)
fnormtol = error ! Function-norm stopping tolerance
scsteptol = error ! Function-norm stopping tolerance

maxl = 2000 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
maxlrst = 0 ! maximum number of restarts
max_niter = 20000
globalstrat = 0

call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
  print*, 'call_kinsol: SUNDIALS_ERROR: FNVINITS returned IER = ', ier
  stop
endif

call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
if (ier .ne. 0) then
   print*, 'call_kinsol: SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
   stop
 endif

call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information

call fkinsetrin('FNORM_TOL', fnormtol, ier)
call fkinsetrin('SSTEP_TOL', scsteptol, ier)
!call fkinsetiin('MAX_NITER', max_niter, ier)

constr = 0.0

do i = 1, (ntot)*(NS-2)  !constraint vector
   constr(i) = 2.0 ! xh > 0
enddo


call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector

! CALL FKINSPTFQMR (MAXL, IER)
call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
!call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG
!CALL FKINDENSE (neq, ier)

if (ier .ne. 0) then
  print*, 'call_kinsol: SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
  call fkinfree ! libera memoria
  stop
endif
call fkinspilssetprec(1, ier) ! preconditiones

do i = 1, neq ! scaling vector
  uscale(i) = 1.0
  fscale(i) = 1.0
enddo

do i = 1, neq ! Initial guess
      x1(i) = x1_old(i)
      xg1(i) = x1(i)  
enddo

call fkinsol(x1, globalstrat, uscale, fscale, ier)         ! Llama a kinsol

print*, 'IER:', ier

if (ier .lt. 0) then
      print*, 'call_kinsol: SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'call_kinsol: Linear Solver returned IER = ', iout(9)
      call fkinfree
!      stop
endif

do i = 1, neq ! output
  x1_old(i) = x1(i)
  xg1_old(i) = x1(i)
enddo

call fkinfree
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol


subroutine call_fkfun(x1_old)
use brush
use string
use string
implicit none
integer i
integer neqs

real*8 x1_old((ntot)*(NS-2))
real*8 x1((ntot)*(NS-2))
real*8 f((ntot)*(NS-2))
integer*4 ier

neqs = (ntot)*(NS-2)
x1 = 0.0
do i = 1,neqs
  x1(i) = x1_old(i)
enddo

call fkfun(x1,f, ier) ! todavia no hay solucion => fkfun 
end

