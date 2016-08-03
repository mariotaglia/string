subroutine allocation
use brush
use mkinsol
use longs
use kai
use partfunc
allocate (pp(ntot))
allocate (avpol(ntot))
allocate (avsol(ntot))
allocate (in1n(cuantas,long))
allocate (Xu(-Xulimit:Xulimit,-Xulimit:Xulimit))
allocate (pro(cuantas,dimx))
allocate (q(dimx))
end
