subroutine allocation
use brush
use mkinsol
use longs
use kai
use string
use partfunc
use MPI

allocate (startx(size))
allocate (endx(size))
allocate (avpol(ntot, NS))
allocate (avsol(ntot, NS))
allocate (in1n(cuantas,long))
allocate (Xu(-Xulimit:Xulimit,-Xulimit:Xulimit))
allocate (pro(cuantas, dimx, NS))
allocate (q(dimx,NS))
end
