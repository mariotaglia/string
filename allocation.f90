subroutine allocation
use brush
use mkinsol
use longs
use kai
use string
use partfunc
use MPI
use maps

allocate (startx(size))
allocate (endx(size))
allocate (avpol(ntot, NS))
allocate (avsol(ntot, NS))
allocate (in1n(cuantas,long))
allocate (Xu(-Xulimit:Xulimit,-Xulimit:Xulimit))
allocate (pro(cuantas, dimx, NS))
allocate (q(dimx,NS))
allocate (imap(dimx,dimy))
allocate (mapx(dimx*dimy))
allocate (mapy(dimx*dimy))
end
