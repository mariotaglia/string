subroutine allocation
use brush
use mkinsol
use longs
use kai

allocate (pp(ntot))
allocate (avpol(ntot))
allocate (avsol(ntot))
allocate (in1n(cuantas,ntot))
allocate (Xu(-Xulimit:Xulimit))
allocate (pro(cuantas))
end
