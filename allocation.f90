subroutine allocation
use brush
use mkinsol
use longs
use kai

allocate (pp(ntot))
allocate (avpol(ntot))
allocate (in1n(cuantas,ntot))
allocate (Xu(-Xulimit:Xulimit))
end
