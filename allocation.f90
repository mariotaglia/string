subroutine allocation
use brush
use mkinsol
use longs
use kai
use string
use partfunc

allocate (LM(NS-2))
allocate (pp((ntot+1)*(NS-2)))
allocate (avpol(ntot, NS))
allocate (avsol(ntot, NS))
allocate (in1n(cuantas,ntot))
allocate (Xu(-Xulimit:Xulimit))
allocate (pro(cuantas, NS))
allocate (q(NS))
allocate (xfirst(ntot))
allocate (xlast(ntot))
allocate (fixLM(NS-2))
end
