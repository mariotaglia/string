subroutine read
use longs
use brush
use bulk
use kai
use string
implicit none
integer i

!     reading in of variables from stdin

read(8,*)nada
read(8,*)NS0, NS

read(8,*)nada
read(8,*)FIX

read(8,*),nada
read(8,*),dimx, dimy
ntot = dimx*dimy

read(8,*),nada
read(8,*)cuantas

read(8,*),nada
read(8,*)long

READ(8,*),nada
read(8,*),nsigma
read(8,*), (sigmas(i), i=1, nsigma)

READ(8,*),nada
read(8,*),error, STEP

READ(8,*),nada
read(8,*),infile

read(8,*)nada
read(8,*)lseg

read(8,*)nada
read(8,*)Xulimit

read(8,*)nada
read(8,*)nst
read(8,*), (sts(i), i=1, nst)
end
