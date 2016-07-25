
module mkinsol
double precision, allocatable :: pp(:)
endmodule

module kai
integer Xulimit
real*8, allocatable :: Xu(:)
real*8 st
end module

module brush

integer, parameter :: ncha_max = 700
real*8 lseg
real*8 error              ! error imposed accuaracy
real*8 infile             ! inputfile control variable for reading input files  value 0,1
CHARACTER nada

real*8 norma
integer cuantas          ! number of polymer configuration or  bound sequences

integer ntot ! lattice sites
real*8, allocatable :: avpol(:) ! volume fraction polymers already adsorbed

integer*1, allocatable :: in1n(:,:)
real*8 sigma

integer iter              ! counts number of iterations
endmodule

module partfunc
real*8 q
endmodule

module layer
real*8, parameter :: delta = 0.5
endmodule

module volume
real*8 vpol, vsol
endmodule


module bulk
real*8 xsolbulk
endmodule

module seed1
integer seed              ! seed for random number generator
endmodule


module longs
integer long            ! length of polymer
endmodule

module pis
real*8 pi
endmodule

module matrices
real*8 tt(3,3),tp(3,3),tm(3,3)
endmodule

module sines
real*8 sitheta,cotheta,siphip,cophip
endmodule
