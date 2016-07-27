subroutine free_energy(cc)
use brush
use layer
use volume
use bulk
use kai
use partfunc
implicit none

integer cc
real*8 F_tot, F_tot2
real*8 F_Mix_s, F_Conf, F_vdW
integer iz, i, j
real*8 sumpi, sumrho
real*8 xtotal(1-Xulimit:ntot+Xulimit)
real*8 mupol

xtotal = 0.0
do i = 1,ntot
xtotal(i) = 1.0-avsol(i, cc)
enddo


F_tot = 0.0
F_tot2 = 0.0

! 1. Mezcla solvente

F_Mix_s = 0.0 

do iz = 1, ntot
  F_Mix_s = F_Mix_s + avsol(iz,cc)*(log(avsol(iz,cc))-1.0)
  F_Mix_s = F_Mix_s - xsolbulk*(log(xsolbulk)-1.0)
enddo      

F_Mix_s = F_Mix_s * delta/vsol
F_tot = F_tot + F_Mix_s

! 6. Entropia interna polimero

F_Conf = 0.0

do i = 1, newcuantas
   F_Conf = F_Conf + (pro(i,cc))*log((pro(i,cc)))*2.0*sigma ! it is 2.0*sigma because wa have brushes on both walls
enddo

F_tot = F_tot + F_Conf 

! 8.vdW ! Ojo, los kai son negativos => atraccion

F_vdW = 0.0

do iz = 1, ntot
 do j = -Xulimit, Xulimit
   F_vdW = F_vdW - 0.5000*delta*xtotal(iz)*xtotal(iz+j)*Xu(j)*st/(vpol*vsol)/(vpol*vsol)
 enddo
enddo

F_tot = F_tot + F_vdW

! 10. Pol-Sup
!      F_eps = 0.0 
!      do iz = 1, dimz
!      F_eps = F_eps - eps(iz)
!      enddo
!      Free_Energy = Free_Energy + F_eps


print*, 'fe: Free energy 1:', F_tot, cc

! Segun Pong

F_tot2 = 0.0

sumrho=0.0
sumpi = 0.0

do iz=1,ntot
            
  sumpi = sumpi+log(avsol(iz,cc))     
  sumpi = sumpi-log(xsolbulk)     

  sumrho = sumrho + ( - avsol(iz,cc) )! sum over  rho_i i=+,-,si
  sumrho = sumrho - ( - xsolbulk )! sum over  rho_i i=+,-,si

enddo
         
sumpi = (delta/vsol)*sumpi
sumrho = (delta/vsol)*sumrho

F_tot2 = -2.0*sigma*log(q(cc)/shift) + sumpi + sumrho -F_vdW  ! It is 2.0*sigma because we have brush on both walls

print*, 'fe: Free energy 2:', F_tot2,cc

! Calcula mupol

mupol = -log(q(cc)/shift)
 

open(unit=20, file='F_tot.dat', access='append')
write(20,*)cc,F_tot
close(20)

open(unit=20, file='F_tot2.dat', access='append')
write(20,*)cc,F_tot2
close(20)

open(unit=20, file='F_mixs.dat', access='append')
write(20,*)cc,F_mix_s
close(20)

open(unit=20, file='mupol.dat', access='append')
write(20,*)cc,mupol
close(20)

open(unit=20, file='F_vdW.dat', access='append')
write(20,*)cc,F_vdW
close(20)

open(unit=20, file='F_conf.dat', access='append')
write(20,*)cc,F_conf
close(20)

end subroutine

