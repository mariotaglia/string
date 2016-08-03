subroutine free_energy
use brush
use layer
use volume
use bulk
use kai
use partfunc
implicit none
integer ii
real*8 F_tot, F_tot2
real*8 F_Mix_s, F_Conf, F_vdW
integer iz, ix, iy, jx, jy, kx,ky, i, j, k
real*8 sumpi, sumrho
real*8 xtotal(1-Xulimit:ntot+Xulimit)
real*8 mupol
integer, external :: imap, mapx, mapy

xtotal = 0.0
do i = 1,ntot
xtotal(i) = 1.0-avsol(i)
enddo


F_tot = 0.0
F_tot2 = 0.0

! 1. Mezcla solvente

F_Mix_s = 0.0 

do iz = 1, ntot
  F_Mix_s = F_Mix_s + avsol(iz)*(dlog(avsol(iz))-1.0)
  F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)
enddo      

F_Mix_s = F_Mix_s * delta/vsol
F_tot = F_tot + F_Mix_s

! 6. Entropia interna polimero

F_Conf = 0.0

do ii = 1, dimx
do i = 1, newcuantas
   F_Conf = F_Conf + (pro(i,ii)/q(ii))*dlog((pro(i,ii))/q(ii))*2.0*sigma ! it is 2.0*sigma because wa have brushes on both walls
enddo
enddo

F_tot = F_tot + F_Conf 

! 8.vdW ! Ojo, los kai son negativos => atraccion

F_vdW = 0.0

do i = 1, ntot
 ix=mapx(i)
 iy=mapy(i)
 do jx = -Xulimit, Xulimit
 do jy = -Xulimit, Xulimit
  kx = ix+jx
  kx= mod(kx-1+50*dimx, dimx) + 1
  ky = ix+jy
  ky= mod(ky-1+50*dimy, dimy) + 1
  k = imap(kx,ky)

  F_vdW = F_vdW - 0.5000*delta*xtotal(i)*xtotal(k)*Xu(jx,jy)*st/(vpol*vsol)/(vpol*vsol)

 enddo
 enddo
enddo

F_tot = F_tot + F_vdW

! 10. Pol-Sup
!      F_eps = 0.0 
!      do iz = 1, dimz
!      F_eps = F_eps - eps(iz)
!      enddo
!      Free_Energy = Free_Energy + F_eps


print*, 'fe: Free energy 1:', F_tot

! Segun Pong

F_tot2 = 0.0

sumrho=0.0
sumpi = 0.0

do iz=1,ntot
            
  sumpi = sumpi+dlog(avsol(iz))     
  sumpi = sumpi-dlog(xsolbulk)     

  sumrho = sumrho + ( - avsol(iz) )! sum over  rho_i i=+,-,si
  sumrho = sumrho - ( - xsolbulk )! sum over  rho_i i=+,-,si

enddo
         
sumpi = (delta/vsol)*sumpi
sumrho = (delta/vsol)*sumrho

F_tot2 = 0.0
do ii = 1, dimx
F_tot2 = F_tot2  -2.0*sigma*dlog(q(ii)/shift)
enddo

F_tot2 = F_tot2 + sumpi + sumrho -F_vdW  ! It is 2.0*sigma because we have brush on both walls

print*, 'fe: Free energy 2:', F_tot2

! Calcula mupol

mupol = -dlog(q(1)/shift)
 

open(unit=20, file='F_tot.dat', access='append')
write(20,*)st,F_tot
close(20)

open(unit=20, file='F_tot2.dat', access='append')
write(20,*)st,F_tot2
close(20)

open(unit=20, file='F_mixs.dat', access='append')
write(20,*)st,F_mix_s
close(20)

open(unit=20, file='mupol.dat', access='append')
write(20,*)st,mupol
close(20)

open(unit=20, file='F_vdW.dat', access='append')
write(20,*)st,F_vdW
close(20)

open(unit=20, file='F_conf.dat', access='append')
write(20,*)st,F_conf
close(20)

end subroutine

