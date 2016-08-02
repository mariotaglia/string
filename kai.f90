!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 1D-coordenadas
! polares usando un metodo de MC
!
!#####################################################################


subroutine kai

use kai
use layer
use brush

implicit none
integer seed
real*8 xmin,xmax,ymin,ymax,zmin,zmax
integer MCsteps ! numero de steps de MC
integer jx,jy
real*8 R,theta,z
real*8 rn
integer i, ii
real*8 rands
real*8 pi
real*8 x1,x2,y1, y2, z1, z2, vect
integer iR, ix,iy,iz, itheta
integer j
real*8 radio
real*8 cutoff
integer rank
rank = 0

cutoff = (float(Xulimit)+0.5)*delta


if(rank.eq.0)print*,'Kai calculation'
open(unit=111, file='kais.dat')

pi=dacos(-1.0d0)          ! pi = arccos(-1) 

Xu = 0.0 ! vector Xu

seed = 1010
MCsteps = 200

      ymax =cutoff
      ymin = -cutoff

      zmax = cutoff
      zmin = -cutoff

      xmax = cutoff
      xmin = -cutoff
    
      do ix = 1, MCsteps
      do iy = 1, MCsteps
      do iz = 1, MCsteps

! coordenadas del segmento (x1,y1,z1) y del punto a integrar (x2,y2,z2)

         x1 = 0.0 ! position of first segment
         y1 = 0.0
         z1 = 0.0

         x2 = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(MCsteps-1)
         y2 = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(MCsteps-1)
         z2 = zmin + (zmax-zmin)*dfloat(iz-1)/dfloat(MCsteps-1)

         R = x2

         vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diferencia
         jx = int(x2/delta) ! j tiene la celda donde cae el punto a integrar
         jy = int(y2/delta) ! j tiene la celda donde cae el punto a integrar

         

!         if((j.gt.Xulimit).or.(j.lt.-Xulimit)) then
!            print*,'kai: error', j, R, x2
!         endif

         if(vect.le.(cutoff)) then ! esta dentro de la esfera del cut-off   
         if(vect.ge.lseg) then ! esta dentro de la esfera del segmento
         Xu(jx,jy) = Xu(jx,jy) + ((lseg/vect)**6) ! incluye el jacobiano R(segmento)
         endif
         endif

      enddo
      enddo
      enddo

      do j = -Xulimit, Xulimit
      Xu(jx,jy) = Xu(jx,jy)/(MCsteps**3)*(2.0*cutoff)**3
     write(111,*)jx,jy,Xu(jx,jy) ! residual size of iteration vector
     enddo

close(111)
end

