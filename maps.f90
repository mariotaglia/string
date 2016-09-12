subroutine makemaps
use brush
use maps
implicit none
integer ix,iy,i

do ix = 1, dimx
do iy = 1, dimy
imap(ix,iy) = ix+dimx*(iy-1)
enddo
enddo

do i = 1, dimx*dimy
mapx(i) = mod((i-1),dimx)+1
mapy(i) = int((i-1)/dimx)+1
enddo

end

