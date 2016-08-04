function imap(ix,iy)
use brush
implicit none
integer ix,iy,imap
imap = ix+dimx*(iy-1)
end function

function mapx(i)
use brush
implicit none
integer mapx, i
mapx = mod((i-1),dimx)+1
endfunction 

function mapy(i)
use brush
implicit none
integer mapy, i
mapy = int((i-1)/dimx)+1
endfunction

