SUBROUTINE gridit_ref (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx)
! Compute the data of one quadrant of the grid based on the radial data using pchip interpolation
! for component 77 and 88 (TM-mode and TE-mode reflection response, respectively).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the reflection response in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx
COMPLEX               :: i
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx, realtemp1, imagtemp1, realtemp2, imagtemp2
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad1, imagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad1, dimagxtotrad1
INTEGER               :: ix, iy, lopos, hipos, minel, ierr, next
INTEGER, DIMENSION(1) :: temp

i=cmplx(0.,1.)
allocate(realxtotrad1(1:nlogx))
allocate(imagxtotrad1(1:nlogx))
do ix=1,nlogx
	realxtotrad1(ix) = real(xtotrad1(ix));
	imagxtotrad1(ix) = aimag(xtotrad1(ix));
end do
allocate(drealxtotrad1(1:nlogx))
allocate(dimagxtotrad1(1:nlogx))
CALL dpchim (nlogx,xpos,realxtotrad1,drealxtotrad1,1,ierr);
CALL dpchim (nlogx,xpos,imagxtotrad1,dimagxtotrad1,1,ierr);
do iy=1,nyh
	y = (iy-1)*dy;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		temp = minloc(abs(xpos-rad));
		minel = temp(1);
		if (xpos(minel)<=rad) then
			lopos = minel;
			hipos = minel + 1;
		else
			hipos = minel;
			lopos = minel - 1;
		endif
		if (lopos == hipos) then
			hipos = hipos + 1;
		endif
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad1(lopos),realxtotrad1(hipos),&
			drealxtotrad1(lopos),drealxtotrad1(hipos),1,rad,realtemp1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad1(lopos),imagxtotrad1(hipos),&
			dimagxtotrad1(lopos),dimagxtotrad1(hipos),1,rad,imagtemp1,next,ierr);
		Ptot(ix,iy) = realtemp1+i*imagtemp1;
	end do
end do
END SUBROUTINE gridit_ref

SUBROUTINE gridit_ref_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 77 and 88 (TM-mode and TE-mode reflection response, respectively).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the reflection response in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, lopos, hipos, minel
INTEGER, DIMENSION(1) :: temp

do iy=1,nyh
	y = (iy-1)*dy;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		temp = minloc(abs(xpos-rad));
		minel = temp(1);
		if (xpos(minel)<=rad) then
			lopos = minel;
			hipos = minel + 1;
		else
			hipos = minel;
			lopos = minel - 1;
		endif
		if (lopos == hipos) then
			hipos = hipos + 1;
		endif
		thisdx = xpos(hipos) - xpos(lopos);
		weightlo = (xpos(hipos)-rad)/thisdx;
		weighthi = (rad-xpos(lopos))/thisdx;
		Ptot(ix,iy) = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
	end do
end do
END SUBROUTINE gridit_ref_lin
