SUBROUTINE gridit_yy (Ptot,xtotrad1,xtotrad2,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using pchip interpolation
! for component 22 (crossline oriented electric receiver, crossline oriented electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2: 2nd Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! gamB:     Small gamma bar squared (zetaV*etaH), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin, i
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx, realtemp1, imagtemp1, realtemp2, imagtemp2
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad1, imagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad1, dimagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad2, imagxtotrad2
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad2, dimagxtotrad2
INTEGER               :: ix, iy, ipos, lopos, hipos, ierr, next

i=cmplx(0.,1.)
allocate(realxtotrad1(1:nlogx))
allocate(imagxtotrad1(1:nlogx))
allocate(realxtotrad2(1:nlogx))
allocate(imagxtotrad2(1:nlogx))
do ix=1,nlogx
	realxtotrad1(ix) = real(xtotrad1(ix));
	imagxtotrad1(ix) = aimag(xtotrad1(ix));
	realxtotrad2(ix) = real(xtotrad2(ix));
	imagxtotrad2(ix) = aimag(xtotrad2(ix));
end do
allocate(drealxtotrad1(1:nlogx))
allocate(dimagxtotrad1(1:nlogx))
allocate(drealxtotrad2(1:nlogx))
allocate(dimagxtotrad2(1:nlogx))
CALL dpchim (nlogx,xpos,realxtotrad1,drealxtotrad1,1,ierr);
CALL dpchim (nlogx,xpos,imagxtotrad1,dimagxtotrad1,1,ierr);
CALL dpchim (nlogx,xpos,realxtotrad2,drealxtotrad2,1,ierr);
CALL dpchim (nlogx,xpos,imagxtotrad2,dimagxtotrad2,1,ierr);
do iy=1,nyh
	y = (iy-1)*dy;
	lopos = 1;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		do ipos=lopos,nlogx
			if (xpos(ipos) > rad .OR. ipos == nlogx) then
				hipos = ipos;
				lopos = ipos -1;
				exit
			endif
		end do
		if (ix==1) then
			fact = -1.0;
		else
			fact = cos(2.0*atan(y/x));
		endif
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad1(lopos),realxtotrad1(hipos),&
			drealxtotrad1(lopos),drealxtotrad1(hipos),1,rad,realtemp1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad1(lopos),imagxtotrad1(hipos),&
			dimagxtotrad1(lopos),dimagxtotrad1(hipos),1,rad,imagtemp1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad2(lopos),realxtotrad2(hipos),&
			drealxtotrad2(lopos),drealxtotrad2(hipos),1,rad,realtemp2,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad2(lopos),imagxtotrad2(hipos),&
			dimagxtotrad2(lopos),dimagxtotrad2(hipos),1,rad,imagtemp2,next,ierr);
		Ptot(ix,iy) = (realtemp1+i*imagtemp1) + fact*(realtemp2+i*imagtemp2);
		if (above==0 .AND. xdirect==1) then
			CALL gin22(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,rad,y,1);
			Ptot(ix,iy) = Ptot(ix,iy) + gin;
		endif
	end do
end do
END SUBROUTINE gridit_yy
	  
SUBROUTINE gridit_yy_lin (Ptot,xtotrad1,xtotrad2,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 22 (crossline oriented electric receiver, crossline oriented electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2: 2nd Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! gamB:     Small gamma bar squared (zetaV*etaH), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos

do iy=1,nyh
	y = (iy-1)*dy;
	lopos = 1;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		do ipos=lopos,nlogx
			if (xpos(ipos) > rad .OR. ipos == nlogx) then
				hipos = ipos;
				lopos = ipos -1;
				exit
			endif
		end do
		thisdx = xpos(hipos) - xpos(lopos);
		weightlo = (xpos(hipos)-rad)/thisdx;
		weighthi = (rad-xpos(lopos))/thisdx;
		if (ix==1) then
			fact = -1.0;
		else
			fact = cos(2.0*atan(y/x));
		endif
		Ptot(ix,iy) = (weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
		Ptot(ix,iy) = Ptot(ix,iy) + fact*(weightlo*xtotrad2(lopos)+weighthi*xtotrad2(hipos));
		if (above==0 .AND. xdirect==1) then
			CALL gin22(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,rad,y,1);
			Ptot(ix,iy) = Ptot(ix,iy) + gin;
		endif
	end do
end do
END SUBROUTINE gridit_yy_lin

SUBROUTINE gridit_yy_fullspace (Ptot,dx,nxh,dy,nyh,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB)
! Compute the data of one quadrant of the grid for a homogeneous fullspace
! for component 22 (crossline oriented electric receiver, crossline oriented electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! gamB:     Small gamma bar squared (zetaV*etaH), complex scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh
COMPLEX               :: gin
REAL                  :: x, y, rad
INTEGER               :: ix, iy

do iy=1,nyh
	y = (iy-1)*dy;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		CALL gin22(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,gamB,rad,y,1);
		Ptot(ix,iy) = Ptot(ix,iy) + gin;
	end do
end do
END SUBROUTINE gridit_yy_fullspace

SUBROUTINE gridit_yz (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using pchip interpolation
! for component 23 (crossline oriented electric receiver, vertical electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin, i
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx, realtemp1, imagtemp1, realtemp2, imagtemp2
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad1, imagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad1, dimagxtotrad1
INTEGER               :: ix, iy, ipos, lopos, hipos, ierr, next

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
	lopos = 1;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		do ipos=lopos,nlogx
			if (xpos(ipos) > rad .OR. ipos == nlogx) then
				hipos = ipos;
				lopos = ipos -1;
				exit
			endif
		end do
		if (ix==1) then
			fact = 1.0;
		else
			fact = sin(atan(y/x));
		endif
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad1(lopos),realxtotrad1(hipos),&
			drealxtotrad1(lopos),drealxtotrad1(hipos),1,rad,realtemp1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad1(lopos),imagxtotrad1(hipos),&
			dimagxtotrad1(lopos),dimagxtotrad1(hipos),1,rad,imagtemp1,next,ierr);
		Ptot(ix,iy) = fact*(realtemp1+i*imagtemp1);
		if (above==0 .AND. xdirect==1) then
			CALL gin23(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,y,1);
			Ptot(ix,iy) = Ptot(ix,iy) - gin;
		endif
	end do
end do
END SUBROUTINE gridit_yz

SUBROUTINE gridit_yz_lin (Ptot,xtotrad1,dx,nxh,dy,nyh,xpos,nlogx,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,above,xdirect)
! Compute the data of one quadrant of the grid based on the radial data using linear interpolation
! for component 23 (crossline oriented electric receiver, vertical electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! xtotrad1: 1st Hankel transformed field-term, complex array of size nlogx
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh, nlogx, above, xdirect
COMPLEX               :: gin
REAL                  :: x, y, rad, radfilo, radfihi, weightlo, weighthi, fact
REAL                  :: thisdx
INTEGER               :: ix, iy, ipos, lopos, hipos

do iy=1,nyh
	y = (iy-1)*dy;
	lopos = 1;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		do ipos=lopos,nlogx
			if (xpos(ipos) > rad .OR. ipos == nlogx) then
				hipos = ipos;
				lopos = ipos -1;
				exit
			endif
		end do
		thisdx = xpos(hipos) - xpos(lopos);
		weightlo = (xpos(hipos)-rad)/thisdx;
		weighthi = (rad-xpos(lopos))/thisdx;
		if (ix==1) then
			fact = 1.0;
		else
			fact = sin(atan(y/x));
		endif
		Ptot(ix,iy) = fact*(weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos));
		if (above==0 .AND. xdirect==1) then
			CALL gin23(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,y,1);
			Ptot(ix,iy) = Ptot(ix,iy) - gin;
		endif
	end do
end do
END SUBROUTINE gridit_yz_lin

SUBROUTINE gridit_yz_fullspace (Ptot,dx,nxh,dy,nyh,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA)
! Compute the data of one quadrant of the grid for a homogeneous fullspace
! for component 23 (crossline oriented electric receiver, vertical electric source).
! 
! Calling arguments:
! Ptot:     After completion of this function, this array will contain one quadrant of the EM-field in the space domain, 
!           complex array of size nxh times nyh
! dx:       Sampling in inline direction (x-direction), real scalar
! nxh:      Number of samples of one quadrant in inline direction (x-direction), integer scalar
! dy:       Sampling in crossline direction (y-direction), real scalar
! nyh:      Number of samples of one quadrant in crossline direction (y-direction), integer scalar
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! etaH:     Material parameter eta for the horizontal direction, complex scalar
! etaV:     Material parameter eta for the vertical direction, complex scalar
! zetaH:    Material parameter zeta for the horizontal direction, complex scalar
! zetaV:    Material parameter zeta for the vertical direction, complex scalar
! gamA:     Small gamma squared (zetaH*etaV), complex scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
COMPLEX,INTENT(IN)    :: etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: dx, dy, zrcv, zsrc
INTEGER,INTENT(IN)    :: nxh, nyh
COMPLEX               :: gin
REAL                  :: x, y, rad
INTEGER               :: ix, iy

do iy=1,nyh
	y = (iy-1)*dy;
	do ix=1,nxh
		x = (ix-1)*dx;
		rad = sqrt(x*x+y*y);
		CALL gin23(gin,zrcv,zsrc,etaH,etaV,zetaH,zetaV,gamA,rad,y,1);
		Ptot(ix,iy) = Ptot(ix,iy) - gin;
	end do
end do
END SUBROUTINE gridit_yz_fullspace
