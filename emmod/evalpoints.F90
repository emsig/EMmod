SUBROUTINE evalpoints (xposnew,xtotrad1new,xtotrad2new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,xtotrad2,nlogx,&
	c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for pchip interpolation and Hankel transformations with two terms.
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! xtotrad2new: After completion of this function, this array contains the values of xtotrad2 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2:    2nd Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew), xtotrad2new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad1, imagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad1, dimagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad2, imagxtotrad2
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad2, dimagxtotrad2
REAL                  :: realtest1, imagtest1, realtest2, imagtest2
REAL                  :: lobound, hibound, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
	if (domark .eqv. .true.) then
		marker(ix) = 1; 
	else
		marker(ix) = 0;
	endif
	domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! compute the 1st derivative for the fixed points
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
! interpolate and test the points in question
do ix=1,nlogx
	if (marker(ix) == 0) then
		lopos = ix-1;
		hipos = ix+1;
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad1(lopos),realxtotrad1(hipos),&
			drealxtotrad1(lopos),drealxtotrad1(hipos),1,xpos(ix),realtest1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad1(lopos),imagxtotrad1(hipos),&
			dimagxtotrad1(lopos),dimagxtotrad1(hipos),1,xpos(ix),imagtest1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad2(lopos),realxtotrad2(hipos),&
			drealxtotrad2(lopos),drealxtotrad2(hipos),1,xpos(ix),realtest2,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad2(lopos),imagxtotrad2(hipos),&
			dimagxtotrad2(lopos),dimagxtotrad2(hipos),1,xpos(ix),imagtest2,next,ierr);
		lim = c1*log10(abs(realxtotrad1(ix)))+c2;
		lobound = log10(abs(realxtotrad1(ix)))-lim*abs(log10(abs(realxtotrad1(ix))));
		hibound = log10(abs(realxtotrad1(ix)))+lim*abs(log10(abs(realxtotrad1(ix))));
		if (real(log10(abs(realtest1)))>=real(lobound) .AND. real(log10(abs(realtest1)))<=real(hibound)) then
			marker(ix) = 1;
		else 
			marker(ix) = 2;
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(imagxtotrad1(ix)))+c2;
			lobound = log10(abs(imagxtotrad1(ix)))-lim*abs(log10(abs(imagxtotrad1(ix))));
			hibound = log10(abs(imagxtotrad1(ix)))+lim*abs(log10(abs(imagxtotrad1(ix))));
			if (real(log10(abs(imagtest1)))>=real(lobound) .AND. real(log10(abs(imagtest1)))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(realxtotrad2(ix)))+c2;
			lobound = log10(abs(realxtotrad2(ix)))-lim*abs(log10(abs(realxtotrad2(ix))));
			hibound = log10(abs(realxtotrad2(ix)))+lim*abs(log10(abs(realxtotrad2(ix))));
			if (real(log10(abs(realtest2)))>=real(lobound) .AND. real(log10(abs(realtest2)))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(imagxtotrad2(ix)))+c2;
			lobound = log10(abs(imagxtotrad2(ix)))-lim*abs(log10(abs(imagxtotrad2(ix))));
			hibound = log10(abs(imagxtotrad2(ix)))+lim*abs(log10(abs(imagxtotrad2(ix))));
			if (real(log10(abs(imagtest2)))>=real(lobound) .AND. real(log10(abs(imagtest2)))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
	endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
	if (marker(ix) == 1) then
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		xtotrad2new(ixs) = xtotrad2(ix);
		ixs = ixs + 1;
	else if (marker(ix) == 2) then
		newit = 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
		xtotrad1new(ixs) = 0.0;
		xtotrad2new(ixs) = 0.0;
		ixs = ixs + 1;
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		xtotrad2new(ixs) = xtotrad2(ix);
		ixs = ixs + 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
		xtotrad1new(ixs) = 0.0;
		xtotrad2new(ixs) = 0.0;
		ixs = ixs + 1;
	endif
end do
nlogxdata = ixs-1;
END SUBROUTINE evalpoints

SUBROUTINE evalpoints_mono (xposnew,xtotrad1new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,nlogx,&
	c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for pchip interpolation and Hankel transformations with one term.
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
REAL, DIMENSION(:), ALLOCATABLE    :: realxtotrad1, imagxtotrad1
REAL, DIMENSION(:), ALLOCATABLE    :: drealxtotrad1, dimagxtotrad1
REAL                  :: realtest1, imagtest1, realtest2, imagtest2
REAL                  :: lobound, hibound, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
	if (domark .eqv. .true.) then
		marker(ix) = 1; 
	else
		marker(ix) = 0;
	endif
	domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! compute the 1st derivative for the fixed points
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
! interpolate and test the points in question
do ix=1,nlogx
	if (marker(ix) == 0) then
		lopos = ix-1;
		hipos = ix+1;
		CALL dchfev(xpos(lopos),xpos(hipos),realxtotrad1(lopos),realxtotrad1(hipos),&
			drealxtotrad1(lopos),drealxtotrad1(hipos),1,xpos(ix),realtest1,next,ierr);
		CALL dchfev(xpos(lopos),xpos(hipos),imagxtotrad1(lopos),imagxtotrad1(hipos),&
			dimagxtotrad1(lopos),dimagxtotrad1(hipos),1,xpos(ix),imagtest1,next,ierr);
		lim = c1*log10(abs(realxtotrad1(ix)))+c2;
		lobound = log10(abs(realxtotrad1(ix)))-lim*abs(log10(abs(realxtotrad1(ix))));
		hibound = log10(abs(realxtotrad1(ix)))+lim*abs(log10(abs(realxtotrad1(ix))));
		if (real(log10(abs(realtest1)))>=real(lobound) .AND. real(log10(abs(realtest1)))<=real(hibound)) then
			marker(ix) = 1;
		else 
			marker(ix) = 2;
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(imagxtotrad1(ix)))+c2;
			lobound = log10(abs(imagxtotrad1(ix)))-lim*abs(log10(abs(imagxtotrad1(ix))));
			hibound = log10(abs(imagxtotrad1(ix)))+lim*abs(log10(abs(imagxtotrad1(ix))));
			if (real(log10(abs(imagtest1)))>=real(lobound) .AND. real(log10(abs(imagtest1)))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
	endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
	if (marker(ix) == 1) then
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		ixs = ixs + 1;
	else if (marker(ix) == 2) then
		newit = 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
		xtotrad1new(ixs) = 0.0;
		ixs = ixs + 1;
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		ixs = ixs + 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
		xtotrad1new(ixs) = 0.0;
		ixs = ixs + 1;
	endif
end do
nlogxdata = ixs-1;
END SUBROUTINE evalpoints_mono

SUBROUTINE evalpoints_lin (xposnew,xtotrad1new,xtotrad2new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,xtotrad2,nlogx,&
	c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for linear interpolation and Hankel transformations with two terms.
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! xtotrad2new: After completion of this function, this array contains the values of xtotrad2 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! xtotrad2:    2nd Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew), xtotrad2new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx), xtotrad2(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
COMPLEX               :: test1, test2, lobound, hibound
REAL                  :: thisdx, weightlo, weighthi, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
	if (domark .eqv. .true.) then
		marker(ix) = 1; 
	else
		marker(ix) = 0;
	endif
	domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! interpolate and test the points in question
do ix=1,nlogx
	if (marker(ix) == 0) then
		lopos = ix-1;
		hipos = ix+1;
		thisdx = xpos(hipos) - xpos(lopos);
		weightlo = (xpos(hipos)-xpos(ix))/thisdx;
		weighthi = (xpos(ix)-xpos(lopos))/thisdx;
		test1 = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
		lim = c1*log10(abs(real(xtotrad1(ix))))+c2;
		lobound = log10(abs(real(xtotrad1(ix))))-lim*abs(log10(abs(real(xtotrad1(ix)))));
		hibound = log10(abs(real(xtotrad1(ix))))+lim*abs(log10(abs(real(xtotrad1(ix)))));
		if (real(log10(abs(real(test1))))>=real(lobound) .AND. real(log10(abs(real(test1))))<=real(hibound)) then
			marker(ix) = 1;
		else 
			marker(ix) = 2;
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(aimag(xtotrad1(ix))))+c2;
			lobound = log10(abs(aimag(xtotrad1(ix))))-lim*abs(log10(abs(aimag(xtotrad1(ix)))));
			hibound = log10(abs(aimag(xtotrad1(ix))))+lim*abs(log10(abs(aimag(xtotrad1(ix)))));
			if (real(log10(abs(aimag(test1))))>=real(lobound) .AND. real(log10(abs(aimag(test1))))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
		if (marker(ix) == 1) then
			test2 = weightlo*xtotrad2(lopos)+weighthi*xtotrad2(hipos);
			lim = c1*log10(abs(real(xtotrad2(ix))))+c2;
			lobound = log10(abs(real(xtotrad2(ix))))-lim*abs(log10(abs(real(xtotrad2(ix)))));
			hibound = log10(abs(real(xtotrad2(ix))))+lim*abs(log10(abs(real(xtotrad2(ix)))));
			if (real(log10(abs(real(test2))))>=real(lobound) .AND. real(log10(abs(real(test2))))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
		if (marker(ix) == 1) then
			lim = c1*log10(abs(aimag(xtotrad2(ix))))+c2;
			lobound = log10(abs(aimag(xtotrad2(ix))))-lim*abs(log10(abs(aimag(xtotrad2(ix)))));
			hibound = log10(abs(aimag(xtotrad2(ix))))+lim*abs(log10(abs(aimag(xtotrad2(ix)))));
			if (real(log10(abs(aimag(test2))))>=real(lobound) .AND. real(log10(abs(aimag(test2))))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
	endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
	if (marker(ix) == 1) then
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		xtotrad2new(ixs) = xtotrad2(ix);
		ixs = ixs + 1;
	else if (marker(ix) == 2) then
		newit = 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
		xtotrad1new(ixs) = 0.0;
		xtotrad2new(ixs) = 0.0;
		ixs = ixs + 1;
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		xtotrad2new(ixs) = xtotrad2(ix);
		ixs = ixs + 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
		xtotrad1new(ixs) = 0.0;
		xtotrad2new(ixs) = 0.0;
		ixs = ixs + 1;
	endif
end do
nlogxdata = ixs-1;
END SUBROUTINE evalpoints_lin

SUBROUTINE evalpoints_lin_mono (xposnew,xtotrad1new,markernew,nlogxnew,nlogxdata,newit,xpos,xtotrad1,nlogx,&
	c1,c2)
! Evaluate which points are necessary to achieve the precision defined with the parameters c1 and c2 
! for linear interpolation and Hankel transformations with one term.
!
! Calling arguments:
! xposnew:     After completion of the function, this array contains the coordinates in the space domain 
!              for all datapoints for which xtotrad1 and xtotrad2 have been computed 
!              as well as all the possible new datapoints, real array of size nlogxnew 
! xtotrad1new: After completion of this function, this array contains the values of xtotrad1 (Hankel transformed EM-field) 
!              and zeros on the location where the Hankel transformation still needs to be computed, complex array of size nlogxnew
! markernew:   After completion of this function, this array hold the indication for which locations in space 
!              the Hankel transformation still needs to be computed, integer array of size nlogxnew
! nlogxnew:    Maximum amount of values in the space domain that are possible. which is  
!              nlogx points of data plus nlogx-1 points that could be added, integer scalar
!              Note: In the worst case, between every point a new datapoint is required.
! nlogxdata:   After completion of this function, this parameter contains the amount of actual datapoints, which is
!              nlogx points of data plus the amount of points that were actually added, integer scalar
!              Note: This value can not be larger than nlogxnew.
! newit:       After completion of this function, this parameter indicates if a new iteration is required, i.e., 
!              if new datapoints have to be added, integer scalar
! xpos:        Logarithmic coordinate vector in the space domain, real array of size nlogx
! xtotrad1:    1st Hankel transformed field-term, complex array of size nlogx
! nlogx:       Amount of samples in logarithmic coordinate vector, integer scalar
! c1:          1st precision parameter, real scalar
! c2:          2nd precision parameter, real scalar

IMPLICIT NONE

COMPLEX,INTENT(INOUT) :: xtotrad1new(nlogxnew)
REAL,INTENT(INOUT)    :: xposnew(nlogxnew)
INTEGER,INTENT(INOUT) :: markernew(nlogxnew), newit, nlogxdata
COMPLEX,INTENT(IN)    :: xtotrad1(nlogx)
REAL   ,INTENT(IN)    :: xpos(nlogx), c1, c2
INTEGER,INTENT(IN)    :: nlogx, nlogxnew
COMPLEX               :: test1, test2, lobound, hibound
REAL                  :: thisdx, weightlo, weighthi, lim
INTEGER, DIMENSION(:), ALLOCATABLE :: marker
INTEGER               :: nlogxfix, nlogxint, ix, ixs, lopos, hipos, ierr, next
LOGICAL               :: domark

! mark elements of xpos with 1 that are kept and  with 0 that are tested
domark = .true.;
allocate(marker(1:nlogx))
do ix=1,nlogx
	if (domark .eqv. .true.) then
		marker(ix) = 1; 
	else
		marker(ix) = 0;
	endif
	domark = .not.domark; 
end do
! the last element of xpos should always be kept
marker(nlogx) = 1;
! count the elements that are kept
nlogxfix = sum(marker);
! compute the elements that need to be tested
nlogxint = nlogx-nlogxfix;
! interpolate and test the points in question
do ix=1,nlogx
	if (marker(ix) == 0) then
		lopos = ix-1;
		hipos = ix+1;
		thisdx = xpos(hipos) - xpos(lopos);
		weightlo = (xpos(hipos)-xpos(ix))/thisdx;
		weighthi = (xpos(ix)-xpos(lopos))/thisdx;
		test1 = weightlo*xtotrad1(lopos)+weighthi*xtotrad1(hipos);
		lim = c1*log10(abs(real(xtotrad1(ix))))+c2;
		lobound = log10(abs(real(xtotrad1(ix))))-lim*abs(log10(abs(real(xtotrad1(ix)))));
		hibound = log10(abs(real(xtotrad1(ix))))+lim*abs(log10(abs(real(xtotrad1(ix)))));
		if (real(log10(abs(real(test1))))>=real(lobound) .AND. real(log10(abs(real(test1))))<=real(hibound)) then
			marker(ix) = 1;
		else
			marker(ix) = 2;
		endif
		! Testing of the imaginary part is only required if no datapoints need to be added for the real part.
		if (marker(ix) == 1) then
			lim = c1*log10(abs(aimag(xtotrad1(ix))))+c2;
			lobound = log10(abs(aimag(xtotrad1(ix))))-lim*abs(log10(abs(aimag(xtotrad1(ix)))));
			hibound = log10(abs(aimag(xtotrad1(ix))))+lim*abs(log10(abs(aimag(xtotrad1(ix)))));
			if (real(log10(abs(aimag(test1))))>=real(lobound) .AND. real(log10(abs(aimag(test1))))<=real(hibound)) then
				marker(ix) = 1;
			else 
				marker(ix) = 2;
			endif
		endif
	endif
end do
! Set up the new coordinate and data vector
ixs = 1;
newit = 0;
do ix = 1,nlogx
	if (marker(ix) == 1) then
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		ixs = ixs + 1;
	else if (marker(ix) == 2) then
		newit = 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix)-xpos(ix-1))/2+xpos(ix-1);
		xtotrad1new(ixs) = 0.0;
		ixs = ixs + 1;
		markernew(ixs) = 1;
		xposnew(ixs) = xpos(ix);
		xtotrad1new(ixs) = xtotrad1(ix);
		ixs = ixs + 1;
		markernew(ixs) = 2;
		xposnew(ixs) = (xpos(ix+1)-xpos(ix))/2+xpos(ix);
		xtotrad1new(ixs) = 0.0;
		ixs = ixs + 1;
	endif
end do
nlogxdata = ixs-1;
END SUBROUTINE evalpoints_lin_mono
