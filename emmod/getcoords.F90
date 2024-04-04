SUBROUTINE getcoords (xpos,startlogx,deltalogx,nlogx)
! Set up logarithmic coordinate vector
!
! Calling arguments:
! xpos:      Logarithmic coordinate vector in the space domain, real array of size nlogx
! startlogx: Exponent of the smallest coordinate point larger than 0.0, real scalar
! deltalogx: Logarithmic spatial sampling, real scalar
! nlogx:     Amount of samples in logarithmic coordinate vector, integer scalar

IMPLICIT NONE

REAL,INTENT(INOUT)    :: xpos(nlogx)
REAL,INTENT(IN)       :: startlogx, deltalogx
INTEGER,INTENT(IN)    :: nlogx
INTEGER               :: n

xpos(1) = 0.0;
do n=2,nlogx
	xpos(n) = 10**(startlogx+(n-1)*deltalogx);
end do

END SUBROUTINE getcoords 
