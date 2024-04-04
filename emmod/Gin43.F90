SUBROUTINE gin43 (Gin, zrcv, zsrc, etaH, etaV, zetaH, gamA, xpos, y, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: RbigA, GA 
REAL                  :: pi, rad, d, signd
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
	if (abs(xpos(ix))<= 10**(-6.0)) then
		rad = 10**(-6.0);
	else
		rad = xpos(ix);
	endif
	RbigA = sqrt(rad**2.0+(etaH*(zrcv-zsrc)**2.0/etaV));

	GA = exp(-1.0*sqrt(gamA)*RbigA)/(4.0*pi*RbigA*sqrt(etaH/etaV));
	Gin(ix) = -1.0*etaH/etaV*(sqrt(gamA)*y(ix)*RbigA**(-1)+y(ix)*RbigA**(-2))*GA;
end do

END SUBROUTINE gin43
