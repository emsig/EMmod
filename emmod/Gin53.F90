SUBROUTINE gin53 (Gin, zrcv, zsrc, etaH, etaV, zetaH, gamA, xpos, x, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), zrcv, zsrc
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
	Gin(ix) = etaH/etaV*(sqrt(gamA)*x(ix)*RbigA**(-1)+x(ix)*RbigA**(-2))*GA;
end do

END SUBROUTINE gin53
