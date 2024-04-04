SUBROUTINE gin61 (Gin, zrcv, zsrc, zetaH, zetaV, gamB, xpos, y, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: zetaH, zetaV, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: RbigB, GB 
REAL                  :: pi, rad, d, signd
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
	if (abs(xpos(ix))<= 10**(-6.0)) then
		rad = 10**(-6.0);
	else
		rad = xpos(ix);
	endif
	RbigB = sqrt(rad**2.0+(zetaH*(zrcv-zsrc)**2.0/zetaV));
	GB = exp(-1.0*sqrt(gamB)*RbigB)/(4.0*pi*RbigB*sqrt(zetaH/zetaV));
	Gin(ix) = zetaH/zetaV*(sqrt(gamB)*y(ix)*RbigB**(-1.0)+y(ix)*RbigB**(-2.0))*GB;
end do

END SUBROUTINE gin61
