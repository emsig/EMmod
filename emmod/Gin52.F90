SUBROUTINE gin52 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, x, y, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: RbigA, RbigB, mult, afac, bfac, term1, term2
REAL                  :: pi, rad
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
	if (abs(xpos(ix))<= 10**(-6.0)) then
		rad = 10**(-6.0);
	else
		rad = xpos(ix);
	endif
	RbigA = sqrt(rad**2.0+(etaH*(zrcv-zsrc)**2.0/etaV));
	RbigB = sqrt(rad**2.0+(zetaH*(zrcv-zsrc)**2.0/zetaV));

	mult = y(ix)*(zrcv-zsrc)/(4.0*pi*sqrt(etaH*zetaH));
	afac = sqrt(gamA)*etaH/(RbigA*etaV);
	bfac = sqrt(gamB)*zetaH/(RbigB*zetaV);
	term1 = afac*mult*exp(-1.0*sqrt(gamA)*RbigA)*(2.0*x(ix)*rad**(-4.0)+sqrt(gamA)*x(ix)/(rad**2.0*RbigA)+x(ix)/(rad**2.0*RbigA**2.0));
	term2 = bfac*mult*exp(-1.0*sqrt(gamB)*RbigB)*(2.0*x(ix)*rad**(-4.0)+sqrt(gamB)*x(ix)/(rad**2.0*RbigB)+x(ix)/(rad**2.0*RbigB**2.0));
	Gin(ix) = -1.0*term1 + term2;
end do

END SUBROUTINE gin52
