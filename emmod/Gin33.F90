SUBROUTINE gin33 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, xpos, nlogx);

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), zrcv, zsrc
COMPLEX               :: Rbig, G, term1, term2
REAL                  :: pi
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
	Rbig = sqrt(xpos(ix)**2.0+etaH*(zrcv-zsrc)**2.0/etaV);
	if (abs(Rbig) <= 10**(-6.0)) then
		Rbig = 10**(-6.0);
	endif
	G = exp(-1.0*sqrt(gamA)*Rbig)/(4.0*pi*Rbig*sqrt(etaH/etaV));
	term1 = (((zrcv-zsrc)*etaH/etaV)**2.0*(sqrt(gamA)/Rbig+2.0/(Rbig**2.0))-etaH/etaV)*(sqrt(gamA)/Rbig+Rbig**(-2.0));
	term2 = ((zrcv-zsrc)*etaH/(etaV*Rbig**2.0))**2.0-etaH*zetaH;
	Gin(ix) = (term1 + term2)*G/etaV;
end do

END SUBROUTINE gin33
