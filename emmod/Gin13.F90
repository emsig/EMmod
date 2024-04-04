SUBROUTINE gin13 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, xpos, x, nlogx);

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), zrcv, zsrc
COMPLEX               :: Rbig, G
REAL                  :: pi, rad
INTEGER               :: ix

pi = 2*acos(0.0);

do ix=1,nlogx
        if (abs(xpos(ix))<= 10**(-6.0)) then
                rad = 10**(-6.0);
        else
                rad = xpos(ix);
        endif
        Rbig = sqrt(rad**2.0+(etaH*(zrcv-zsrc)**2.0/etaV));
	!Rbig = sqrt(xpos(ix)**2.0+etaH*(zrcv-zsrc)**2.0/etaV);
	!if (abs(Rbig) <= 10**(-6.0)) then
	!	Rbig = 10**(-6.0);
	!endif
	G = exp(-1.0*sqrt(gamA)*Rbig)/(4.0*pi*Rbig*sqrt(etaH/etaV));
	Gin(ix) = -1.0*x(ix)*etaH/etaV*(zsrc-zrcv)*(gamA/(Rbig**2.0)+3.0*sqrt(gamA)/(Rbig**3.0)+3.0/(Rbig**4))*G/etaV;
end do

END SUBROUTINE gin13
