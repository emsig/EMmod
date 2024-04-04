SUBROUTINE gin51 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, x, y, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: RbigA, RbigB, fac, term1, term2, term3, term4
REAL                  :: pi, rad, d, signd
INTEGER               :: ix

pi = 2*acos(0.0);

d = zrcv-zsrc;
if (d < 0) then
        signd = -1.0;
else
        signd = 1.0;
endif


do ix=1,nlogx
	if (abs(xpos(ix))<= 10**(-6.0)) then
		rad = 10**(-6.0);
	else
		rad = xpos(ix);
	endif
	RbigA = sqrt(rad**2.0+(etaH*(zrcv-zsrc)**2.0/etaV));
	RbigB = sqrt(rad**2.0+(zetaH*(zrcv-zsrc)**2.0/zetaV));

	fac = -1.0*(zrcv-zsrc)*etaH/(4.0*pi*sqrt(etaH*zetaH)*etaV)*exp(-1.0*sqrt(gamA)*RbigA);
	term1 = fac*(sqrt(gamA)*x(ix)**2.0/(rad**2.0*RbigA**3.0)-(rad**(-2.0)-2.0*x(ix)**2.0*rad**(-4.0)&
		&-sqrt(gamA)*x(ix)**2.0/(rad**2.0*RbigA))*sqrt(gamA)/RbigA);
	term2 = signd/(4.0*pi)*exp(sqrt(etaH*zetaH)*abs(zrcv-zsrc))*(rad**(-2.0)-2.0*x(ix)**2.0*rad**(-4.0));
	fac = -1.0*(zrcv-zsrc)*zetaH/(4.0*pi*sqrt(etaH*zetaH)*zetaV)*exp(-1.0*sqrt(gamB)*RbigB);
	term3 = fac*(sqrt(gamB)*y(ix)**2.0/(rad**2.0*RbigB**3.0)-(rad**(-2.0)-2.0*y(ix)**2.0*rad**(-4.0)&
		&-sqrt(gamB)*y(ix)**2.0/(rad**2.0*RbigB))*sqrt(gamB)/RbigB);
	term4 = signd/(4.0*pi)*exp(sqrt(etaH*zetaH)*abs(zrcv-zsrc))*(rad**(-2.0)-2.0*y(ix)**2.0*rad**(-4.0));
	Gin(ix) = term1 + term2 + term3 + term4;
end do

END SUBROUTINE gin51
