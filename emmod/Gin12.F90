SUBROUTINE gin12 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, x, y, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), y(nlogx), zrcv, zsrc
COMPLEX               :: RbigA, RbigB, GA, GB, term1, term2, term3, term4
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

	GA = exp(-1.0*sqrt(gamA)*RbigA)/(4.0*pi*RbigA*sqrt(etaH/etaV));
	GB = exp(-1.0*sqrt(gamB)*RbigB)/(4.0*pi*RbigB*sqrt(zetaH/zetaV));

	term1 = ((3.0*x(ix)*y(ix)/(RbigA**2))*(1.0/(etaV*RbigA**2)&
		&+sqrt(gamA)/(etaV*RbigA))+zetaH*x(ix)*y(ix)/(RbigA**2))*GA;
	term2 = -1.0*zetaH*x(ix)*y(ix)/(rad*rad)*(GA-GB);
	term3 = -1.0*sqrt(zetaH)*(2*x(ix)*y(ix)/(rad*rad))*((exp(-1.0*sqrt(gamA)*RbigA)&
		&-exp(-1.0*sqrt(gamB)*RbigB))/(4.0*pi*sqrt(etaH)*rad*rad));
	Gin(ix) = term1 + term2 + term3;
end do

END SUBROUTINE gin12
