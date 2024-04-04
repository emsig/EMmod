SUBROUTINE gin11 (Gin, zrcv, zsrc, etaH, etaV, zetaH, zetaV, gamA, gamB, xpos, x, nlogx)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nlogx
COMPLEX,INTENT(INOUT) :: Gin(nlogx)
COMPLEX, INTENT(IN)   :: etaH, etaV, zetaH, zetaV, gamA, gamB
REAL   ,INTENT(IN)    :: xpos(nlogx), x(nlogx), zrcv, zsrc
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
	!if (abs(RbigA) <= 10**(-6.0)) then
	!	RbigA = 10**(-6.0);
	!endif
	RbigB = sqrt(rad**2.0+(zetaH*(zrcv-zsrc)**2.0/zetaV));
	!if (abs(RbigB) <= 10**(-6.0)) then
	!	RbigB = 10**(-6.0);
	!endif

	GA = exp(-1.0*sqrt(gamA)*RbigA)/(4.0*pi*RbigA*sqrt(etaH/etaV));
	GB = exp(-1.0*sqrt(gamB)*RbigB)/(4.0*pi*RbigB*sqrt(zetaH/zetaV));
	term1 = ((3.0*x(ix)*x(ix)/(RbigA**2.0)-1.0)*(1.0/(etaV*RbigA**2.0)&
		&+sqrt(gamA)/(etaV*RbigA))+zetaH*x(ix)*x(ix)/(RbigA**2.0))*GA;
	term2 = -1.0*zetaH*GB;
	term3 = -1.0*zetaH*x(ix)*x(ix)/(rad*rad)*(GA-GB);
	term4 = -1.0*sqrt(zetaH)*(2.0*x(ix)*x(ix)/(rad*rad)-1.0)*((exp(-1.0*sqrt(gamA)*RbigA)&
		&-exp(-1.0*sqrt(gamB)*RbigB))/(4.0*pi*sqrt(etaH)*rad*rad));
	Gin(ix) = term1 + term2 + term3 + term4;
end do

END SUBROUTINE gin11
