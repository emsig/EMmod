SUBROUTINE hankeltransxxm (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 41 (inline oriented magnetic receiver, inline oriented electric source) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalxxm, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalxxm, complex array of size corel
! corel:   Amount of points in the wavenumber domain (61 times the amount of integration subdomains nd), integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size corel
! xpos:    Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:   Amount of samples in logarithmic coordinate vector, integer scalar
! nd:      Amount of integration subdomains, integer scalar
! kmax:    Largest wavenumber, real scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: marker(nlogx), corel, nlogx, nd
COMPLEX,INTENT(INOUT) :: xtotrad(nlogx)
COMPLEX,INTENT(IN)    :: temptot(corel)
REAL   ,INTENT(IN)    :: cor(corel), xpos(corel), kmax
COMPLEX               :: integrand(corel), tempin(61), temp
REAL                  :: kintel, pi, BesJ2
INTEGER               :: ikx, n, ind, il

pi = 2*acos(0.0);
kintel = kmax/nd;
do n=1,nlogx
	if (marker(n) == 2) then
		integrand = 0;
		do ikx=1,corel
			integrand(ikx) = temptot(ikx)*BesJ2(cor(ikx)*xpos(n))*cor(ikx);
		end do
		xtotrad(n) = 0.0;
		do ind =1,nd
			temp = 0.0;
			do il=1,61
				tempin(il) = integrand((ind-1)*61+il);
			end do
			call zqk61n(tempin,(ind-1)*kintel,ind*kintel,temp)
			xtotrad(n) = xtotrad(n) + temp;
		end do
		xtotrad(n) = xtotrad(n)/(4*pi);
	endif
end do
END SUBROUTINE hankeltransxxm

SUBROUTINE hankeltransxym (xtotrad1,xtotrad2,marker,temptot1,temptot2,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 42 (inline oriented magnetic receiver, crossline oriented electric source) 
!
! Calling arguments:
! xtotrad1: Hankel Transform of the 1st output term of Ptotalxym, complex array of size nlogx
! xtotrad2: Hankel Transform of the 2nd output term of Ptotalxym, complex array of size nlogx
! marker:   Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot1: 1st output term of Ptotalxym, complex array of size corel
! temptot2: 2nd output term of Ptotalxym, complex array of size corel
! corel:    Amount of points in the wavenumber domain (61 times the amount of integration subdomains nd), integer scalar
! cor:      Coordinates in the wavenumber domain, real array of size corel
! xpos:     Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:    Amount of samples in logarithmic coordinate vector, integer scalar
! nd:       Amount of integration subdomains, integer scalar
! kmax:     Largest wavenumber, real scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: marker(nlogx), corel, nlogx, nd
COMPLEX,INTENT(INOUT) :: xtotrad1(nlogx), xtotrad2(nlogx)
COMPLEX,INTENT(IN)    :: temptot1(corel), temptot2(corel)
REAL   ,INTENT(IN)    :: cor(corel), xpos(nlogx), kmax
COMPLEX               :: integrand1(corel), integrand2(corel), tempin1(61), tempin2(61), temp
REAL                  :: kintel, pi, BesJ0, BesJ2
INTEGER               :: ikx, n, ind, il

pi = 2*acos(0.0);
kintel = kmax/nd;
do n=1,nlogx
	if (marker(n) == 2) then
		integrand1 = 0;
		integrand2 = 0;
		do ikx=1,corel
			integrand1(ikx) = temptot1(ikx)*BesJ0(cor(ikx)*xpos(n))*cor(ikx);
			integrand2(ikx) = temptot2(ikx)*BesJ2(cor(ikx)*xpos(n))*cor(ikx);
		end do
		xtotrad1(n) = 0.0;
		xtotrad2(n) = 0.0;
		do ind =1,nd
			temp = 0.0;
			do il=1,61
				tempin1(il) = integrand1((ind-1)*61+il);
				tempin2(il) = integrand2((ind-1)*61+il);
			end do
			call zqk61n(tempin1,(ind-1)*kintel,ind*kintel,temp)
			xtotrad1(n) = xtotrad1(n) + temp;
			call zqk61n(tempin2,(ind-1)*kintel,ind*kintel,temp)
			xtotrad2(n) = xtotrad2(n) + temp;
		end do
		xtotrad1(n) = xtotrad1(n)/(4*pi);
		xtotrad2(n) = xtotrad2(n)/(4*pi);
	endif
end do
END SUBROUTINE hankeltransxym

SUBROUTINE hankeltransxzm (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 43 (inline oriented magnetic receiver, vertical electric source) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalxzm, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalxzm, complex array of size corel
! corel:   Amount of points in the wavenumber domain (61 times the amount of integration subdomains nd), integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size corel
! xpos:    Logarithmic coordinate vector in the space domain, real array of size nlogx
! nlogx:   Amount of samples in logarithmic coordinate vector, integer scalar
! nd:      Amount of integration subdomains, integer scalar
! kmax:    Largest wavenumber, real scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: marker(nlogx), corel, nlogx, nd
COMPLEX,INTENT(INOUT) :: xtotrad(nlogx)
COMPLEX,INTENT(IN)    :: temptot(corel)
REAL   ,INTENT(IN)    :: cor(corel), xpos(nlogx), kmax
COMPLEX               :: integrand(corel), tempin(61), temp
REAL                  :: kintel, pi, BesJ1
INTEGER               :: ikx, n, ind, il

pi = 2*acos(0.0);
kintel = kmax/nd;
do n=1,nlogx
	if (marker(n) == 2) then
		integrand = 0;
		do ikx=1,corel
			integrand(ikx) = temptot(ikx)*BesJ1(cor(ikx)*xpos(n))*cor(ikx)*cor(ikx);
		end do
		xtotrad(n) = 0.0;
		do ind =1,nd
			temp = 0.0;
			do il=1,61
				tempin(il) = integrand((ind-1)*61+il);
			end do
			call zqk61n(tempin,(ind-1)*kintel,ind*kintel,temp)
			xtotrad(n) = xtotrad(n) + temp;
		end do
		xtotrad(n) = xtotrad(n)/(2*pi);
	endif
end do
END SUBROUTINE hankeltransxzm
