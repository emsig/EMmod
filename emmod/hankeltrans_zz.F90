SUBROUTINE hankeltranszx (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 31 (vertical electric receiver, inline oriented electric source) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalzx, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalzx, complex array of size corel
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
END SUBROUTINE hankeltranszx

SUBROUTINE hankeltranszy (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 32 (vertical electric receiver, crossline oriented electric source) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalzy, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalzy, complex array of size corel
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
END SUBROUTINE hankeltranszy

SUBROUTINE hankeltranszz (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 33 (vertical electric receiver, vertical electric source) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalzz, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalzz, complex array of size corel
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
REAL                  :: kintel, pi, BesJ0
INTEGER               :: ikx, n, ind, il

pi = 2*acos(0.0);
kintel = kmax/nd;
do n=1,nlogx
	if (marker(n) == 2) then
		integrand = 0;
		do ikx=1,corel
			integrand(ikx) = temptot(ikx)*BesJ0(cor(ikx)*xpos(n))*cor(ikx)*cor(ikx)*cor(ikx);
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
END SUBROUTINE hankeltranszz
