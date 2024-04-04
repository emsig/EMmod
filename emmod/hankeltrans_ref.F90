SUBROUTINE hankeltransref (xtotrad,marker,temptot,corel,cor,xpos,nlogx,nd,kmax)
! Perform Hankel Transformation for component 77 and 88 (TM-mode and TE-mode reflection response, respectively) 
!
! Calling arguments:
! xtotrad: Hankel Transform of the output term of Ptotalref, complex array of size nlogx
! marker:  Indicates for which points in space the Hankel Transformation needs to be computed, integer array of size nlogx
! temptot: Output term of Ptotalref, complex array of size corel
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
REAL                  :: kintel, pi, Besj0
INTEGER               :: ikx, n, ind, il

pi = 2*acos(0.0);
kintel = kmax/nd;
do n=1,nlogx
	if (marker(n) == 2) then
		integrand = 0;
		do ikx=1,corel
			integrand(ikx) = temptot(ikx)*BesJ0(cor(ikx)*xpos(n))*cor(ikx);
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
END SUBROUTINE hankeltransref
