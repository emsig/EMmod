SUBROUTINE Ptotalref (Ptot, Wup, Rp, nk);
! calculate the reflection response at the receiver level
!
! Calling arguments:
! Ptot: Total reflection response, complex array of size nk
! Wup:  Propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Rp:   Upgoing global reflection coefficient, complex array of size nk
! nk:   Amount of wavenumbers, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: Wup(nk)
COMPLEX,INTENT(IN)    :: Rp(nk)
INTEGER               :: ikx

do ikx=1,nk
	Ptot(ikx) = Rp(ikx)*Wup(ikx);
end do

END SUBROUTINE Ptotalref
