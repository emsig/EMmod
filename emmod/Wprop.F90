SUBROUTINE Wprop (Wup, Wdown, Gam, nkx, cor, nz, z, zrcv, zrcvlay)
! Calculate the propagation of the EM-field from the layer boundaries to the receiver position as given in equation 74.
!
! Calling arguments:
! Wup:     Propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nkx
! Wdown:   Propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nkx
! Gam:     Vertical wavenumber Gamma, complex array of size nkx
! nkx:     Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nkx
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcv:    Depth of receiver, real scalar
! zrcvlay: The number of the layer where the receivers are located, integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nkx, nz, zrcvlay
COMPLEX,INTENT(INOUT) :: Wup(nkx), Wdown(nkx)
COMPLEX,INTENT(IN)    :: Gam(nkx,nz)
REAL   ,INTENT(IN)    :: cor(nkx), z(nz), zrcv
REAL                  :: dp, dm
INTEGER               :: ikx, zrcvlayfort

zrcvlayfort = zrcvlay + 1;
dp = z(zrcvlayfort+1)-zrcv;
dm = zrcv-z(zrcvlayfort);

do ikx=1,nkx
	Wup(ikx) = exp(-1.0*Gam(ikx,zrcvlayfort)*dp);
	Wdown(ikx) = exp(-1.0*Gam(ikx,zrcvlayfort)*dm);
end do

END SUBROUTINE Wprop
