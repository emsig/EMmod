      
SUBROUTINE Rplus (Rp, nk, cor, etaH, Gam, nz, z, zrcvlay, zsrclay, nlayers, above)
! Calculate global reflection coefficient of layers below the receivers using equations 64, 65, A-11 and A-12
!
! Calling arguments:
! Rp:      Upgoing global reflection coefficient, complex array of size nk times nlayers
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! etaH:    Material parameter eta for the horizontal direction, complex array of size nz
! Gam:     Gamma, the vertical wavenumber, complex array of size nk times nz
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zrcvlay: The number of the layer where the receivers are located, integer scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar
!
! Note: To compute Rpbar replace etaH with zetaH and Gamma with Gamma bar.

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay, zsrclay, nlayers, above
COMPLEX,INTENT(INOUT) :: Rp(nk,nlayers)
COMPLEX,INTENT(IN)    :: etaH(nz), Gam(nk,nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: kx
INTEGER               :: iz, izout, ikx, zrcvlayfort, zsrclayfort
COMPLEX               :: term1, term2 
COMPLEX, DIMENSION(:), ALLOCATABLE :: rloc, Rptemp

ALLOCATE(rloc (nk))

zrcvlayfort = zrcvlay + 1;
zsrclayfort = zsrclay + 1;
if (above == 0) then ! source and receivers in the same layer
	do iz=nz-1,zsrclayfort,-1
		do ikx=1,nk
			kx = cor(ikx);
			rloc(ikx) = (etaH(iz+1)*Gam(ikx,iz)-etaH(iz)*Gam(ikx,iz+1))/(etaH(iz+1)*Gam(ikx,iz)+etaH(iz)*Gam(ikx,iz+1));
			if (iz == nz-1) then
				term1 = cmplx(0.0,0.0);
				term2 = cmplx(0.0,0.0);
			else
				term1 = Rp(ikx,1)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
				term2 = rloc(ikx)*Rp(ikx,1)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
			endif
			Rp(ikx,1) = (rloc(ikx)+term1)/(1.0+term2);
		end do
	end do
elseif (above == -1) then ! receivers below source layer
	ALLOCATE(Rptemp (nk))
	izout = nlayers + 1;
	if (zrcvlayfort == nz) then
		izout = izout - 1;
		do ikx=1,nk
			Rp(ikx,izout) = cmplx(0.0,0.0);
		end do
	endif
	do iz=nz-1,zsrclayfort,-1
		if (iz <= zrcvlayfort) then
			izout = izout - 1;
		endif
		do ikx=1,nk
			kx = cor(ikx);
			rloc(ikx) = (etaH(iz+1)*Gam(ikx,iz)-etaH(iz)*Gam(ikx,iz+1))/(etaH(iz+1)*Gam(ikx,iz)+etaH(iz)*Gam(ikx,iz+1));
			if (iz == nz-1) then
				term1 = cmplx(0.0,0.0);
				term2 = cmplx(0.0,0.0);
			else
				term1 = Rptemp(ikx)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
				term2 = rloc(ikx)*Rptemp(ikx)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
			endif
			Rptemp(ikx) = (rloc(ikx)+term1)/(1.0+term2);
			! The global reflection coefficient is given back for all layers between and including the source- and the 
			! receiver-layer
			if (iz <= zrcvlayfort) then
				Rp(ikx,izout) = Rptemp(ikx);
			endif
		end do
	end do
elseif (above == 1) then ! receivers above source layer
	ALLOCATE(Rptemp (nk))
	izout = nlayers + 1;
	do iz=nz-1,zrcvlayfort,-1
		if (iz <= zsrclayfort) then
			izout = izout - 1;
		endif
		do ikx=1,nk
			kx = cor(ikx);
			rloc(ikx) = (etaH(iz+1)*Gam(ikx,iz)-etaH(iz)*Gam(ikx,iz+1))/(etaH(iz+1)*Gam(ikx,iz)+etaH(iz)*Gam(ikx,iz+1));
			if (iz == nz-1) then
				term1 = cmplx(0.0,0.0);
				term2 = cmplx(0.0,0.0);
			else
				term1 = Rptemp(ikx)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
				term2 = rloc(ikx)*Rptemp(ikx)*exp(-2.0*Gam(ikx,iz+1)*(z(iz+2)-z(iz+1)));
			endif
			Rptemp(ikx) = (rloc(ikx)+term1)/(1.0+term2);
			! The global reflection coefficient is given back for all layers between and including the source- and the 
			! receiver-layer
			if (iz <= zsrclayfort) then
				Rp(ikx,izout) = Rptemp(ikx);
			endif
		end do
	end do
endif

END SUBROUTINE Rplus
