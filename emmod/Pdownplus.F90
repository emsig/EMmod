SUBROUTINE Pdownplus (Pdown, Rp, Rm, Gam, nk, cor, nz, z, zsrc, zrcvlay, zsrclay, nlayers, above);
! Calculate Pdown+
! for receivers in the source layer: equation 82 
! for receivers below the source layer: equations 103 and 104
!
! Calling arguments:
! Pdown:   Downgoing field, complex array of size nk
! Rp:      Upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:      Downgoing global reflection coefficient, complex array of size nk times nlayers
! Gam:     Gamma, the vertical wavenumber, complex array of size nk times nz
! nk:      Amount of wavenumbers, integer scalar
! cor:     Coordinates in the wavenumber domain, real array of size nk
! nz:      Amount of layers, integer scalar
! z:       Depth of interfaces, real array of size nz
! zsrc:    Depth of source, real scalar
! zrcvlay: The number of the layer where the receivers are located, integer scalar
! zsrclay: The number of the layer where the source is located, integer scalar
! nlayers: Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:   Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!          in the same layer (above=0), integer scalar

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zrcvlay, zsrclay, nlayers, above
COMPLEX,INTENT(INOUT) :: Pdown(nk), Rp(nk,nlayers), Rm(nk,nlayers)
COMPLEX,INTENT(IN)    :: Gam(nk,nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: zsrc, ds, dp, dm
INTEGER               :: iz, ikx, zsrclayfort, zrcvlayfort
COMPLEX               :: Ms

zsrclayfort = zsrclay + 1;
zrcvlayfort = zrcvlay + 1;
ds = z(zsrclayfort+1)-z(zsrclayfort);
dp = z(zsrclayfort+1)-zsrc;
dm = zsrc-z(zsrclayfort);
if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			Pdown(ikx) = 0.0;
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			Ms = 1.0;
			Pdown(ikx) = Rm(ikx,1)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
		else ! If the source and the receivers are in any layer in between
			Ms = 1.0-Rm(ikx,1)*Rp(ikx,1)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
			Pdown(ikx) = Rm(ikx,1)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
				+Rp(ikx,1)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
		endif
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		Pdown(ikx) = 0.0;
	end do
elseif (above == -1) then ! receivers below the source layer
	! First compute P_{s+1}
	do ikx=1,nk
		if (zsrclayfort == 1) then ! If the source is in the top most layer
			Ms = 1.0;
			Pdown(ikx) = (1+Rp(ikx,1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
		else
			Ms = 1.0-Rm(ikx,1)*Rp(ikx,1)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
			Pdown(ikx) = (1+Rp(ikx,1))*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
				+Rm(ikx,1)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
		endif
		if (zsrclayfort+1 < nz) then ! if the source is in the last but one layer, Rp(:,:,2) is zero.
			Pdown(ikx) = Pdown(ikx)/(Ms*(1+Rp(ikx,2)&
				*exp(-2.0*Gam(ikx,zsrclayfort+1)*(z(zsrclayfort+2)-z(zsrclayfort+1)))));
		else
			Pdown(ikx) = Pdown(ikx)/Ms;
		end if
	end do
	if (nlayers > 2) then
		! Second compute P for all other layers
		do iz = 3,nlayers
			do ikx=1,nk
				Pdown(ikx) = Pdown(ikx)*(1+Rp(ikx,iz-1))&
					*exp(-1.0*Gam(ikx,zsrclayfort+iz-2)*(z(zsrclayfort+iz-1)-z(zsrclayfort+iz-2)));
				!if (zrcvlayfort+iz <= nz) then ! If the receiver is NOT in the last layer
				!if (zsrclayfort+iz <= nz) then ! If the receiver is NOT in the last layer
				if (zsrclayfort+iz-1 /= nz) then ! If the receiver is NOT in the last layer
					Pdown(ikx) = Pdown(ikx)/(1+Rp(ikx,iz)&
						*exp(-2.0*Gam(ikx,zsrclayfort+iz-1)*(z(zsrclayfort+iz)-z(zsrclayfort+iz-1))));
				endif
			end do
		end do
	endif
endif

END SUBROUTINE Pdownplus
