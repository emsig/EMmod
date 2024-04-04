SUBROUTINE Pupplus (Pup, Rp, Rm, Gam, nk, cor, nz, z, zsrc, zrcvlay, zsrclay, nlayers, above);
! Calculate Pup+
! for receivers in the source layer: equation 81 
! for receivers above the source layer: equations 95 and 96
!
! Calling arguments:
! Pup:     Upgoing field, complex array of size nk
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
COMPLEX,INTENT(INOUT) :: Pup(nk), Rp(nk,nlayers), Rm(nk,nlayers)
COMPLEX,INTENT(IN)    :: Gam(nk,nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz)
REAL                  :: zsrc, ds, dp, dm
INTEGER               :: iz, ikx, zsrclayfort, zrcvlayfort
COMPLEX               :: Ms

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
ds = z(zsrclayfort+1)-z(zsrclayfort);
dp = z(zsrclayfort+1)-zsrc;
dm = zsrc-z(zsrclayfort);
if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			Ms = 1.0;
			Pup(ikx) = Rp(ikx,1)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dp));
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			Pup(ikx) = 0.0;
		else ! If the source and the receivers are in any layer in between
			Ms = 1.0-Rm(ikx,1)*Rp(ikx,1)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
			Pup(ikx) = Rp(ikx,1)/Ms*(exp(-1.0*Gam(ikx,zsrclayfort)*dp)&
				+Rm(ikx,1)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dm)));
		endif
		! Note that if the source and the receivers are in the top most layer, Rm is zero. 
		! The if statement is introduced to avoid the exponential function to produce values
		! like infinity. Even if Rm is zero, a product with infinity causes problems. 
	end do
elseif (above == 1) then ! receivers above the source layer
	! First compute P_{s-1}
	do ikx=1,nk
		if (zsrclayfort == nz) then ! If the source is in the lower most layer
			Ms = 1.0;
			Pup(ikx) = (1+Rm(ikx,nlayers))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm));
		else
			Ms = 1.0-Rm(ikx,nlayers)*Rp(ikx,nlayers)*exp(-2.0*Gam(ikx,zsrclayfort)*ds);
			Pup(ikx) = (1+Rm(ikx,nlayers))*(exp(-1.0*Gam(ikx,zsrclayfort)*dm)&
				+Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zsrclayfort)*(ds+dp)));
		endif
		Pup(ikx) = Pup(ikx)/(Ms*(1+Rm(ikx,nlayers-1)&
			*exp(-2.0*Gam(ikx,zsrclayfort-1)*(z(zsrclayfort)-z(zsrclayfort-1)))));
	end do
	if (nlayers > 2) then
		! Second compute P for all other layers
		do iz = 1,nlayers-2
			do ikx=1,nk
				Pup(ikx) = Pup(ikx)*(1+Rm(ikx,iz+1))&
					*exp(-1.0*Gam(ikx,zrcvlayfort+iz)*(z(zrcvlayfort+iz+1)-z(zrcvlayfort+iz)));
				if (zrcvlayfort+iz-1 /= 1) then ! If the receivers are NOT in the first layer
					Pup(ikx) = Pup(ikx)/(1+Rm(ikx,iz)&
						*exp(-2.0*Gam(ikx,zrcvlayfort+iz-1)*(z(zrcvlayfort+iz)-z(zrcvlayfort+iz-1)))); 
					! Note that if the receivers are in the first layer, Rm is zero and this line is just 
					! a division by 1. The if statement is introduced to avoid the exponential function
					! to cause a problem, because z(zrcvlayfort+iz-1) has only a dummy value.
				endif
			end do
		end do
	endif
elseif (above == -1) then ! receivers below the source layer: Pupplus is not used
	do ikx=1,nk
		Pup(ikx) = 0.0;
	end do
endif

END SUBROUTINE Pupplus
