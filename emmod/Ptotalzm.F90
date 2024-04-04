SUBROUTINE Ptotalzzm (Ptot, nxh, nyh);
! Calculate total field for the 63 component (vertical magnetic receiver, vertical electric source)
!
! Calling arguments:
! Ptot: Total field, complex array of size nk
! nxh:  Quarter of the size of the receiver grid in inline direction in the space domain
! nxh:  Quarter of the size of the receiver grid in crossline direction in the space domain 
!
! Note: This component is zero everywhere.

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nxh, nyh
COMPLEX,INTENT(INOUT) :: Ptot(nxh,nyh)
INTEGER               :: ikx, iky

do iky=1,nyh
	do ikx=1,nxh
		Ptot(ikx,iky) = 0.0;
	end do
end do

END SUBROUTINE Ptotalzzm

SUBROUTINE Ptotalzxm (Ptot, Pdownbar, Pupbar, Wupbar, Wdownbar, GammaB, Rpbar, Rmbar, zetaH, zetaV, nk, cor, &
	nz, z, zrcv, zsrc, zsrclay, zrcvlay, gamB, nlayers, above, xdirect);
! Calculate total field for the 61 component (vertical magnetic receiver, inline oriented electric source)
!
! Calling arguments:
! Ptot:     Total field, complex array of size nk
! Pdownbar: TE-mode part of the downgoing field, complex array of size nk
! Pupbar:   TE-mode part of the upgoing field, complex array of size nk
! Wupbar:   TE-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdownbar: TE-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! GammaB:   Vertical wavenumber Gamma Bar, complex array of size nk
! Rpbar:    TE-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rmbar:    TE-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! zetaH:    Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:    Material parameter zeta for the vertical direction, complex array of size nz
! nk:       Amount of wavenumbers, integer scalar
! cor:      Coordinates in the wavenumber domain, real array of size nk
! nz:       Amount of layers, integer scalar
! z:        Depth of interfaces, real array of size nz
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! zsrclay:  The number of the layer where the source is located, integer scalar
! zrcvlay:  The number of the layer where the receivers are located, integer scala
! gamB:     Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:  Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: GammaB(nk,nz), Pupbar(nk), Pdownbar(nk), Wupbar(nk), Wdownbar(nk)
COMPLEX, INTENT(IN)   :: Rpbar(nk,nlayers), Rmbar(nk,nlayers), zetaH(nz), zetaV(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) then
	signd = -1.0;
else 
	signd = 1.0;
endif
if (d < 0) d = -1.0*d; ! replaces abs(d)

if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zsrclayfort)*GammaB(ikx,zsrclayfort));
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx));
			else
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pdownbar(ikx)*Wdownbar(ikx));
			else
				Ptot(ikx) = factor*(Pdownbar(ikx)*Wdownbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+Pdownbar(ikx)*Wdownbar(ikx));
			else
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+Pdownbar(ikx)*Wdownbar(ikx)&
					&+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		endif
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zrcvlayfort)*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			Ptot(ikx) = factor*Pupbar(ikx)*Wupbar(ikx);
		else
			Ptot(ikx) = factor*Pupbar(ikx)*(Wupbar(ikx)+Rmbar(ikx,1)&
				*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdownbar(ikx));
		endif
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zrcvlayfort)*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			Ptot(ikx) = factor*Pdownbar(ikx)*Wdownbar(ikx);
		else
			Ptot(ikx) = factor*Pdownbar(ikx)*(Wdownbar(ikx)+Rpbar(ikx,nlayers)&
				*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wupbar(ikx));
		endif
	end do
endif

END SUBROUTINE Ptotalzxm

SUBROUTINE Ptotalzym (Ptot, Pdownbar, Pupbar, Wupbar, Wdownbar, GammaB, Rpbar, Rmbar, zetaH, zetaV, nk, cor, &
	nz, z, zrcv, zsrc, zsrclay, zrcvlay, gamB, nlayers, above, xdirect);
! Calculate total field for the 62 component (vertical magnetic receiver, crossline oriented electric source)
!
! Calling arguments:
! Ptot:     Total field, complex array of size nk
! Pdownbar: TE-mode part of the downgoing field, complex array of size nk
! Pupbar:   TE-mode part of the upgoing field, complex array of size nk
! Wupbar:   TE-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdownbar: TE-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! GammaB:   Vertical wavenumber Gamma Bar, complex array of size nk
! Rpbar:    TE-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rmbar:    TE-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! zetaH:    Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:    Material parameter zeta for the vertical direction, complex array of size nz
! nk:       Amount of wavenumbers, integer scalar
! cor:      Coordinates in the wavenumber domain, real array of size nk
! nz:       Amount of layers, integer scalar
! z:        Depth of interfaces, real array of size nz
! zrcv:     Depth of receivers, real scalar
! zsrc:     Depth of source, real scalar
! zsrclay:  The number of the layer where the source is located, integer scalar
! zrcvlay:  The number of the layer where the receivers are located, integer scala
! gamB:     Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:  Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:    Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!           in the same layer (above=0), integer scalar
! xdirect:  Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: GammaB(nk,nz), Pupbar(nk), Pdownbar(nk), Wupbar(nk), Wdownbar(nk)
COMPLEX, INTENT(IN)   :: Rpbar(nk,nlayers), Rmbar(nk,nlayers), zetaH(nz), zetaV(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) then
	signd = -1.0;
else 
	signd = 1.0;
endif
if (d < 0) d = -1.0*d; ! replaces abs(d)

if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zsrclayfort)*GammaB(ikx,zsrclayfort));
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx));
			else
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pdownbar(ikx)*Wdownbar(ikx));
			else
				Ptot(ikx) = factor*(Pdownbar(ikx)*Wdownbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+Pdownbar(ikx)*Wdownbar(ikx));
			else
				Ptot(ikx) = factor*(Pupbar(ikx)*Wupbar(ikx)+Pdownbar(ikx)*Wdownbar(ikx)&
				&+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		endif
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zrcvlayfort)*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			Ptot(ikx) = factor*Pupbar(ikx)*Wupbar(ikx);
		else
			Ptot(ikx) = factor*Pupbar(ikx)*(Wupbar(ikx)+Rmbar(ikx,1)&
				*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdownbar(ikx));
		endif
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor = zetaH(zsrclayfort)/(2.0*zetaV(zrcvlayfort)*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			Ptot(ikx) = factor*Pdownbar(ikx)*Wdownbar(ikx);
		else
			Ptot(ikx) = factor*Pdownbar(ikx)*(Wdownbar(ikx)+Rpbar(ikx,nlayers)&
				*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wupbar(ikx));
		endif
	end do
endif

END SUBROUTINE Ptotalzym
