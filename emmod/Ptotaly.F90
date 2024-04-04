SUBROUTINE Ptotalyy (Ptot1, Ptot2, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, Wup, Wdown, Wupbar, Wdownbar, &
	Gam, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, nk, cor, nz, z, zrcv, zsrc, &
	zsrclay, zrcvlay, gamA, gamB, nlayers, above, xdirect);
! Calculate total field for the 22 component (crossline oriented electric receiver, crossline oriented electric source)
!
! Calling arguments:  
! Ptot1:        1st term of the total field, complex array of size nk
! Ptot2:        2nd term of the total field, complex array of size nk
! Pdownmin:     TM-mode part of the downgoing field, complex array of size nk
! Pupmin:       TM-mode part of the upgoing field, complex array of size nk
! Pdownplusbar: TE-mode part of the downgoing field, complex array of size nk
! Pupplusbar:   TE-mode part of the upgoing field, complex array of size nk
! Wup:          TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:        TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Wupbar:       TE-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdownbar:     TE-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:          Vertical wavenumber Gamma, complex array of size nk
! GammaB:       Vertical wavenumber Gamma Bar, complex array of size nk
! Rp:           TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rpbar:        TE-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:           TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! Rmbar:        TE-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:         Material parameter eta for the horizontal direction, complex array of size nz
! etaV:         Material parameter eta for the vertical direction, complex array of size nz
! zetaH:        Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:        Material parameter zeta for the vertical direction, complex array of size nz
! nk:           Amount of wavenumbers, integer scalar
! cor:          Coordinates in the wavenumber domain, real array of size nk
! nz:           Amount of layers, integer scalar
! z:            Depth of interfaces, real array of size nz
! zrcv:         Depth of receivers, real scalar
! zsrc:         Depth of source, real scalar
! zsrclay:      The number of the layer where the source is located, integer scalar
! zrcvlay:      The number of the layer where the receivers are located, integer scala
! gamA:         Small gamma squared (zetaH*etaV), complex array of size nz
! gamB:         Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:      Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:        Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!               in the same layer (above=0), integer scalar
! xdirect:      Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot1(nk), Ptot2(nk)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), GammaB(nk,nz)
COMPLEX,INTENT(IN)    :: Pdownmin(nk), Pupmin(nk), Pdownplusbar(nk), Pupplusbar(nk)
COMPLEX,INTENT(IN)    :: Wup(nk), Wdown(nk), Wupbar(nk), Wdownbar(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nlayers), Rpbar(nk,nlayers), Rm(nk,nlayers), Rmbar(nk,nlayers) 
COMPLEX,INTENT(IN)    :: etaH(nz), etaV(nz), zetaH(nz), zetaV(nz), gamA(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
COMPLEX               :: temp1, temp2
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort
COMPLEX               :: factor1, factor2

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) d = -1.0*d; ! replaces abs(d)

if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		factor1 = Gam(ikx,zsrclayfort)/(2.0*etaH(zsrclayfort));
		factor2 = zetaH(zsrclayfort)/(2.0*GammaB(ikx,zsrclayfort));
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = factor1*(Pupmin(ikx)*Wup(ikx));
				temp2 = factor2*(Pupplusbar(ikx)*Wupbar(ikx));
			else
				temp1 = factor1*(Pupmin(ikx)*Wup(ikx)-exp(-1.0*Gam(ikx,zsrclayfort)*d));
				temp2 = factor2*(Pupplusbar(ikx)*Wupbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = factor1*(Pdownmin(ikx)*Wdown(ikx));
				temp2 = factor2*(Pdownplusbar(ikx)*Wdownbar(ikx));
			else
				temp1 = factor1*(Pdownmin(ikx)*Wdown(ikx)-exp(-1.0*Gam(ikx,zsrclayfort)*d));
				temp2 = factor2*(Pdownplusbar(ikx)*Wdownbar(ikx)+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = factor1*(Pupmin(ikx)*Wup(ikx)+Pdownmin(ikx)*Wdown(ikx));
				temp2 = factor2*(Pupplusbar(ikx)*Wupbar(ikx)+Pdownplusbar(ikx)*Wdownbar(ikx));
			else
				temp1 = factor1*(Pupmin(ikx)*Wup(ikx)+Pdownmin(ikx)*Wdown(ikx)&
					&-exp(-1.0*Gam(ikx,zsrclayfort)*d));
				temp2 = factor2*(Pupplusbar(ikx)*Wupbar(ikx)+Pdownplusbar(ikx)*Wdownbar(ikx)&
					&+exp(-1.0*GammaB(ikx,zsrclayfort)*d));
			endif
		endif
		Ptot1(ikx) = temp1 - temp2;
		Ptot2(ikx) = temp1 + temp2;
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor1 = Gam(ikx,zrcvlayfort)/(2.0*etaH(zrcvlayfort));
		factor2 = zetaH(zsrclayfort)/(2.0*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			temp1 = factor1*Pupmin(ikx)*Wup(ikx);
			temp2 = factor2*Pupplusbar(ikx)*Wupbar(ikx);
		else
			temp1 = factor1*Pupmin(ikx)*(Wup(ikx)&
				-Rm(ikx,1)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
			temp2 = factor2*Pupplusbar(ikx)*(Wupbar(ikx)&
				+Rmbar(ikx,1)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdownbar(ikx));
		endif
		Ptot1(ikx) = temp1 - temp2;
		Ptot2(ikx) = temp1 + temp2;
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor1 = Gam(ikx,zrcvlayfort)/(2.0*etaH(zrcvlayfort));
		factor2 = zetaH(zsrclayfort)/(2.0*GammaB(ikx,zsrclayfort));
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			temp1 = factor1*Pdownmin(ikx)*Wdown(ikx);
			temp2 = factor2*Pdownplusbar(ikx)*Wdownbar(ikx);
		else
			temp1 = factor1*Pdownmin(ikx)*(Wdown(ikx)&
				-Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
			temp2 = factor2*Pdownplusbar(ikx)*(Wdownbar(ikx)&
				+Rpbar(ikx,nlayers)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wupbar(ikx));
		endif
		Ptot1(ikx) = (-1.0)*temp1 - temp2;
		Ptot2(ikx) = (-1.0)*temp1 + temp2;
	end do
endif

END SUBROUTINE Ptotalyy

SUBROUTINE Ptotalyz (Ptot, Pdownplus, Pupplus, Wup, Wdown, Gam, Rp, Rm, etaH, etaV, nk, cor, nz, z, zrcv, zsrc, &
	zsrclay, zrcvlay, gamA, nlayers, above, xdirect);
! Calculate total field for the 23 component (crossline oriented electric receiver, vertical electric source)
!
! Calling arguments:
! Ptot:         Total field, complex array of size nk
! Pdownplus:    TM-mode part of the downgoing field, complex array of size nk
! Pupplus:      TM-mode part of the upgoing field, complex array of size nk
! Wup:          TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:        TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:          Vertical wavenumber Gamma, complex array of size nk
! Rp:           TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:           TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:         Material parameter eta for the horizontal direction, complex array of size nz
! etaV:         Material parameter eta for the vertical direction, complex array of size nz
! nk:           Amount of wavenumbers, integer scalar
! cor:          Coordinates in the wavenumber domain, real array of size nk
! nz:           Amount of layers, integer scalar
! z:            Depth of interfaces, real array of size nz
! zrcv:         Depth of receivers, real scalar
! zsrc:         Depth of source, real scalar
! zsrclay:      The number of the layer where the source is located, integer scalar
! zrcvlay:      The number of the layer where the receivers are located, integer scala
! gamA:         Small gamma squared (zetaH*etaV), complex array of size nz
! nlayers:      Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:        Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!               in the same layer (above=0), integer scalar
! xdirect:      Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), Pupplus(nk), Pdownplus(nk), Wup(nk), Wdown(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nlayers), Rm(nk,nlayers), etaH(nz), etaV(nz), gamA(nz)
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
		factor = 1.0/(2.0*etaV(zsrclayfort));
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupplus(ikx)*Wup(ikx));
			else
				Ptot(ikx) = factor*(Pupplus(ikx)*Wup(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(-1.0*Pdownplus(ikx)*Wdown(ikx));
			else
				Ptot(ikx) = factor*(-1.0*Pdownplus(ikx)*Wdown(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor*(Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx));
			else
				Ptot(ikx) = factor*(Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx)&
					&-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		endif
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor = etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(2.0*etaH(zrcvlayfort)*etaV(zsrclayfort)*Gam(ikx,zsrclayfort));
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			Ptot(ikx) = factor*Pupplus(ikx)*Wup(ikx);
		else
			Ptot(ikx) = factor*Pupplus(ikx)*(Wup(ikx)-Rm(ikx,1)&
				*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
		endif
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor = -1.0*etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(2.0*etaH(zrcvlayfort)*etaV(zsrclayfort)*Gam(ikx,zsrclayfort));
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			Ptot(ikx) = factor*Pdownplus(ikx)*Wdown(ikx);
		else
			Ptot(ikx) = factor*Pdownplus(ikx)*(Wdown(ikx)-Rp(ikx,nlayers)&
				*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
		endif
	end do
endif

END SUBROUTINE Ptotalyz
