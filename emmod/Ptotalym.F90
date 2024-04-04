SUBROUTINE Ptotalyxm (Ptot1, Ptot2, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, &
	Gam, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, nk, cor, nz, z, zrcv, zsrc, &
	zsrclay, zrcvlay, gamA, gamB, nlayers, above, xdirect);
! Calculate total field for the 51 component (crossline oriented magnetic receiver, inline oriented electric source)
!
! Calling arguments:
! Ptot1:       1st term of the total field, complex array of size nk
! Ptot2:       2nd term of the total field, complex array of size nk
! Pdownplus:   TM-mode part of the downgoing field, complex array of size nk
! Pupplus:     TM-mode part of the upgoing field, complex array of size nk
! Pdownminbar: TE-mode part of the downgoing field, complex array of size nk
! Pupminbar:   TE-mode part of the upgoing field, complex array of size nk
! Wup:         TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:       TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Wupbar:      TE-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdownbar:    TE-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:         Vertical wavenumber Gamma, complex array of size nk
! GammaB:      Vertical wavenumber Gamma Bar, complex array of size nk
! Rp:          TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rpbar:       TE-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:          TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! Rmbar:       TE-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:        Material parameter eta for the horizontal direction, complex array of size nz
! etaV:        Material parameter eta for the vertical direction, complex array of size nz
! zetaH:       Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:       Material parameter zeta for the vertical direction, complex array of size nz
! nk:          Amount of wavenumbers, integer scalar
! cor:         Coordinates in the wavenumber domain, real array of size nk
! nz:          Amount of layers, integer scalar
! z:           Depth of interfaces, real array of size nz
! zrcv:        Depth of receivers, real scalar
! zsrc:        Depth of source, real scalar
! zsrclay:     The number of the layer where the source is located, integer scalar
! zrcvlay:     The number of the layer where the receivers are located, integer scala
! gamA:        Small gamma squared (zetaH*etaV), complex array of size nz
! gamB:        Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:     Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:       Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!              in the same layer (above=0), integer scalar
! xdirect:     Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot1(nk), Ptot2(nk)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), GammaB(nk,nz)
COMPLEX,INTENT(IN)    :: Pdownplus(nk), Pupplus(nk), Pdownminbar(nk), Pupminbar(nk)
COMPLEX,INTENT(IN)    :: Wup(nk), Wdown(nk), Wupbar(nk), Wdownbar(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nlayers), Rpbar(nk,nlayers), Rm(nk,nlayers), Rmbar(nk,nlayers) 
COMPLEX,INTENT(IN)    :: etaH(nz), etaV(nz), zetaH(nz), zetaV(nz), gamA(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
COMPLEX               :: temp1, temp2
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor1, factor2

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
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = Pupplus(ikx)*Wup(ikx);
				temp2 = Pupminbar(ikx)*Wupbar(ikx);
			else
				temp1 = Pupplus(ikx)*Wup(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = -1.0*Pdownplus(ikx)*Wdown(ikx);
				temp2 = -1.0*Pdownminbar(ikx)*Wdownbar(ikx);
			else
				temp1 = -1.0*Pdownplus(ikx)*Wdown(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = -1.0*Pdownminbar(ikx)*Wdownbar(ikx)+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)-Pdownminbar(ikx)*Wdownbar(ikx);
			else
				temp1 = Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx)&
					&-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)-Pdownminbar(ikx)*Wdownbar(ikx)&
					&+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		endif
		Ptot1(ikx) = temp1 - temp2;
		Ptot2(ikx) = temp1 + temp2;
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor1 = etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(etaH(zrcvlayfort)*Gam(ikx,zsrclayfort));
		factor2 = 1.0;
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			temp1 = factor1*Pupplus(ikx)*Wup(ikx);
			temp2 = factor2*Pupminbar(ikx)*Wupbar(ikx);
		else
			temp1 = factor1*Pupplus(ikx)*(Wup(ikx)&
				-Rm(ikx,1)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
			temp2 = factor2*Pupminbar(ikx)*(Wupbar(ikx)&
				+Rmbar(ikx,1)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdownbar(ikx));
		endif
		Ptot1(ikx) = temp1 - temp2;
		Ptot2(ikx) = temp1 + temp2;
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor1 = -1.0*etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(etaH(zrcvlayfort)*Gam(ikx,zsrclayfort));
		factor2 = 1.0;
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			temp1 = factor1*Pdownplus(ikx)*Wdown(ikx);
			temp2 = factor2*Pdownminbar(ikx)*Wdownbar(ikx);
		else
			temp1 = factor1*Pdownplus(ikx)*(Wdown(ikx)&
				-Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
			temp2 = factor2*Pdownminbar(ikx)*(Wdownbar(ikx)&
				+Rpbar(ikx,nlayers)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wupbar(ikx));
		endif
		Ptot1(ikx) = temp1 - temp2;
		Ptot2(ikx) = temp1 + temp2;
	end do
endif

END SUBROUTINE Ptotalyxm

SUBROUTINE Ptotalyym (Ptot, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, &
	Gam, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, nk, cor, nz, z, zrcv, zsrc, &
	zsrclay, zrcvlay, gamA, gamB, nlayers, above, xdirect);
! Calculate total field for the 52 component (crossline oriented magnetic receiver, crossline oriented electric source)
!
! Calling arguments:
! Ptot:        Total field, complex array of size nk
! Pdownplus:   TM-mode part of the downgoing field, complex array of size nk
! Pupplus:     TM-mode part of the upgoing field, complex array of size nk
! Pdownminbar: TE-mode part of the downgoing field, complex array of size nk
! Pupminbar:   TE-mode part of the upgoing field, complex array of size nk
! Wup:         TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:       TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Wupbar:      TE-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdownbar:    TE-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:         Vertical wavenumber Gamma, complex array of size nk
! GammaB:      Vertical wavenumber Gamma Bar, complex array of size nk
! Rp:          TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rpbar:       TE-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:          TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! Rmbar:       TE-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:        Material parameter eta for the horizontal direction, complex array of size nz
! etaV:        Material parameter eta for the vertical direction, complex array of size nz
! zetaH:       Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:       Material parameter zeta for the vertical direction, complex array of size nz
! nk:          Amount of wavenumbers, integer scalar
! cor:         Coordinates in the wavenumber domain, real array of size nk
! nz:          Amount of layers, integer scalar
! z:           Depth of interfaces, real array of size nz
! zrcv:        Depth of receivers, real scalar
! zsrc:        Depth of source, real scalar
! zsrclay:     The number of the layer where the source is located, integer scalar
! zrcvlay:     The number of the layer where the receivers are located, integer scala
! gamA:        Small gamma squared (zetaH*etaV), complex array of size nz
! gamB:        Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:     Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:       Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!              in the same layer (above=0), integer scalar
! xdirect:     Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: Gam(nk,nz), GammaB(nk,nz)
COMPLEX,INTENT(IN)    :: Pdownplus(nk), Pupplus(nk), Pdownminbar(nk), Pupminbar(nk)
COMPLEX,INTENT(IN)    :: Wup(nk), Wdown(nk), Wupbar(nk), Wdownbar(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nlayers), Rpbar(nk,nlayers), Rm(nk,nlayers), Rmbar(nk,nlayers) 
COMPLEX,INTENT(IN)    :: etaH(nz), etaV(nz), zetaH(nz), zetaV(nz), gamA(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
COMPLEX               :: temp1, temp2
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor1, factor2

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
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = Pupplus(ikx)*Wup(ikx);
				temp2 = Pupminbar(ikx)*Wupbar(ikx);
			else
				temp1 = Pupplus(ikx)*Wup(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = -1.0*Pdownplus(ikx)*Wdown(ikx);
				temp2 = -1.0*Pdownminbar(ikx)*Wdownbar(ikx);
			else
				temp1 = -1.0*Pdownplus(ikx)*Wdown(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = -1.0*Pdownminbar(ikx)*Wdownbar(ikx)+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				temp1 = Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)-Pdownminbar(ikx)*Wdownbar(ikx);
			else
				temp1 = Pupplus(ikx)*Wup(ikx)-Pdownplus(ikx)*Wdown(ikx)-signd*exp(-1.0*Gam(ikx,zsrclayfort)*d);
				temp2 = Pupminbar(ikx)*Wupbar(ikx)-Pdownminbar(ikx)*Wdownbar(ikx)+signd*exp(-1.0*GammaB(ikx,zsrclayfort)*d);
			endif
		endif
		Ptot(ikx) = temp1 + temp2;
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor1 = etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(etaH(zrcvlayfort)*Gam(ikx,zsrclayfort));
		factor2 = 1.0;
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			Ptot(ikx) = factor1*Pupplus(ikx)*Wup(ikx)&
				+factor2*Pupminbar(ikx)*Wupbar(ikx);
		else
			Ptot(ikx) = factor1*Pupplus(ikx)*(Wup(ikx)&
				-Rm(ikx,1)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx))&
				+factor2*Pupminbar(ikx)*(Wupbar(ikx)&
				+Rmbar(ikx,1)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdownbar(ikx));
		endif
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor1 = -1.0*etaH(zsrclayfort)*Gam(ikx,zrcvlayfort)/(etaH(zrcvlayfort)*Gam(ikx,zsrclayfort));
		factor2 = 1.0;
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			Ptot(ikx) = factor1*Pdownplus(ikx)*Wdown(ikx)&
				+factor2*Pdownminbar(ikx)*Wdownbar(ikx);
		else
			Ptot(ikx) = factor1*Pdownplus(ikx)*(Wdown(ikx)&
				-Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx))&
				+factor2*Pdownminbar(ikx)*(Wdownbar(ikx)&
				+Rpbar(ikx,nlayers)*exp(-1.0*GammaB(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wupbar(ikx));
		endif
	end do
endif

END SUBROUTINE Ptotalyym

SUBROUTINE Ptotalyzm (Ptot, Pdownplus, Pupplus, Wup, Wdown, &
	Gam, Rp, Rm, etaH, etaV, zetaH, zetaV, nk, cor, nz, z, zrcv, zsrc, &
	zsrclay, zrcvlay, gamA, gamB, nlayers, above, xdirect);
! Calculate total field for the 53 component (crossline oriented magnetic receiver, vertical electric source)
!
! Calling arguments:
! Ptot:        Total field, complex array of size nk
! Pdownplus:   TM-mode part of the downgoing field, complex array of size nk
! Pupplus:     TM-mode part of the upgoing field, complex array of size nk
! Wup:         TM-mode propagation of the EM-field from the lower interface of the layer to the receiver, complex array of size nk
! Wdown:       TM-mode propagation of the EM-field from the upper interface of the layer to the receiver, complex array of size nk
! Gam:         Vertical wavenumber Gamma, complex array of size nk
! Rp:          TM-mode upgoing global reflection coefficient, complex array of size nk times nlayers
! Rm:          TM-mode downgoing global reflection coefficient, complex array of size nk times nlayers
! etaH:        Material parameter eta for the horizontal direction, complex array of size nz
! etaV:        Material parameter eta for the vertical direction, complex array of size nz
! zetaH:       Material parameter zeta for the horizontal direction, complex array of size nz
! zetaV:       Material parameter zeta for the vertical direction, complex array of size nz
! nk:          Amount of wavenumbers, integer scalar
! cor:         Coordinates in the wavenumber domain, real array of size nk
! nz:          Amount of layers, integer scalar
! z:           Depth of interfaces, real array of size nz
! zrcv:        Depth of receivers, real scalar
! zsrc:        Depth of source, real scalar
! zsrclay:     The number of the layer where the source is located, integer scalar
! zrcvlay:     The number of the layer where the receivers are located, integer scala
! gamA:        Small gamma squared (zetaH*etaV), complex array of size nz
! gamB:        Small gamma bar squared (zetaV*etaH), complex array of size nz
! nlayers:     Amount of layers between the source and the receivers (including source and receiver layer), integer scalar
! above:       Indicates if the receivers are above the source (above=1), below the source (above=-1) or 
!              in the same layer (above=0), integer scalar
! xdirect:     Indicates if the direct field is computed in the space domain (xdirect=1) or in the wavenumber domain (xdirect=0)

IMPLICIT NONE

INTEGER,INTENT(IN)    :: nk, nz, zsrclay, zrcvlay, nlayers, above, xdirect
COMPLEX,INTENT(INOUT) :: Ptot(nk)
COMPLEX,INTENT(IN)    :: Gam(nk,nz)
COMPLEX,INTENT(IN)    :: Pdownplus(nk), Pupplus(nk)
COMPLEX,INTENT(IN)    :: Wup(nk), Wdown(nk)
COMPLEX,INTENT(IN)    :: Rp(nk,nlayers), Rm(nk,nlayers) 
COMPLEX,INTENT(IN)    :: etaH(nz), etaV(nz), zetaH(nz), zetaV(nz), gamA(nz), gamB(nz)
REAL   ,INTENT(IN)    :: cor(nk), z(nz), zsrc, zrcv
REAL                  :: d
INTEGER               :: ikx, zsrclayfort, zrcvlayfort, signd
COMPLEX               :: factor1

zsrclayfort = zsrclay + 1; ! First array index in C is 0, but in Fortran it is 1!
zrcvlayfort = zrcvlay + 1;
d = zrcv-zsrc;
if (d < 0) then
	signd = -1.0;
else 
	signd = 1.0;
endif
if (d < 0) d = -1.0*d; ! replaces abs(d)

! Implemented is G_{zx}^{em}. To get G_{xz}^{me}, source and receivers are interchanged and a minus-sign is applied. 
! To compensate for lateral symmetries, another minus sign needs to be applied. 

if (above == 0) then ! receivers in the source layer
	do ikx=1,nk
		factor1 = etaH(zsrclayfort)/(2.0*etaV(zsrclayfort)*Gam(ikx,zsrclayfort));
		if (zsrclayfort == 1) then ! If the source and the receivers are in the top most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor1*(Pupplus(ikx)*Wup(ikx));
			else
				Ptot(ikx) = factor1*(Pupplus(ikx)*Wup(ikx)+exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		elseif (zsrclayfort == nz) then ! If the source and the receivers are in the lower most layer
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor1*(Pdownplus(ikx)*Wdown(ikx));
			else
				Ptot(ikx) = factor1*(Pdownplus(ikx)*Wdown(ikx)+exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		else ! If the source and the receivers are in any layer in between. 
			if (xdirect == 1) then ! The direct field is computed in the space domain
				Ptot(ikx) = factor1*(Pupplus(ikx)*Wup(ikx)+Pdownplus(ikx)*Wdown(ikx));
			else
				Ptot(ikx) = factor1*(Pupplus(ikx)*Wup(ikx)+Pdownplus(ikx)*Wdown(ikx)&
					&+exp(-1.0*Gam(ikx,zsrclayfort)*d));
			endif
		endif
	end do
elseif (above == 1) then ! receivers above the source layer: Pdownplus is not used
	do ikx=1,nk
		factor1 = etaH(zsrclayfort)/(2.0*etaV(zrcvlayfort)*Gam(ikx,zsrclayfort));
		if (zrcvlayfort == 1) then ! If the receivers are in the top most layer
			Ptot(ikx) = factor1*Pupplus(ikx)*Wup(ikx);
		else
			Ptot(ikx) = factor1*Pupplus(ikx)*(Wup(ikx)&
				+Rm(ikx,1)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wdown(ikx));
		endif
	end do
elseif (above == -1) then ! receivers below the source layer
	do ikx=1,nk
		factor1 = etaH(zsrclayfort)/(2.0*etaV(zrcvlayfort)*Gam(ikx,zsrclayfort));
		if (zrcvlayfort == nz) then ! If the receivers are in the lower most layer
			Ptot(ikx) = factor1*Pdownplus(ikx)*Wdown(ikx);
		else
			Ptot(ikx) = factor1*Pdownplus(ikx)*(Wdown(ikx)&
				+Rp(ikx,nlayers)*exp(-1.0*Gam(ikx,zrcvlayfort)*(z(zrcvlayfort+1)-z(zrcvlayfort)))*Wup(ikx));
		endif
	end do
endif

END SUBROUTINE Ptotalyzm
