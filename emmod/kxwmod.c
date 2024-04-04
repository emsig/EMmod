#include "emmod.h"

void kxwmod(complex *cdata, REAL freq, int nx, REAL dx, int ny, REAL dy, int nz, REAL *z, REAL *econdV, REAL *econdH, REAL *epermV, REAL *epermH, REAL *mpermV, REAL *mpermH, REAL zsrc, REAL zrcv, int zrcv_layer, int zsrc_layer, int above, int *component, int nd, REAL kmax, REAL startlogx, REAL deltalogx, int nlogx, REAL c1, REAL c2, int maxpt, int dopchip, int fullspace, int xdirect, int verbose)
{
	int 	*marker, *markernew, ikx, iky, iz, nlayers, temp, ind, ipos, corel, nlogxnew, nlogxdata, newit, ixs, itcount, component1, component2, componentold, msource, nocomp, nxh, nyh;
	REAL	om, eperm0, mperm0;
	REAL    kintel, a, b;
	REAL    *cor, *tempcor, *xpos, *xposnew;
	complex *Gamma, *GammaB, comptemp;
	complex *etaV, *etaH, *zetaV, *zetaH, *gamA, *gamB;
	complex *Rp, *Rpbar, *Rm, *Rmbar, *Pupplus, *Pupplusbar, *Pupmin, *Pupminbar, *Pdownplus, *Pdownplusbar, *Pdownmin, *Pdownminbar, *Wup, *Wupbar, *Wdown, *Wdownbar, *temptot, *temptot1, *temptot2, *xtotrad1, *xtotrad2, *xtotrad1new, *xtotrad2new, *Ptot, *Gin; 

	/* First only one quadrant is computed and subsequently copied to the other quadrants. 
           The variables nxh and nyh contain half plus one element of the size of the complete setup. */
	nxh = nx/2+1;
	nyh = ny/2+1;

	/* The variable nlayers contains the amount of layers between the receiver and the source
           including the receiver and the source layer. */
	nlayers = abs(zsrc_layer - zrcv_layer) + 1;
	if (verbose) vmess("Number of layers between ind ncluding source and receiver layers (nlayers) = %d", nlayers);

	/* Determine if source is electric or magnetic. */
	component1 = floor(component[0]/10);
	component2 = component[0] - (component1)*10;
	if (component2 >= 4 && component2 <=6) {
		msource = 1;
		if (verbose) vmess("The source is magnetic.");
	} else {
		msource = 0;
		if (verbose) vmess("The source is electric.");
	}

	/* If the source is magnetic, switch the components using reciprocity. */
	if (msource == 1) {
		componentold = component[0];
		if (component1 >= 1 && component1 <= 3){
			component[0] = (component1+3)*10;
		} else {
			component[0] = (component1-3)*10;
		}
		if (component2 >= 1 && component2 <= 3) {
			component[0] = component[0] + component2+3;
		} else {
			component[0] = component[0] + component2-3;
		}
		if (verbose) vmess("Switching from component %d to component %d.", componentold, component[0]);
	}

	/* Initialize the variable nocomp. 
	   It will be set to 1 if the entered component does not exist. */
	nocomp = 0;

	/* Calculate eta, zeta and gamma for each layer. */
	etaV = (complex *)calloc(nz,sizeof(complex));
	etaH = (complex *)calloc(nz,sizeof(complex));
	zetaV = (complex *)calloc(nz,sizeof(complex));
	zetaH = (complex *)calloc(nz,sizeof(complex));
	gamA = (complex *)calloc(nz,sizeof(complex));
	gamB = (complex *)calloc(nz,sizeof(complex));
	
	om = freq*2.0*M_PI; // angular frequency
	eperm0 = 8.854187817e-12; // electric permittivity in free space
	mperm0 = 4.0e-7*M_PI; // magnetic permeability in free space

	// These relations are given in the introduction section
	for (iz=0; iz<nz; iz++) {
		etaH[iz].r = econdH[iz];
		etaH[iz].i = om*epermH[iz]*eperm0;
		etaV[iz].r = econdV[iz];
		etaV[iz].i = om*epermV[iz]*eperm0;
		zetaH[iz].r = 0.0;
		zetaH[iz].i = om*mpermH[iz]*mperm0;
		zetaV[iz].r = 0.0;
		zetaV[iz].i = om*mpermV[iz]*mperm0;
	}

	/* If the source is magnetic, switch eta and zeta because 
	   reciprocity is used to computed the EM-field due to a magnetic source. */
	if (msource == 1) {
		for (iz=0; iz<nz; iz++) {
			comptemp.r = etaH[iz].r;
			comptemp.i = etaH[iz].i;
			etaH[iz].r = -1.0*zetaH[iz].r;
			etaH[iz].i = -1.0*zetaH[iz].i;
			zetaH[iz].r = -1.0*comptemp.r;
			zetaH[iz].i = -1.0*comptemp.i;
			comptemp.r = etaV[iz].r;
			comptemp.i = etaV[iz].i;
			etaV[iz].r = -1.0*zetaV[iz].r;
			etaV[iz].i = -1.0*zetaV[iz].i;
			zetaV[iz].r = -1.0*comptemp.r;
			zetaV[iz].i = -1.0*comptemp.i;
		}
	}
	
	for (iz=0; iz<nz; iz++) {
		/* Note: gamA and gamB are given as gamA squared and gamB squared, respectively, in the introduction. */
		gamA[iz].r = zetaH[iz].r*etaV[iz].r - zetaH[iz].i*etaV[iz].i;
		gamA[iz].i = zetaH[iz].r*etaV[iz].i + zetaH[iz].i*etaV[iz].r;
		gamB[iz].r = zetaV[iz].r*etaH[iz].r - zetaV[iz].i*etaH[iz].i;
		gamB[iz].i = zetaV[iz].r*etaH[iz].i + zetaV[iz].i*etaH[iz].r;
	}

	/* Determine on which points in the wavenumber domain the EM-field needs to be computed
           in order to integrate it with the 61-point Gauss-Kronrod rule (zqk61n.f). */
	corel = 61*nd;
	cor = (REAL *)calloc(corel,sizeof(REAL));
	tempcor = (REAL *)calloc(61,sizeof(REAL));
	kintel = kmax/nd; 
	/* 61 points are necessary for each of the nd subdomains. */
	for (ind=0; ind<nd; ind++) {
		a = ind*kintel; // lower integration boundary
		b = (ind+1)*kintel; // upper integration boundary
		zqk61_setup_grid_(tempcor, &a, &b);
		for (ipos=0; ipos<61; ipos++) {
			cor[ind*61+ipos] = tempcor[ipos];
		}
	}

	/* The computation of the different components begins. */
	if (component[0] == 11) {
		/* This component has slightly more comments than the other components. 
		   It serves as an example, because all other components work similarly. */
		if (verbose) vmess("Computing now component %d", component[0]);

		/* If the medium is layered (fullspace=0) the EM-field needs to be 
		   computed in the wavenumber domain. */
		if (fullspace == 0) {
			/* Allocate memory */
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupmin = (complex *)calloc(corel,sizeof(complex));
			Pdownmin = (complex *)calloc(corel,sizeof(complex));
			Pupplusbar = (complex *)calloc(corel,sizeof(complex));
			Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
			temptot1 = (complex *)calloc(corel,sizeof(complex));
			temptot2 = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));
		
			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			/* compute field propagaters in the layer where the source and the receivers are located as in eq. 74 */
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			/* compute the reflection coming from below the receiver R+ => Pu using equations 64, 65, A-11 and A-12 */
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the reflection coming from above the receiver R- => Pd using equations 64, 65, A-11 and A-12 */
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming below the receiver Pu using equations 81, 95 and 96 */
			pupmin_(Pupmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming above the receiver Pd using equations 82, 103 and 104 */
			pdownmin_(Pdownmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming below the receiver Pu using equations 81, 95 and 96 */
			pupplus_(Pupplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming above the receiver Pd using equations 82, 103 and 104 */
			pdownplus_(Pdownplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);

			/* compute total field using equations B-1, B-13 and B-16 */
			ptotalxx_(temptot1, temptot2, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			/* In the while-loop the points in space for which the Hankel-transformation has to be carried out
			   are determined. The loop runs until no more points need to be added for the specified precision 
			   (newit=0) or until the maximum amount of points (maxpt) has been reached. */
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform using equation 105 (cosine is used in function gridit_xx.F90) */
				hankeltransxx_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
		
				free(marker);
				/* The variable nlogxnew has the size of the maximum amount of points after this iteration
				   is completed. That is nlogx points from the previous iteration plus nlogx-1 points that 
				   could be added (In the worst case, between every point a new datapoint is required.). */
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				xtotrad2new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) { // in case pchip interpolation has been chosen
					evalpoints_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				} else { // in case linear interpolation has been chosen
					evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/

				free(xpos);
				free(xtotrad1);
				free(xtotrad2);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						xtotrad2[ixs].r = xtotrad2new[ikx].r;
						xtotrad2[ixs].i = xtotrad2new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxx_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				free(xtotrad2new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data using the cosine term from equation 105
			if (dopchip==1) { // in case pchip interpolation has been chosen
				gridit_xx_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else { // in case linear interpolation has been chosen
				gridit_xx_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		/* If the medium is a homogeneous fullspace (fullspace=1) 
		   the EM-field is computed directly in the space domain. */
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xx_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		/* Uncomment this piece of code to write intermediate results to the harddisk. */
		/*{
			int nwrite;
			FILE *fp_out = fopen("intermediate.bin","w");
			nwrite = fwrite( &Ptot[0].r, sizeof(complex), nxh*nyh, fp_out);
			assert(nwrite == nxh*nyh);
			fclose(fp_out);
		}*/

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		/* Free memory */
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupmin);
			free(Pdownmin);
			free(Pupplusbar);
			free(Pdownplusbar);
		}
	}else if (component[0] == 12 || component[0] == 21) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupmin = (complex *)calloc(corel,sizeof(complex));
			Pdownmin = (complex *)calloc(corel,sizeof(complex));
			Pupplusbar = (complex *)calloc(corel,sizeof(complex));
			Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));
	
			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}
	
			/* compute field propagaters in the layer where the source and the receivers are located */
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			/* compute the reflection coming from below the receiver R+ => Pu */
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the reflection coming from above the receiver R- => Pd */
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming below the receiver Pu */
			pupmin_(Pupmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming above the receiver Pd */
			pdownmin_(Pdownmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming below the receiver Pu */
			pupplus_(Pupplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			/* compute the field at receiver level coming above the receiver Pd */
			pdownplus_(Pdownplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);

			// compute total field
			ptotalxy_(temptot, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransxy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_xy_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_xy_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xy_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupmin);
			free(Pdownmin);
			free(Pupplusbar);
			free(Pdownplusbar);
		}
	}else if (component[0] == 13) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalxz_(temptot, Pdownplus, Pupplus, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransxz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_xz_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_xz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xz_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupplus);
			free(Pdownplus);
		}
	}else if (component[0] == 22) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupmin = (complex *)calloc(corel,sizeof(complex));
			Pdownmin = (complex *)calloc(corel,sizeof(complex));
			Pupplusbar = (complex *)calloc(corel,sizeof(complex));
			Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
			temptot1 = (complex *)calloc(corel,sizeof(complex));
			temptot2 = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}
	
			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field
			ptotalyy_(temptot1, temptot2, Pdownmin, Pupmin, Pdownplusbar, Pupplusbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransyy_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				xtotrad2new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				free(xtotrad2);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						xtotrad2[ixs].r = xtotrad2new[ikx].r;
						xtotrad2[ixs].i = xtotrad2new[ikx].i;
						ixs++;
					}
			 	}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransyy_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				free(xtotrad2new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_yy_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_yy_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_yy_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupmin);
			free(Pdownmin);
			free(Pupplusbar);
			free(Pdownplusbar);
		}
	}else if (component[0] == 23) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalyz_(temptot, Pdownplus, Pupplus, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransyz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransyz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_yz_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_yz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_yz_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupplus);
			free(Pdownplus);
		}
	}else if (component[0] == 31) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupmin = (complex *)calloc(corel,sizeof(complex));
			Pdownmin = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalzx_(temptot, Pdownmin, Pupmin, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltranszx_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltranszx_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_zx_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_zx_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_zx_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer]);
		}
	
		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupmin);
			free(Pdownmin);
		}
	}else if (component[0] == 32) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupmin = (complex *)calloc(corel,sizeof(complex));
			Pdownmin = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownmin, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalzy_(temptot, Pdownmin, Pupmin, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltranszy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltranszy_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_zy_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_zy_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_zy_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupmin);
			free(Pdownmin);
		}
	}else if (component[0] == 33) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdown + Pup 
			ptotalzz_(temptot, Pdownplus, Pupplus, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamA, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltranszz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);

				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltranszz_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
	
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_zz_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_zz_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_zz_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupplus);
			free(Pdownplus);
		}
	}else if (component[0] == 41) {
		//The implementation is G_{xx}^{em}
		//Interchange of source and receiver position and a minus sign in ptotalxxm leads to G_{xx}^{me}
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			Pupminbar = (complex *)calloc(corel,sizeof(complex));
			Pdownminbar = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalxxm_(temptot, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);
		
			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransxxm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxxm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_xxm_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_xxm_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xxm_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupplus);
			free(Pdownplus);
			free(Pupminbar);
			free(Pdownminbar);
		}
	}else if (component[0] == 42) {
		//The implementation is G_{xx}^{em}
		//Interchange of source and receiver position and a minus sign in ptotalxxm leads to G_{xx}^{me}
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			Pupminbar = (complex *)calloc(corel,sizeof(complex));
			Pdownminbar = (complex *)calloc(corel,sizeof(complex));
			temptot1 = (complex *)calloc(corel,sizeof(complex));
			temptot2 = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalxym_(temptot1, temptot2, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransxym_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				xtotrad2new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				free(xtotrad2);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						xtotrad2[ixs].r = xtotrad2new[ikx].r;
						xtotrad2[ixs].i = xtotrad2new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxym_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				free(xtotrad2new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_xym_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_xym_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xym_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}
		
		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupplus);
			free(Pdownplus);
			free(Pupminbar);
			free(Pdownminbar);
		}
	}else if (component[0] == 43) {
		//The implementation is G_{xx}^{em}
		//Interchange of source and receiver position and a minus sign in ptotalxxm leads to G_{xx}^{me}
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalxzm_(temptot, Pdownplus, Pupplus, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransxzm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransxzm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_xzm_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_xzm_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_xzm_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupplus);
			free(Pdownplus);
		}
	}else if (component[0] == 51) {
		//The implementation is G_{xx}^{em}
		//Interchange of source and receiver position and a minus sign in ptotalxxm leads to G_{xx}^{me}
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			Pupminbar = (complex *)calloc(corel,sizeof(complex));
			Pdownminbar = (complex *)calloc(corel,sizeof(complex));
			temptot1 = (complex *)calloc(corel,sizeof(complex));
			temptot2 = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalyxm_(temptot1, temptot2, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransyxm_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				xtotrad2new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_(xposnew,xtotrad1new,xtotrad2new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,xtotrad2,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				free(xtotrad2);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				xtotrad2 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						xtotrad2[ixs].r = xtotrad2new[ikx].r;
						xtotrad2[ixs].i = xtotrad2new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransyxm_(xtotrad1,xtotrad2,marker,temptot1,temptot2,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				free(xtotrad2new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_yxm_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_yxm_lin_(Ptot, xtotrad1, xtotrad2, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_yxm_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupplus);
			free(Pdownplus);
			free(Pupminbar);
			free(Pdownminbar);
		}
	}else if (component[0] == 52) {
		//The implementation is G_{yy}^{em}
		//Interchange of source and receiver position
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			Pupminbar = (complex *)calloc(corel,sizeof(complex));
			Pdownminbar = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupmin_(Pupminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownmin_(Pdownminbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalyym_(temptot, Pdownplus, Pupplus, Pdownminbar, Pupminbar, Wup, Wdown, Wupbar, Wdownbar, Gamma, GammaB, Rp, Rpbar, Rm, Rmbar, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransyym_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransyym_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_yym_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_yym_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_yym_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamA[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wdownbar);
			free(Wup);
			free(Wupbar);
			free(Rp);
			free(Rpbar);
			free(Rm);
			free(Rmbar);
			free(Pupplus);
			free(Pdownplus);
			free(Pupminbar);
			free(Pdownminbar);
		}
	}else if (component[0] == 53) {
		//The implementation is G_{xx}^{em}
		//Interchange of source and receiver position and a minus sign in ptotalxxm leads to G_{xx}^{me}
		above = -above;
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdown = (complex *)calloc(corel,sizeof(complex));
			Wup = (complex *)calloc(corel,sizeof(complex));
			Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rm = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplus = (complex *)calloc(corel,sizeof(complex));
			Pdownplus = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			Gamma = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zsrc, &zsrc_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rm, &corel, cor, etaH, Gamma, &nz, z, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplus, Rp, Rm, Gamma, &corel, cor, &nz, z, &zrcv, &zsrc_layer, &zrcv_layer, &nlayers, &above);
			// compute total field
			ptotalyzm_(temptot, Pdownplus, Pupplus, Wup, Wdown, Gamma, Rp, Rm, etaH, etaV, zetaH, zetaV, &corel, cor, &nz, z, &zsrc, &zrcv, &zrcv_layer, &zsrc_layer, gamA, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltransyzm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltransyzm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_yzm_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			} else {
				gridit_yzm_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_yzm_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &etaH[zsrc_layer], &etaV[zsrc_layer], &zetaH[zsrc_layer], &gamA[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdown);
			free(Wup);
			free(Rp);
			free(Rm);
			free(Pupplus);
			free(Pdownplus);
		}
	}else if (component[0] == 61) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplusbar = (complex *)calloc(corel,sizeof(complex));
			Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));
	
			for (iz=0; iz<nz; iz++) {
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalzxm_(temptot, Pdownplusbar, Pupplusbar, Wupbar, Wdownbar, GammaB, Rpbar, Rmbar, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltranszxm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltranszxm_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_zxm_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_zxm_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_zxm_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = -cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdownbar);
			free(Wupbar);
			free(Rpbar);
			free(Rmbar);
			free(Pupplusbar);
			free(Pdownplusbar);
		}
	}else if (component[0] == 62) {
		if (verbose) vmess("Computing now component %d", component[0]);
		
		if (fullspace == 0) {
			// Allocate memory
			Wdownbar = (complex *)calloc(corel,sizeof(complex));
			Wupbar = (complex *)calloc(corel,sizeof(complex));
			Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Rmbar = (complex *)calloc(corel*nlayers,sizeof(complex));
			Pupplusbar = (complex *)calloc(corel,sizeof(complex));
			Pdownplusbar = (complex *)calloc(corel,sizeof(complex));
			temptot = (complex *)calloc(corel,sizeof(complex));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
			GammaB = (complex *)calloc(nz*corel,sizeof(complex));

			for (iz=0; iz<nz; iz++) {
				gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
			}

			// compute field propagaters in the layer where the source and the receivers are located
			wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
			// compute the reflection coming from below the receiver R+ => Pu
			rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the reflection coming from above the receiver R- => Pd
			rmin_(Rmbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming below the receiver Pu
			pupplus_(Pupplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute the field at receiver level coming above the receiver Pd
			pdownplus_(Pdownplusbar, Rpbar, Rmbar, GammaB, &corel, cor, &nz, z, &zsrc, &zrcv_layer, &zsrc_layer, &nlayers, &above);
			// compute total field: incident + Pdownmin + Pupmin
			ptotalzym_(temptot, Pdownplusbar, Pupplusbar, Wupbar, Wdownbar, GammaB, Rpbar, Rmbar, zetaH, zetaV, &corel, cor, &nz, z, &zrcv, &zsrc, &zsrc_layer, &zrcv_layer, gamB, &nlayers, &above, &xdirect);

			/* Set up coordinate vector*/
			getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

			/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
			newit = 1;
			/* The Hankel transformation will be evaluated where the marker is 2.*/
			marker = (int *)calloc(nlogx,sizeof(int));
			for (ikx=0; ikx<nlogx; ikx++) {
				marker[ikx] = 2;
			}
			itcount = 0;
			while ((newit == 1) && (nlogx < maxpt)) {
				if (verbose) vmess("Iteration %d", itcount);
				/* Compute the Hankel Transform */
				hankeltranszym_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				free(marker);
				/* In the worst case, between every point a new datapoint is required.*/
				nlogxnew = 2*nlogx-1;
				markernew = (int *)calloc(nlogxnew,sizeof(int));
				xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
				xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
				/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
				nlogxdata = 0;
				/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
				if (dopchip==1) {
					evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				} else {
					evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
				}
				if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
				/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
				The marker is set to 2 at these locations.*/
				free(xpos);
				free(xtotrad1);
				nlogx = nlogxdata;
				marker = (int *)calloc(nlogx,sizeof(int));
				xpos = (REAL *)calloc(nlogx,sizeof(REAL));
				xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
				ixs = 0;
				for (ikx=0; ikx<nlogxnew; ikx++) {
					if (markernew[ikx]!=0) { 
						marker[ixs] = markernew[ikx];
						xpos[ixs] = xposnew[ikx];
						xtotrad1[ixs].r = xtotrad1new[ikx].r;
						xtotrad1[ixs].i = xtotrad1new[ikx].i;
						ixs++;
					}
				}
				if (nlogx>maxpt) {
					if (verbose) vwarn("More points are needed to achieve required precision,");
					if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
					hankeltranszym_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
				}
				free(markernew);
				free(xposnew);
				free(xtotrad1new);
				itcount = itcount + 1;
			}
			if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

			// Compute the field for one quadrant based on the radial data
			if (dopchip==1) {
				gridit_zym_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			} else {
				gridit_zym_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer], &above, &xdirect);
			}
		} else {
			Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
			gridit_zym_fullspace_(Ptot, &dx, &nxh, &dy, &nyh, &zrcv, &zsrc, &zetaH[zsrc_layer], &zetaV[zsrc_layer], &gamB[zsrc_layer]);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = -Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		if (fullspace == 0) {
			free(Wdownbar);
			free(Wupbar);
			free(Rpbar);
			free(Rmbar);
			free(Pupplusbar);
			free(Pdownplusbar);
		}
	}else if (component[0] == 63) {
		if (verbose) vmess("Computing now component %d", component[0]);
		// Allocate memory
		Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		// compute total field: incident + Pdownmin + Pupmin
		ptotalzzm_(Ptot, &nxh, &nyh);

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

	}else if (component[0] == 77) {
		if (verbose) vmess("Computing now component %d: TM reflection response", component[0]);
		// Allocate memory
		Wup = (complex *)calloc(corel,sizeof(complex));
		Wdown = (complex *)calloc(corel,sizeof(complex));
		Rp = (complex *)calloc(corel*nlayers,sizeof(complex));
		temptot = (complex *)calloc(corel,sizeof(complex));
		xpos = (REAL *)calloc(nlogx,sizeof(REAL));
		xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
		Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));
		
		Gamma = (complex *)calloc(nz*corel,sizeof(complex));
		
		for (iz=0; iz<nz; iz++) {
			gammarad_(&Gamma[iz*corel], &corel, cor, &etaV[iz], &etaH[iz], &gamA[iz]);
		}

		// compute field propagaters in the layer where the source and the receivers are located
		wprop_(Wup, Wdown, Gamma, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
		// compute the reflection coming from below the receiver R+ => Pu
		temp = 0;
		rplus_(Rp, &corel, cor, etaH, Gamma, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &temp);
		// compute total field
		ptotalref_(temptot, Wup, Rp, &corel);

		/* Set up coordinate vector*/
		getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

		/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
		newit = 1;
		/* The Hankel transformation will be evaluated where the marker is 2.*/
		marker = (int *)calloc(nlogx,sizeof(int));
		for (ikx=0; ikx<nlogx; ikx++) {
			marker[ikx] = 2;
		}
		itcount = 0;
		while ((newit == 1) && (nlogx < maxpt)) {
			if (verbose) vmess("Iteration %d", itcount);
			/* Compute the Hankel Transform */
			hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
			free(marker);
			/* In the worst case, between every point a new datapoint is required.*/
			nlogxnew = 2*nlogx-1;
			markernew = (int *)calloc(nlogxnew,sizeof(int));
			xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
			xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
			/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
			nlogxdata = 0;
			/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
			if (dopchip==1) {
				evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
			} else {
				evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
			}
			if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
			/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
			The marker is set to 2 at these locations.*/
			free(xpos);
			free(xtotrad1);
			nlogx = nlogxdata;
			marker = (int *)calloc(nlogx,sizeof(int));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			ixs = 0;
			for (ikx=0; ikx<nlogxnew; ikx++) {
				if (markernew[ikx]!=0) { 
					marker[ixs] = markernew[ikx];
					xpos[ixs] = xposnew[ikx];
					xtotrad1[ixs].r = xtotrad1new[ikx].r;
					xtotrad1[ixs].i = xtotrad1new[ikx].i;
					ixs++;
				}
			}
			if (nlogx>maxpt) {
				if (verbose) vwarn("More points are needed to achieve required precision,");
				if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
				hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
			} else {
			}
			free(markernew);
			free(xposnew);
			free(xtotrad1new);
			itcount = itcount + 1;
		}
		if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

		// Compute the field for one quadrant based on the radial data
		if (dopchip==1) {
			gridit_ref_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
		} else {
			gridit_ref_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}

		// Free memory
		free(Wup);
		free(Wdown);
		free(Rp);
	}else if (component[0] == 88) {
		if (verbose) vmess("Computing now component %d: TE reflection response", component[0]);
		// Allocate memory
		Wupbar = (complex *)calloc(corel,sizeof(complex));
		Wdownbar = (complex *)calloc(corel,sizeof(complex));
		Rpbar = (complex *)calloc(corel*nlayers,sizeof(complex));
		temptot = (complex *)calloc(corel,sizeof(complex));
		xpos = (REAL *)calloc(nlogx,sizeof(REAL));
		xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
		Ptot = (complex *)calloc(nxh*nyh,sizeof(complex));

		GammaB = (complex *)calloc(nz*corel,sizeof(complex));

		for (iz=0; iz<nz; iz++) {
			gammarad_(&GammaB[iz*corel], &corel, cor, &zetaV[iz], &zetaH[iz], &gamB[iz]);
		}

		// compute field propagaters in the layer where the source and the receivers are located
		wprop_(Wupbar, Wdownbar, GammaB, &corel, cor, &nz, z, &zrcv, &zrcv_layer);
		// compute the reflection coming from below the receiver R+ => Pu
		temp = 0;
		rplus_(Rpbar, &corel, cor, zetaH, GammaB, &nz, z, &zrcv_layer, &zsrc_layer, &nlayers, &temp);
		// compute total field
		ptotalref_(temptot, Wupbar, Rpbar, &corel);

		/* Set up coordinate vector*/
		getcoords_(xpos,&startlogx,&deltalogx,&nlogx);

		/* If newit is set to 0 by the evalpoints subroutine no more points are added to the coordinate vector and the optimization process stops.*/
		newit = 1;
		/* The Hankel transformation will be evaluated where the marker is 2.*/
		marker = (int *)calloc(nlogx,sizeof(int));
		for (ikx=0; ikx<nlogx; ikx++) {
			marker[ikx] = 2;
		}
		itcount = 0;
		while ((newit == 1) && (nlogx < maxpt)) {
			if (verbose) vmess("Iteration %d", itcount);
			/* Compute the Hankel Transform */
			hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
			free(marker);
			/* In the worst case, between every point a new datapoint is required.*/
			nlogxnew = 2*nlogx-1;
			markernew = (int *)calloc(nlogxnew,sizeof(int));
			xposnew = (REAL *)calloc(nlogxnew,sizeof(REAL));
			xtotrad1new = (complex *)calloc(nlogxnew,sizeof(complex));
			/* nlogxdata will show how many of the nlogxnew elements actually contain data*/
			nlogxdata = 0;
			/* Evaluate datapoints, decide on which points the Hankel transformation still needs to be computed*/
			if (dopchip==1) {
				evalpoints_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
			} else {
				evalpoints_lin_mono_(xposnew,xtotrad1new,markernew,&nlogxnew,&nlogxdata,&newit,xpos,xtotrad1,&nlogx,&c1,&c2);
			}
			if (verbose) vmess("Added %d new datapoints out of %d possible datapoints.",nlogxdata-nlogx,nlogxnew-nlogx);
			/* Set up the new vectors which contain the new datapoints. xtotrad is 0.0 where a new datapoint should be.
			The marker is set to 2 at these locations.*/
			free(xpos);
			free(xtotrad1);
			nlogx = nlogxdata;
			marker = (int *)calloc(nlogx,sizeof(int));
			xpos = (REAL *)calloc(nlogx,sizeof(REAL));
			xtotrad1 = (complex *)calloc(nlogx,sizeof(complex));
			ixs = 0;
			for (ikx=0; ikx<nlogxnew; ikx++) {
				if (markernew[ikx]!=0) { 
					marker[ixs] = markernew[ikx];
					xpos[ixs] = xposnew[ikx];
					xtotrad1[ixs].r = xtotrad1new[ikx].r;
					xtotrad1[ixs].i = xtotrad1new[ikx].i;
					ixs++;
				}
			}
			if (nlogx>maxpt) {
				if (verbose) vwarn("More points are needed to achieve required precision,");
				if (verbose) vwarn("but maximum amount of points (%d) is reached",maxpt);
				hankeltransref_(xtotrad1,marker,temptot,&corel,cor,xpos,&nlogx,&nd,&kmax);
			}
			free(markernew);
			free(xposnew);
			free(xtotrad1new);
			itcount = itcount + 1;
		}
		if (verbose) vmess("Using %d datapoints in the space domain.",nlogx);

		// Compute the field for one quadrant based on the radial data
		if (dopchip==1) {
			gridit_ref_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
		} else {
			gridit_ref_lin_(Ptot, xtotrad1, &dx, &nxh, &dy, &nyh, xpos, &nlogx);
		}

		// copy quadrants to other parts
		for (iky = 0; iky<nyh; iky++) {
			// 1st quadrant
			for (ikx = 0; ikx<nxh; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(nxh-ikx-1)+(nyh-iky-1)*nxh].i;
			}
			// 2nd quadrant
			for (ikx = nxh; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].r;
				cdata[ikx+iky*nx].i = Ptot[(ikx-nxh+1)+(nyh-iky-1)*nxh].i;
			}
		}
		for (iky = nyh; iky<ny; iky++) {
			// 3rd and 4th quadrant
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = cdata[ikx+(ny-iky)*nx].r;
				cdata[ikx+iky*nx].i = cdata[ikx+(ny-iky)*nx].i;
			}
		}
	
		// Free memory
		free(Wupbar);
		free(Wdownbar);
		free(Rpbar);
	}else { // in case the user entered a component that does not exist
		if (verbose) vmess("The component %d is not implemented in the current version of emmod.",component[0]);
		nocomp = 1;
	}

	/* If the source is magnetic, reciprocity is used. In order to get the EM-fields right, 
	   a negative sign needs to be added. */
	if (msource == 1 && nocomp != 1) {
		for (iky = 0; iky<ny; iky++) {
			for (ikx = 0; ikx<nx; ikx++) {
				cdata[ikx+iky*nx].r = -1.0*cdata[ikx+iky*nx].r;
				cdata[ikx+iky*nx].i = -1.0*cdata[ikx+iky*nx].i;
			}
		}
	}

	/* Free memory */
	if (nocomp != 1) { // if the component entered exists
		if (fullspace == 0) { // if the model is layered
			if (component[0] == 31 || component[0] == 13 || component[0] == 32 || component[0] == 33 || component[0] == 23 || component[0] == 43 || component[0] == 53 || component[0] == 77) { // components containing only the TM-mode
				free(Gamma);
				free(Ptot);
			}else if (component[0] == 61 || component[0] == 62 || component[0] == 88) { // components containing onlhy the TE-mode
				free(GammaB);
				free(Ptot);
			}else if (component[0] == 11 || component[0] == 12 || component[0] == 21 || component[0] == 22 || component[0] == 41 || component[0] == 42 || component[0] == 51 || component[0] == 52) { // components containing TE- as well as TM-modes
				free(Gamma);
				free(GammaB);
				free(Ptot);
			}
		} else { // if the model is a homogeneous fullspace
			free(Ptot);
		}
	}
	
	return;
}
