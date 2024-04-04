#include "emmod.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" EMmod - modeling EM responses in a layered VTI medium",
" 								",
" Required parameters:",
" ",
"   zsrc ..................... depth of source",
"   zrcv ..................... depth of receiver",
"   nx ....................... number of x-samples",
"   ny ....................... number of y-samples",
"   dx ....................... x-sampling",
"   dy ....................... y-sampling",
"   econdV ................... vertical electric conductivity",
"   econdH ................... horizontal electric conductivity",
"   epermV ................... vertical electric permittivity",
"   epermH ................... horizontal electric permittivity",
"   mpermV ................... vertical magnetic permeability",
"   mpermH ................... horizontal magnetic permeability",
"   z=0  ..................... depth of layer",
" 							        ",
" Optional parameters:",
" ",
"   file_out ................. output file (default EMmod_output.bin)",
"   writebin ................. write output in binary format (1) or ASCI (0)",
"                              (default is 1: binary format)",
"   freq=0.5 ................. modelling frequency in Hz",
"   verbose=0 ................ silent option; >0 display info",
"   component=11 ............. receiver and source orientation",
"   nd=1000 .................. number of integration domains",
"   startlogx=-6 ............. first integration point in space",
"   deltalogx=0.025 .......... logarithmic sampling rate of",
"                              integration points in space",
"                              at first iteration",
"   nlogx=512 ................ amount of integration points in",
"                              space at first iteration",
"   kmax=0.628625 ............ largest wavenumber to be integrated",
"   c1=0.0 ................... first precision parameter",
"   c2=0.01 .................. second precision parameter",
"   maxpt=500 ................ maximum amount of integration points in space",
"   dopchip=0 ................ pchip interpolation (1) or linear interpolation (0)",
"   xdirect=1 ................ direct field in space domain (1)",
"                              or in the wavenumber domain (0)",
"",
" Code:   Juerg Hunziker ",
"         Jan Thorbecke ",
" Theory: Evert Slob", 
"         Juerg Hunziker",
" ",
" Last update: 11th of October 2013",
" ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	FILE	*fp_out;
	int	iz, nrz, nrc1, nrc2, nre1, nre2, nrm1, nrm2, fullspace, xdirect;
	int     nx, ny, ix, iy, verbose, above, *component, nd, nlogx;
	int     zrcv_layer, zsrc_layer, maxpt, dopchip, writebin;
	REAL   	kmax, dx, dy, fxpos, fypos, freq, *z;
	REAL   	*econdV, *econdH, *epermV, *epermH, *mpermV, *mpermH;
	REAL   	sumecondV, sumecondH, sumepermV, sumepermH, summpermV, summpermH;
	REAL   	avecondV, avecondH, avepermV, avepermH, avmpermV, avmpermH;
	REAL   	sumallpar;
	double 	t0, t1, t2;
	complex *data;
	REAL	zrcv, zsrc, startlogx, deltalogx, c1, c2;
	char    *file_out;
	size_t  nwrite;

/* ======================= Reading and determining parameters ==================== */

	/* For a description of the parameters to be read, see documentation 
	   at beginning of this file. */

	initargs(argc, argv);
	requestdoc(1);

	/* Time when program was started. */
	t0=wallclock_time();

	/* Print information to screen (verbose=1) or not (verbose=0)? */
	if(!getparint("verbose", &verbose)) verbose = 0;

	/* Read parameters zsrc, zrcv, dx, dy, nx and ny. */
	if(!getparfloat("zsrc", &zsrc)) verr("zsrc (source depth) must be specified.");
	if(!getparfloat("zrcv", &zrcv)) verr("zrcv (receiver depth) must be specified."); 
	if(!getparfloat("dx", &dx)) verr("dx (x-sampling) must be specified.");
	if(!getparfloat("dy", &dy)) verr("dy (y-sampling) must be specified.");
	if(!getparint("nx", &nx)) verr("nx (number of x-samples) must be specified.");
	if(!getparint("ny", &ny)) verr("ny (number of y-samples) must be specified.");
	/* If an uneven number is entered for nx or ny, it is made even. */
	if (fabs(nx/2.0-ceil(nx/2.0))>=0.000001) {nx = nx+1;}
	if (fabs(ny/2.0-ceil(ny/2.0))>=0.000001) {ny = ny+1;}
	if (verbose) vmess("source depth in meters = %f", zsrc);
	if (verbose) vmess("Receiver depth in meters = %f", zrcv);
	if (verbose) vmess("Sampling in x-direction in meters = %f", dx);
	if (verbose) vmess("Sampling in y-direction in meters = %f", dy);
	if (verbose) vmess("Number of samples in x-direction = %d", nx);
	if (verbose) vmess("Number of samples in y-direction = %d", ny);

	/* Read path and filename for output file. */
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("Parameter file_out not found, writing to file EMmod_output.bin");
		file_out = "EMmod_output.bin";
	}

	/* Write data in binary or ASCII format? */
	if(!getparint("writebin", &writebin)) writebin = 1;
	if (writebin==1) {
		if (verbose) vmess("Writing data in binary format (writebin=1).");
	} else {
		if (verbose) vmess("Writing data in ASCII format (writebin=0).");
	}

	/* Check if all the material parameters have the same amount of values given. */
	nrz  = countparval("z");
	if(nrz==0) verr("Depth of interfaces unspecified.");
	nrc1 = countparval("econdV");
	if(nrc1==0) verr("Vertical conductivities unspecified.");
	nrc2 = countparval("econdH");
	if(nrc2==0) verr("Horizontal conductivities unspecified.");
	nre1 = countparval("epermV");
	if(nre1==0) verr("Vertical electric permittivity unspecified.");
	nre2 = countparval("epermH");
	if(nre2==0) verr("Horizontal electric permittivity unspecified.");
	nrm1 = countparval("mpermV");
	if(nrm1==0) verr("Vertical magnetic permeability unspecified.");
	nrm2 = countparval("mpermH");
	if(nrm2==0) verr("Horizontal magnetic permeability unspecified.");
	assert(nrz==nrc1);
	assert(nrz==nrc2);
	assert(nrz==nre1);
	assert(nrz==nre2);
	assert(nrz==nrm1);
	assert(nrz==nrm2);

	/* Read the material parameters. */	
	z  = (REAL *)calloc(nrz,sizeof(REAL));
	econdV = (REAL *)calloc(nrz,sizeof(REAL));
	econdH = (REAL *)calloc(nrz,sizeof(REAL));
	epermV = (REAL *)calloc(nrz,sizeof(REAL));
	epermH = (REAL *)calloc(nrz,sizeof(REAL));
	mpermV = (REAL *)calloc(nrz,sizeof(REAL));
	mpermH = (REAL *)calloc(nrz,sizeof(REAL));
	getparfloat("z", z);
	getparfloat("econdV", econdV);
	getparfloat("econdH", econdH);
	getparfloat("epermV", epermV);
	getparfloat("epermH", epermH);
	getparfloat("mpermV", mpermV);
	getparfloat("mpermH", mpermH);

	/* Check if medium is a homogeneous fullspace. If that is the case, 
	   the EM-field is computed directly in the space-domain.
           Note: Also a stack of layers with the same material parameters is 
                 treated as a homogeneous fullspace. */
	sumecondV = 0.0;
	sumecondH = 0.0;
	sumepermV = 0.0;
	sumepermH = 0.0;
	summpermV = 0.0;
	summpermH = 0.0;
        for (iz=0; iz<nrz; iz++) {
		sumecondV = sumecondV + econdV[iz];
		sumecondH = sumecondH + econdH[iz];
		sumepermV = sumepermV + epermV[iz];
		sumepermH = sumepermH + epermH[iz];
		summpermV = summpermV + mpermV[iz];
		summpermH = summpermH + mpermH[iz];
	}
	avecondV = sumecondV/nrz;
	avecondH = sumecondH/nrz;
	avepermV = sumepermV/nrz;
	avepermH = sumepermH/nrz;
	avmpermV = summpermV/nrz;
	avmpermH = summpermH/nrz;
	sumecondV = 0.0;
	sumecondH = 0.0;
	sumepermV = 0.0;
	sumepermH = 0.0;
	summpermV = 0.0;
	summpermH = 0.0;
	for (iz=0; iz<nrz; iz++) {
		sumecondV = sumecondV+(econdV[iz]-avecondV)*(econdV[iz]-avecondV);
		sumecondH = sumecondH+(econdH[iz]-avecondH)*(econdH[iz]-avecondH);
		sumepermV = sumepermV+(epermV[iz]-avepermV)*(epermV[iz]-avepermV);
		sumepermH = sumepermH+(epermH[iz]-avepermH)*(epermH[iz]-avepermH);
		summpermV = summpermV+(mpermV[iz]-avmpermV)*(mpermV[iz]-avmpermV);
		summpermH = summpermH+(mpermH[iz]-avmpermH)*(mpermH[iz]-avmpermH);
	}
	sumallpar = sumecondV + sumecondH + sumepermV + sumepermH + summpermV + summpermH;
	if (sumallpar <= 10e-12) {
		if (verbose) vmess("Configuration is a homogeneous fullspace.");
		fullspace = 1;
	} else {
		if (verbose) vmess("Configuration is a stack of layers.");
		fullspace = 0;
	}

	/* The variable component is initialized as array of size 1, because then it 
	   is handled with a pointer, which ensures that emmod.c retrieves any changes
	   done to this variable in kxwmod.c. */
	component = (int *)calloc(1,sizeof(int));
	if(!getparint("component", component)) component[0]=11;

	if(!getparfloat("freq", &freq)) freq = 0.5;
	if (verbose) vmess("Frequency in Hertz = %f", freq);

	/* Parameters related to the integration. */
	if(!getparint("nd", &nd)) nd = 1000;
	if (verbose) vmess("Number of integration domains (nd) = %d", nd);
	if(!getparfloat("kmax", &kmax)) kmax = 0.628625;
	if (verbose) vmess("Largest wavenumber to be integrated in 1/m (kmax) = %f", kmax);
	if(!getparfloat("startlogx", &startlogx)) startlogx = -6;
	if (verbose) vmess("First integration point in space in meters (startlogx) = 10^%f", startlogx);
	if(!getparfloat("deltalogx", &deltalogx)) deltalogx = 0.025;
	if (verbose) vmess("Logarithmic sampling rate of integration points in space at first iteration (deltalogx) = %f", deltalogx);
	if(!getparint("nlogx", &nlogx)) nlogx = 512;
	if (verbose) vmess("Amount of integration points in space at first iteration (nlogx) = %d", nlogx);
	if(!getparfloat("c1", &c1)) c1 = 0.0;
	if (verbose) vmess("First precision parameter (c1) = %f", c1);
	if(!getparfloat("c2", &c2)) c2 = 0.01;
	if (verbose) vmess("Second precision parameter (c2) = %f", c2);
	if(!getparint("maxpt", &maxpt)) maxpt = 500;
	if (verbose) vmess("Maximum amount of integration points in space (maxpt) = %d", maxpt);

	/* Pchip interpolation or linear interpolation? */
	if(!getparint("dopchip", &dopchip)) dopchip = 0;
	if (dopchip==1) {
		if (verbose) vmess("Using pchip interpolation (dopchip=1).");
	} else {
		if (verbose) vmess("Using linear interpolation (dopchip=0).");
	}

	/* Determine where the source level is in relation to the receiver level.
           Note: If zsrc or zrcv are on a layer interface, the layer above the interface
                 is chosen. */
	zrcv_layer = -1;
	zsrc_layer = -1;
	if (zrcv<=z[1]) zrcv_layer = 0;
	if (zsrc<=z[1]) zsrc_layer = 0;
	for (iz=2; iz<nrz; iz++) {
		if (zrcv<=z[iz] && zrcv>z[iz-1]) zrcv_layer = iz-1;
		if (zsrc<=z[iz] && zsrc>z[iz-1]) zsrc_layer = iz-1;
	}
	/* If zrcv_layer or zsrc_layer remain -1 after the loop, 
	   the receivers or the source must be below all interfaces. */
	if (zrcv_layer == -1) zrcv_layer = nrz-1;
	if (zsrc_layer == -1) zsrc_layer = nrz-1;

	/* Determine if the receivers are above the source (above=1), 
           below the source (above=-1) or in the same layer as the source (above=0). */
	if (zrcv_layer == zsrc_layer) {
		above=0;
		if (verbose) vmess("Receivers are in the same layer as the source (above=0).");
	} else if (zrcv_layer < zsrc_layer) {
		above=1;
		if (verbose) vmess("Receivers are in a layer above the source (above=1).");
	} else {
		above=-1;
		if (verbose) vmess("Receivers are in a layer below the source (above=-1).");
	}

	/* Is the direct field computed in the space domain (xdirect=1) or in the 
           wavenumber domain (xdirect=0)? */
	if(!getparint("xdirect", &xdirect)) xdirect = 1;
	if (above==0 && fullspace==0) {
		if (xdirect==1) {
			if (verbose) vmess("The direct field is computed in the space domain (xdirect=1).");
		} else {
			if (verbose) vmess("The direct field is computed in the wavenumber domain (xdirect=0).");
		}
	}

/* ============ start modeling in wavenumber domain =============== */

	/* Allocate matrix for output data. */
	data = (complex *)calloc(nx*ny,sizeof(complex));

	/* Time when the actual modeling starts. */	
	t1=wallclock_time();

	/* Call kxwmod with all the loaded parameters to do the actual modeling. */
	kxwmod(data, freq, nx, dx, ny, dy, nrz, z, econdV, econdH, epermV, epermH, mpermV, mpermH, zsrc, zrcv, zrcv_layer, zsrc_layer, above, component, nd, kmax, startlogx, deltalogx, nlogx, c1, c2, maxpt, dopchip, fullspace, xdirect, verbose);

	/* Write modeling result to disk if the component that was entered is a component that
           was actually computed. */
	if (component[0]==11 || component[0]==12 || component[0]==13 || component[0]==21 || component[0]==22 || component[0]==23 || component[0]==31 || component[0]==32 || component[0]==33 || component[0]==41 || component[0]==42 || component[0]==43 || component[0]==51 || component[0]==52 || component[0]==53 || component[0]==61 || component[0]==62 || component[0]==63 || component[0]==77 || component[0]==88) {
		if (writebin==1) { // write data in binary format
			if (file_out==NULL) fp_out=stdout;
			else fp_out = fopen(file_out,"w");
			if (fp_out == NULL) verr("error in creating output file");
			nwrite = fwrite( &data[0].r, sizeof(complex), nx*ny, fp_out);
			assert(nwrite == nx*ny);
    			fclose(fp_out);
		} else { // write data in ASCII format
			if (file_out==NULL) fp_out=stdout;
			else fp_out = fopen(file_out,"w");
			if (fp_out == NULL) verr("error in creating output file");
			fprintf(fp_out,"x y z real_data imag_data\n");
			fypos = -ny/2*dy;
			for (iy=0; iy<ny; iy++) {
				fxpos = -nx/2*dx;
				for (ix=0; ix<nx; ix++) {
					fprintf(fp_out,"%f %f %f %e %e\n",fxpos,fypos,zrcv,data[nx*iy+ix].r,data[nx*iy+ix].i);
					fxpos = fxpos+dx;
				}
				fypos = fypos+dy;
			}
    			fclose(fp_out);
		}
	}

/* ============ close files and free memory =============== */

	/* Time when modeling is completed. */
	t2=wallclock_time();

	/* Show runtimes. */
	if (verbose) {
		vmess("emmod runtimes: ");
		vmess(" - initial time=%.3f", t1-t0);
		vmess(" - kxwmod  time=%.3f", t2-t1);
		vmess(" - total   time=%.3f", t2-t0);
	}
		
	/* Deallocate matrix data. */
	if( data ) free(data);

	exit ( 0 );
}

