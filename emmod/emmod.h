#include "par.h"
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#if defined(DOUBLE)
#define REAL double
#define getparfloat getpardouble
#else
#define REAL float
#endif

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)


#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    REAL r,i;
} complex;
#endif/* complex */

/* C Functions */

double wallclock_time(void);

void kxwmod(complex *data, REAL freq, int nx, REAL dx, int ny, REAL dy, int nz, REAL *z, REAL *econdV, REAL *econdH, REAL *epermV, REAL *epermH, REAL *mpermV, REAL *mpermH, REAL zsrc, REAL zrcv, int zrcv_layer, int zsrc_layer, int above, int *component, int nd, REAL kmax, REAL startlogx, REAL deltalogx, int nlogx, REAL c1, REAL c2, int maxpt, int dopchip, int fullspace, int xdirect, int verbose);

/* Fortran Functions*/

void gammarad_(complex *gam, int *nkx, REAL *dkx, complex *etaV, complex *etaH, complex *gamA);

void wprop (complex *Wup, complex *Wdown, complex *Gamma, int *nk, REAL cor, int *nz, REAL *z, REAL *zrcv, int *zrcv_layer);

void rplus_(complex *Rp, int *corel, REAL *cor, complex *etaH, complex *Gamma, int *nz, REAL *z, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void rmin_(complex *Rm, int *corel, REAL *cor, complex *etaH, complex *Gamma, int *nz, REAL *z, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void pupmin_(complex *Pupmin, complex *Rp, complex *Rm, complex *Gamma, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void pdownmin_(complex *Pdownmin, complex *Rp, complex *Rm, complex *Gamma, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void pupplus_(complex *Pupplus, complex *Rp, complex *Rm, complex *Gamma, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void pdownplus_(complex *Pdownplus, complex *Rp, complex *Rm, complex *Gamma, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, int *zrcv_layer, int *zsrc_layer, int *nlayers, int *above);

void ptotalxx_(complex *temptot1, complex *temptot2, complex *Pdownmin, complex *Pupmin, complex *Pdownplusbar, complex *Pupplusbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalxy_(complex *temptot, complex *Pdownmin, complex *Pupmin, complex *Pdownplusbar, complex *Pupplusbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rbpar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalxz_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, int *nlayers, int *above, int *xdirect);

void ptotalyy_(complex *temptot1, complex *temptot2, complex *Pdownmin, complex *Pupmin, complex *Pdownplusbar, complex *Pupplusbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalyz_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, int *nlayers, int *above, int *xdirect);

void ptotalzx_(complex *temptot, complex *Pdownmin, complex *Pupmin, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, int *nlayers, int *above, int *xdirect);

void ptotalzy_(complex *temptot, complex *Pdownmin, complex *Pupmin, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, int *nlayers, int *above, int *xdirect);

void ptotalzz_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamA, int *nlayers, int *above, int *xdirect);

void ptotalxxm_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Pdownminbar, complex *Pupminbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalxym_(complex *temptot1, complex *temptot2, complex *Pdownplus, complex *Pupplus, complex *Pdownminbar, complex *Pupminbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalxzm_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalyxm_(complex *temptot1, complex *temptot2, complex *Pdownplus, complex *Pupplus, complex *Pdownminbar, complex *Pupminbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalyym_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Pdownminbar, complex *Pupminbar, complex *Wup, complex *Wdown, complex *Wupbar, complex *Wdownbar, complex *Gamma, complex *GammaB, complex *Rp, complex *Rpbar, complex *Rm, complex *Rmbar, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalyzm_(complex *temptot, complex *Pdownplus, complex *Pupplus, complex *Wup, complex *Wdown, complex *Gamma, complex *Rp, complex *Rm, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zsrc, REAL *zrcv, int *zrcv_layer, int *zsrc_layer, complex *gamA, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalzxm_(complex *temptot, complex *Pdownplusbar, complex *Pupplusbar, complex *Wupbar, complex *Wdownbar, complex *GammaB, complex *Rpbar, complex *Rmbar, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalzym_(complex *temptot, complex *Pdownplusbar, complex *Pupplusbar, complex *Wupbar, complex *Wdownbar, complex *GammaB, complex *Rpbar, complex *Rmbar, complex *zetaH, complex *zetaV, int *corel, REAL *cor, int *nz, REAL *z, REAL *zrcv, REAL *zsrc, int *zsrc_layer, int *zrcv_layer, complex *gamB, int *nlayers, int *above, int *xdirect);

void ptotalzzm_(complex *Ptot, int *nxh, int *nyh);

void ptotalref_(complex *temptot, complex *Wup, complex *Rp, int *corel);

void getcoords_(REAL *xpos, REAL *startlogx, REAL *deltalogx, int *nlogx);

void hankeltransxx_(complex *xtotrad1, complex *xtotrad2, int *marker, complex *temptot1, complex *temptot2, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransxy_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransxz_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransyy_(complex *xtotrad1, complex *xtotrad2, int *marker, complex *temptot1, complex *temptot2, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransyz_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltranszx_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltranszy_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltranszz_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransxxm_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransxym_(complex *xtotrad1, complex *xtotrad2, int *marker, complex *temptot1, complex *temptot2, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransxzm_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransyxm_(complex *xtotrad1, complex *xtotrad2, int *marker, complex *temptot1, complex *temptot2, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransyym_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransyzm_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltranszxm_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltranszym_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void hankeltransref_(complex *xtotrad1, int *marker, complex *temptot, int *corel, REAL *cor, REAL *xpos, int *nlogx, int *nd, REAL *kmax);

void evalpoints_(REAL *xposnew, complex *xtotrad1new, complex *xtotrad2new, int *markernew, int *nlogxnew, int *nlogxdata, int *newit, REAL *xpos, complex *xtotrad1, complex *xtotrad2, int *nlogx, REAL *c1, REAL *c2);

void evalpoints_mono_(REAL *xposnew, complex *xtotrad1new, int *markernew, int *nlogxnew, int *nlogxdata, int *newit, REAL *xpos, complex *xtotrad1, int *nlogx, REAL *c1, REAL *c2);

void evalpoints_lin_(REAL *xposnew, complex *xtotrad1new, complex *xtotrad2new, int *markernew, int *nlogxnew, int *nlogxdata, int *newit, REAL *xpos, complex *xtotrad1, complex *xtotrad2, int *nlogx, REAL *c1, REAL *c2);

void evalpoints_lin_mono_(REAL *xposnew, complex *xtotrad1new, int *markernew, int *nlogxnew, int *nlogxdata, int *newit, REAL *xpos, complex *xtotrad1, int *nlogx, REAL *c1, REAL *c2);

void gridit_xx_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xx_lin_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xx_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_xy_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xy_ilin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xy_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_xz_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_xz_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_xz_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA);

void gridit_yy_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yy_lin_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yy_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_yz_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_yz_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_yz_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA);

void gridit_zx_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zx_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zx_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA);

void gridit_zy_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zy_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zy_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA);

void gridit_zz_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zz_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, int *above, int *xdirect);

void gridit_zz_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA);

void gridit_xxm_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xxm_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xxm_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_xym_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xym_lin_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_xym_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_xzm_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA, int *above, int *xdirect);

void gridit_xzm_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA, int *above, int *xdirect);

void gridit_xzm_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA);

void gridit_yxm_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yxm_lin_(complex *Ptot, complex *xtotrad1, complex *xtotrad2, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yxm_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_yym_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yym_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB, int *above, int *xdirect);

void gridit_yym_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *zetaV, complex *gamA, complex *gamB);

void gridit_yzm_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA, int *above, int *xdirect);

void gridit_yzm_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA, int *above, int *xdirect);

void gridit_yzm_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *etaH, complex *etaV, complex *zetaH, complex *gamA);

void gridit_zxm_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB, int *above, int *xdirect);

void gridit_zxm_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB, int *above, int *xdirect);

void gridit_zxm_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB);

void gridit_zym_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB, int *above, int *xdirect);

void gridit_zym_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB, int *above, int *xdirect);

void gridit_zym_fullspace_(complex *Ptot, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *zrcv, REAL *zsrc, complex *zetaH, complex *zetaV, complex *gamB);

void gridit_ref_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx);

void gridit_ref_lin_(complex *Ptot, complex *xtotrad1, REAL *dx, int *nxh, REAL *dy, int *nyh, REAL *xpos, int *nlogx);
