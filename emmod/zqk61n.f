*DECK ZQK61
      SUBROUTINE ZQK61N (F, A, B, RESULT)
C***BEGIN PROLOGUE  ZQK61
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      DOUBLE PRECISION (QK61-S, ZQK61-D)
C***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rule
C        Standard fortran subroutine
C        Double precision version
C
C
C        PARAMETERS
C         ON ENTRY
C           F      - Double precision complex
C                    Function subprogram defining the integrand
C                    function F(X). The actual name for F needs to be
C                    declared E X T E R N A L in the calling program.
C
C           A      - Double precision
C                    Lower limit of integration
C
C           B      - Double precision
C                    Upper limit of integration
C
C         ON RETURN
C           RESULT - Double precision complex
C                    Approximation to the integral I
C                    RESULT is computed by applying the 61-point
C                    Kronrod rule (RESK) obtained by optimal addition of
C                    abscissae to the 30-point Gauss rule (RESG).
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  NONE
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  ZQK61
C
      implicit none
      complex*16 f,fc,fsum,fval1,fval2,fv1,fv2,resg,resk,reskh,result
      DOUBLE PRECISION A,DABSC,ABSERR,B,CENTR,DHLGTH,
     1  EPMACH,HLGTH,RESABS,RESASC,UFLOW,WG,WGK,XGK
      INTEGER J,JTW,JTWM1
C
      DIMENSION FV1(30),FV2(30),XGK(31),WGK(31),WG(15), F(61)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
C           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
C           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
C                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
C                   GAUSS RULE
C                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
C                   TO THE 30-POINT GAUSS RULE
C
C           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
C
C           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
C
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      SAVE WG, XGK, WGK
      DATA WG  (  1) / 0.0079681924 9616660561 5465883474 674 D0 /
      DATA WG  (  2) / 0.0184664683 1109095914 2302131912 047 D0 /
      DATA WG  (  3) / 0.0287847078 8332336934 9719179611 292 D0 /
      DATA WG  (  4) / 0.0387991925 6962704959 6801936446 348 D0 /
      DATA WG  (  5) / 0.0484026728 3059405290 2938140422 808 D0 /
      DATA WG  (  6) / 0.0574931562 1761906648 1721689402 056 D0 /
      DATA WG  (  7) / 0.0659742298 8218049512 8128515115 962 D0 /
      DATA WG  (  8) / 0.0737559747 3770520626 8243850022 191 D0 /
      DATA WG  (  9) / 0.0807558952 2942021535 4694938460 530 D0 /
      DATA WG  ( 10) / 0.0868997872 0108297980 2387530715 126 D0 /
      DATA WG  ( 11) / 0.0921225222 3778612871 7632707087 619 D0 /
      DATA WG  ( 12) / 0.0963687371 7464425963 9468626351 810 D0 /
      DATA WG  ( 13) / 0.0995934205 8679526706 2780282103 569 D0 /
      DATA WG  ( 14) / 0.1017623897 4840550459 6428952168 554 D0 /
      DATA WG  ( 15) / 0.1028526528 9355884034 1285636705 415 D0 /
C
      DATA XGK (  1) / 0.9994844100 5049063757 1325895705 811 D0 /
      DATA XGK (  2) / 0.9968934840 7464954027 1630050918 695 D0 /
      DATA XGK (  3) / 0.9916309968 7040459485 8628366109 486 D0 /
      DATA XGK (  4) / 0.9836681232 7974720997 0032581605 663 D0 /
      DATA XGK (  5) / 0.9731163225 0112626837 4693868423 707 D0 /
      DATA XGK (  6) / 0.9600218649 6830751221 6871025581 798 D0 /
      DATA XGK (  7) / 0.9443744447 4855997941 5831324037 439 D0 /
      DATA XGK (  8) / 0.9262000474 2927432587 9324277080 474 D0 /
      DATA XGK (  9) / 0.9055733076 9990779854 6522558925 958 D0 /
      DATA XGK ( 10) / 0.8825605357 9205268154 3116462530 226 D0 /
      DATA XGK ( 11) / 0.8572052335 4606109895 8658510658 944 D0 /
      DATA XGK ( 12) / 0.8295657623 8276839744 2898119732 502 D0 /
      DATA XGK ( 13) / 0.7997278358 2183908301 3668942322 683 D0 /
      DATA XGK ( 14) / 0.7677774321 0482619491 7977340974 503 D0 /
      DATA XGK ( 15) / 0.7337900624 5322680472 6171131369 528 D0 /
      DATA XGK ( 16) / 0.6978504947 9331579693 2292388026 640 D0 /
      DATA XGK ( 17) / 0.6600610641 2662696137 0053668149 271 D0 /
      DATA XGK ( 18) / 0.6205261829 8924286114 0477556431 189 D0 /
      DATA XGK ( 19) / 0.5793452358 2636169175 6024932172 540 D0 /
      DATA XGK ( 20) / 0.5366241481 4201989926 4169793311 073 D0 /
      DATA XGK ( 21) / 0.4924804678 6177857499 3693061207 709 D0 /
      DATA XGK ( 22) / 0.4470337695 3808917678 0609900322 854 D0 /
      DATA XGK ( 23) / 0.4004012548 3039439253 5476211542 661 D0 /
      DATA XGK ( 24) / 0.3527047255 3087811347 1037207089 374 D0 /
      DATA XGK ( 25) / 0.3040732022 7362507737 2677107199 257 D0 /
      DATA XGK ( 26) / 0.2546369261 6788984643 9805129817 805 D0 /
      DATA XGK ( 27) / 0.2045251166 8230989143 8957671002 025 D0 /
      DATA XGK ( 28) / 0.1538699136 0858354696 3794672743 256 D0 /
      DATA XGK ( 29) / 0.1028069379 6673703014 7096751318 001 D0 /
      DATA XGK ( 30) / 0.0514718425 5531769583 3025213166 723 D0 /
      DATA XGK ( 31) / 0.0000000000 0000000000 0000000000 000 D0 /
C
      DATA WGK (  1) / 0.0013890136 9867700762 4551591226 760 D0 /
      DATA WGK (  2) / 0.0038904611 2709988405 1267201844 516 D0 /
      DATA WGK (  3) / 0.0066307039 1593129217 3319826369 750 D0 /
      DATA WGK (  4) / 0.0092732796 5951776342 8441146892 024 D0 /
      DATA WGK (  5) / 0.0118230152 5349634174 2232898853 251 D0 /
      DATA WGK (  6) / 0.0143697295 0704580481 2451432443 580 D0 /
      DATA WGK (  7) / 0.0169208891 8905327262 7572289420 322 D0 /
      DATA WGK (  8) / 0.0194141411 9394238117 3408951050 128 D0 /
      DATA WGK (  9) / 0.0218280358 2160919229 7167485738 339 D0 /
      DATA WGK ( 10) / 0.0241911620 7808060136 5686370725 232 D0 /
      DATA WGK ( 11) / 0.0265099548 8233310161 0601709335 075 D0 /
      DATA WGK ( 12) / 0.0287540487 6504129284 3978785354 334 D0 /
      DATA WGK ( 13) / 0.0309072575 6238776247 2884252943 092 D0 /
      DATA WGK ( 14) / 0.0329814470 5748372603 1814191016 854 D0 /
      DATA WGK ( 15) / 0.0349793380 2806002413 7499670731 468 D0 /
      DATA WGK ( 16) / 0.0368823646 5182122922 3911065617 136 D0 /
      DATA WGK ( 17) / 0.0386789456 2472759295 0348651532 281 D0 /
      DATA WGK ( 18) / 0.0403745389 5153595911 1995279752 468 D0 /
      DATA WGK ( 19) / 0.0419698102 1516424614 7147541285 970 D0 /
      DATA WGK ( 20) / 0.0434525397 0135606931 6831728117 073 D0 /
      DATA WGK ( 21) / 0.0448148001 3316266319 2355551616 723 D0 /
      DATA WGK ( 22) / 0.0460592382 7100698811 6271735559 374 D0 /
      DATA WGK ( 23) / 0.0471855465 6929915394 5261478181 099 D0 /
      DATA WGK ( 24) / 0.0481858617 5708712914 0779492298 305 D0 /
      DATA WGK ( 25) / 0.0490554345 5502977888 7528165367 238 D0 /
      DATA WGK ( 26) / 0.0497956834 2707420635 7811569379 942 D0 /
      DATA WGK ( 27) / 0.0504059214 0278234684 0893085653 585 D0 /
      DATA WGK ( 28) / 0.0508817958 9874960649 2297473049 805 D0 /
      DATA WGK ( 29) / 0.0512215478 4925877217 0656282604 944 D0 /
      DATA WGK ( 30) / 0.0514261285 3745902593 3862879215 781 D0 /
      DATA WGK ( 31) / 0.0514947294 2945156755 8340433647 099 D0 /
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           DABSC  - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 30-POINT GAUSS RULE
C           RESK   - RESULT OF THE 61-POINT KRONROD RULE
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
C                    OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  ZQK61
C      EPMACH = D1MACH(4)
      EPMACH = 2.0**(-52)
C      UFLOW = D1MACH(1)
      UFLOW = 2.0**(-1022)
C
C      CENTR = 0.5D+00*(B+A)
      HLGTH = 0.5D+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0D+00
C      FC = F(CENTR)
      FC = F(61)
      RESK = WGK(31)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,15
        JTW = J*2
        DABSC = HLGTH*XGK(JTW)
C        FVAL1 = F(CENTR-DABSC)
C        FVAL2 = F(CENTR+DABSC)
        FVAL1 = F(JTW*2-1)
        FVAL2 = F(JTW*2)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J=1,15
        JTWM1 = J*2-1
        DABSC = HLGTH*XGK(JTWM1)
C        FVAL1 = F(CENTR-DABSC)
C        FVAL2 = F(CENTR+DABSC)
        FVAL1 = F(JTWM1*2-1)
        FVAL2 = F(JTWM1*2)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  15    CONTINUE
      RESKH = RESK*0.5D+00
      RESASC = WGK(31)*ABS(FC-RESKH)
      DO 20 J=1,30
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     1  ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5D+02)*RESABS,ABSERR)
      RETURN
      END
