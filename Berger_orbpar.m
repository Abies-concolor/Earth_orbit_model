%This function was translated from the original FORTRAN code
%provided at the NASA GISS website http://aom.giss.nasa.gov/
%specific link no longer available as of September 2013.
%The original FORTRAN code is due to Gary L. Russell and is provided here
%as comments for reference. 

%The algorithm itself is provided by Berger (1978) - see references below. 


function [ecc,obliq,omega_bar] = Berger_orbpar(yr)

%       SUBROUTINE ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
% C****
% C**** ORBPAR calculates the three orbital parameters as a function of
% C**** YEAR.  The source of these calculations is: Andre L. Berger,
% C**** 1978, "Long-Term Variations of Daily Insolation and Quaternary
% C**** Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
% C**** Berger, May 1978, "A Simple Algorithm to Compute Long Term
% C**** Variations of Daily Insolation", published by Institut
% C**** D'Astronomie de Geophysique, Universite Catholique de Louvain,
% C**** Louvain-la Neuve, No. 18.
% C****
% C**** Tables and equations refer to the first reference (JAS).  The
% C**** corresponding table or equation in the second reference is
% C**** enclosed in parentheses.  The coefficients used in this
% C**** subroutine are slightly more precise than those used in either
% C**** of the references.  The generated orbital parameters are precise
% C**** within plus or minus 1000000 years from present.
% C****
% C**** Input:  YEAR   = years A.D. are positive, B.C. are negative
% C**** Output: ECCEN  = eccentricity of orbital ellipse
% C****         OBLIQ  = latitude of Tropic of Cancer in radians
% C****         OMEGVP = longitude of perihelion =
% C****                = spatial angle from vernal equinox to perihelion
% C****                  in radians with sun as angle vertex
% C****
%       IMPLICIT REAL*8 (A-H,O-Z)
%       PARAMETER (TWOPI=6.283185307179586477, PIz180=TWOPI/360.)
twopi = 2*pi;
piz180 = twopi/360;

%REAL*8 TABLE1(3,47),TABLE4(3,19),TABLE5(3,78)
%C**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
t1 = [-2462.2214466, 31.609974, 251.9025;
    -857.3232075, 32.620504, 280.8325;
    -629.3231835, 24.172203, 128.3057;
    -414.2804924, 31.983787, 292.7252;
    -311.7632587, 44.828336,  15.3747;
    308.9408604, 30.973257, 263.7951;
    -162.5533601, 43.668246, 308.4258;
    -116.1077911, 32.246691, 240.0099;
    101.1189923, 30.599444, 222.9725;
    -67.6856209, 42.681324, 268.7809;
    24.9079067, 43.836462, 316.7998;
    22.5811241, 47.439436, 319.6024;
    -21.1648355, 63.219948, 143.8050;
    -15.6549876, 64.230478, 172.7351;
    15.3936813,  1.010530,  28.9300;
    14.6660938,  7.437771, 123.5968;
    -11.7273029, 55.782177,  20.2082;
    10.2742696,   .373813,  40.8226;
    6.4914588, 13.218362, 123.4722;
    5.8539148, 62.583231, 155.6977;
    -5.4872205, 63.593761, 184.6277;
    -5.4290191, 76.438310, 267.2772;
    5.1609570, 45.815258,  55.0196;
    5.0786314,  8.448301, 152.5268;
    -4.0735782, 56.792707,  49.1382;
    3.7227167, 49.747842, 204.6609;
    3.3971932, 12.058272,  56.5233;
    -2.8347004, 75.278220, 200.3284;
    -2.6550721, 65.241008, 201.6651;
    -2.5717867, 64.604291, 213.5577;
    -2.4712188,  1.647247,  17.0374;
    2.4625410,  7.811584, 164.4194;
    2.2464112, 12.207832,  94.5422;
    -2.0755511, 63.856665, 131.9124;
    -1.9713669, 56.155990,  61.0309;
    -1.8813061, 77.448840, 296.2073;
    -1.8468785,  6.801054, 135.4894;
    1.8186742, 62.209418, 114.8750;
    1.7601888, 20.656133, 247.0691;
    -1.5428851, 48.344406, 256.6114;
    1.4738838, 55.145460,  32.1008;
    -1.4593669, 69.000539, 143.6804;
    1.4192259, 11.071350,  16.8784;
    -1.1818980, 74.291298, 160.6835;
    1.1756474, 11.047742,  27.5932;
    -1.1316126,  0.636717, 348.1074;
    1.0896928, 12.844549,  82.6496];


% C**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
t4 = [ .01860798,  4.207205,  28.620089;
    .01627522,  7.346091, 193.788772;
    -.01300660, 17.857263, 308.307024;
    .00988829, 17.220546, 320.199637;
    -.00336700, 16.846733, 279.376984;
    .00333077,  5.199079,  87.195000;
    -.00235400, 18.231076, 349.129677;
    .00140015, 26.216758, 128.443387;
    .00100700,  6.359169, 154.143880;
    .00085700, 16.210016, 291.269597;
    .00064990,  3.065181, 114.860583;
    .00059900, 16.583829, 332.092251;
    .00037800, 18.493980, 296.414411;
    -.00033700,  6.190953, 145.769910;
    .00027600, 18.867793, 337.237063;
    .00018200, 17.425567, 152.092288;
    -.00017400,  6.186001, 126.839891;
    -.00012400, 18.417441, 210.667199;
    .00001250,  0.667863,  72.108838];

% C**** Table 5 (3).  General precession in longitude: psi
t5 = [ 7391.0225890, 31.609974, 251.9025;
    2555.1526947, 32.620504, 280.8325;
    2022.7629188, 24.172203, 128.3057;
    -1973.6517951,  0.636717, 348.1074;
    1240.2321818, 31.983787, 292.7252;
    953.8679112,  3.138886, 165.1686;
    -931.7537108, 30.973257, 263.7951;
    872.3795383, 44.828336,  15.3747;
    606.3544732,  0.991874,  58.5749;
    -496.0274038,  0.373813,  40.8226;
    456.9608039, 43.668246, 308.4258;
    346.9462320, 32.246691, 240.0099;
    -305.8412902, 30.599444, 222.9725;
    249.6173246,  2.147012, 106.5937;
    -199.1027200, 10.511172, 114.5182;
    191.0560889, 42.681324, 268.7809;
    -175.2936572, 13.650058, 279.6869;
    165.9068833,  0.986922,  39.6448;
    161.1285917,  9.874455, 126.4108;
    139.7878093, 13.013341, 291.5795;
    -133.5228399,  0.262904, 307.2848;
    117.0673811,  0.004952,  18.9300;
    104.6907281,  1.142024, 273.7596;
    95.3227476, 63.219948, 143.8050;
    86.7824524,  0.205021, 191.8927;
    86.0857729,  2.151964, 125.5237;
    70.5893698, 64.230478, 172.7351;
    -69.9719343, 43.836462, 316.7998;
    -62.5817473, 47.439436, 319.6024;
    61.5450059,  1.384343,  69.7526;
    -57.9364011,  7.437771, 123.5968;
    57.1899832, 18.829299, 217.6432;
    -57.0236109,  9.500642,  85.5882;
    -54.2119253,  0.431696, 156.2147;
    53.2834147,  1.160090,  66.9489;
    52.1223575, 55.782177,  20.2082;
    -49.0059908, 12.639528, 250.7568;
    -48.3118757,  1.155138,  48.0188;
    -45.4191685,  0.168216,   8.3739;
    -42.2357920,  1.647247,  17.0374;
    -34.7971099, 10.884985, 155.3409;
    34.4623613,  5.610937,  94.1709;
    -33.8356643, 12.658184, 221.1120;
    33.6689362,  1.010530,  28.9300;
    -31.2521586,  1.983748, 117.1498;
    -30.8798701, 14.023871, 320.5095;
    28.4640769,  0.560178, 262.3602;
    -27.1960802,  1.273434, 336.2148;
    27.0860736, 12.021467, 233.0046;
    -26.3437456, 62.583231, 155.6977;
    24.7253740, 63.593761, 184.6277;
    24.6732126, 76.438310, 267.2772;
    24.4272733,  4.280910,  78.9281;
    24.0127327, 13.218362, 123.4722;
    21.7150294, 17.818769, 188.7132;
    -21.5375347,  8.359495, 180.1364;
    18.1148363, 56.792707,  49.1382;
    -16.9603104,  8.448301, 152.5268;
    -16.1765215,  1.978796,  98.2198;
    15.5567653,  8.863925,  97.4808;
    15.4846529,  0.186365, 221.5376;
    15.2150632,  8.996212, 168.2438;
    14.5047426,  6.771027, 161.1199;
    -14.3873316, 45.815258,  55.0196;
    13.1351419, 12.002811, 262.6495;
    12.8776311, 75.278220, 200.3284;
    11.9867234, 65.241008, 201.6651;
    11.9385578, 18.870667, 294.6547;
    11.7030822, 22.009553,  99.8233;
    11.6018181, 64.604291, 213.5577;
    -11.2617293, 11.498094, 154.1631;
    -10.4664199,  0.578834, 232.7153;
    10.4333970,  9.237738, 138.3034;
    -10.2377466, 49.747842, 204.6609;
    10.1934446,  2.147012, 106.5938;
    -10.1280191,  1.196895, 250.4676;
    10.0289441,  2.133898, 332.3345;
    -10.0034259,  0.173168,  27.3039];

ym1950 = yr-1950;

% C**** Obliquity from Table 1 (2):
% C****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
% C****   OBLIQD = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)

SUMC = 0;
for i = 1:47
    ARG    = piz180*(ym1950*t1(i,2)/3600+t1(i,3));
    SUMC   = SUMC + t1(i,1)*cos(ARG);
end
obliq = 23.320556 + SUMC/3600;
obliq  = obliq*piz180;


% C**** Eccentricity from Table 4 (1):
% C****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
% C****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
% C****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
% C****
ESINPI = 0;
ECOSPI = 0;

for i = 1:19
    ARG    = piz180*(ym1950*t4(i,2)/3600+t4(i,3));
    ESINPI = ESINPI + t4(i,1)*sin(ARG);
    ECOSPI = ECOSPI + t4(i,1)*cos(ARG);
end
ECCEN  = sqrt(ESINPI*ESINPI+ECOSPI*ECOSPI);
ecc = ECCEN; %First function output


% C**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
% C****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
% C****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
% C****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
% C****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
% C****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
% C****
PIE    = atan2(ESINPI,ECOSPI);
FSINFD = 0;
for i = 1:78
    ARG    = piz180*(ym1950*t5(i,2)/3600+t5(i,3));
    FSINFD = FSINFD + t5(i,1)*sin(ARG);
end
PSI    = piz180*(3.392506+(ym1950*50.439273+FSINFD)/3600);
OMEGVP = mod(PIE+PSI+.5*twopi, twopi);
omega_bar = OMEGVP;

% C**** ORBPAR.FOR    Calculate ORBital PARameters    2002/09/11
% C****
% C**** Output: ECCEN = eccentricity of the Earth's orbit about the Sun
% C****         OBLIQ = obliquity = latitude of Tropic of Cancer (degrees)
% C****        OMEGVP = longitude of perihelion =
% C****               = spatial angle from vernal equinox to perihelion
% C****                 with Sun as angle vertex
% C****
%       IMPLICIT REAL*8 (A-H,M,O-Z)
%       PARAMETER (TWOPI=6.283185307179586477, PIz180=TWOPI/360.)
%       CHARACTER ARG*80
% C****
%       NARGS = IARGC()
%       IF(NARGS.le.0)  GO TO 800
% C****
% C**** Decode command line arguments
% C****
%       CALL GETARG (1,ARG)
%       READ (ARG,*,ERR=810) IYMIN
%       IF(ABS(IYMIN).gt.5000000) GO TO 811
%       IYMAX = IYMIN
%       IYINC = 1000
%       IF(NARGS.le.1)  GO TO 150
%       CALL GETARG (2,ARG)
%       READ (ARG,*,ERR=812) IYMAX
%       IF(NARGS.le.2)  GO TO 150
%       CALL GETARG (3,ARG)
%       READ (ARG,*,ERR=813) IYINC
%       IF(IYINC.eq.0)  IYINC = 1000
% C**** Limit calculations to 1001 years
%   150 IF((IYINC.gt.0 .and. IYMAX.gt.IYMIN+1000*IYINC) .or.
%      *   (IYINC.lt.0 .and. IYMAX.lt.IYMIN+1000*IYINC))
%      *                     IYMAX =  IYMIN+1000*IYINC
%       IF(IYMAX.gt. 5000000)  IYMAX =  5000000
%       IF(IYMAX.lt.-5000000)  IYMAX = -5000000
% C****
% C**** Loop over years
% C****
%       WRITE (6,920)
%       DO 210 IYEAR = IYMIN,IYMAX,IYINC
% C**** Determine orbital parameters
%       YEAR = IYEAR
%       CALL ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
%       OBLIQ  = OBLIQ /PIz180
%       OMEGVP = OMEGVP/PIz180
%   210 WRITE (6,921) IYEAR, ECCEN,OBLIQ,OMEGVP
%       GO TO 999
% C****
%   800 WRITE (6,*)
%      *'Usage: ORBPAR Ymin(A.D.) [Ymax Yinc]                2002/09/11'
%       WRITE (6,*)
%      *'       Calculate Earth''s orbital parameters for range of years'
%       WRITE (6,*)
%      *'       Parameters are precise + or - 1000000 years from present'
%       GO TO 999
%   810 WRITE (0,*) 'Unable to decipher minimum year from first' //
%      *            ' command line argument: ',TRIM(ARG)
%       STOP 810
%   811 WRITE (0,*) 'Years are limited to be between 5000000 BC and' //
%      *            ' 5000000 AD.'
%       STOP 811
%   812 WRITE (0,*) 'Unable to decipher maximum year from second' //
%      *            ' command line argument: ',TRIM(ARG)
%       STOP 812
%   813 WRITE (0,*) 'Unable to decipher yearly increment from third' //
%      *            ' command line argument: ',TRIM(ARG)
%       STOP 813
% C****
%   900 FORMAT ('Content-type: text/plain' /)
%   920 FORMAT ('Orbital Parmameters'
%      *     // '                                     Long. of',
%      *      / '   Year     Eccentri    Obliquity    Perihel.',
%      *      / '  (A.D.)      city      (degrees)    (degrees)',
%      *      / '  ------    --------    ---------    --------')
%   921 FORMAT (I8,F11.6,F13.4,F13.3)
%   999 END