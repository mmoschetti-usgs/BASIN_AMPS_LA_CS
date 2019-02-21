%% nshmp-haz Ground Motion Model (GMM) calculator example script

% =========================================================================
% This script provides instruction on how to access ground motion models
% (GMMs) implemented in the nshmp-haz library.
%
% See the README for instructions on how to set up Matlab to use nshmp-haz.
% =========================================================================

clear
clc

% Specify path to nshmp-haz library:
nshmpHaz = '/Users/ashumway/USGS_Computer_Codes/nshmp-haz/nshmp-haz/dist/nshmp-haz.jar';

% Make Matlab aware of nshmp-haz by adding it to the 'dynamic' classpath:
javaaddpath(nshmpHaz);

% Alternatively, one can add nshmp-haz to the faster 'static' classpath by
% saving a file with the name `javaclasspath.txt` to the Matlab preferences
% directory, as specified by the `prefdir` command, and with single line:
%
%     /path/to/repository/nshmp-haz/dist/nshmp-haz.jar
%
% Although the static classpath is generally a little faster, you must
% restart Matlab any time nshmp-haz.jar is rebuilt with updated code.

import org.opensha2.etc.*

% =========================================================================
% Single model ground motion calculation:

% Initialize calculator:
hazMat = HazMat.init(nshmpHaz);

% Note that hazMat is stateless and reusable and should therefore be
% initialized external to this script if doing many calculations.

%% ASK_14

%% PGA

gmm = 'ASK_14';
imt = 'PGA';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_PGA_150 = gm;

gmm = 'ASK_14';
imt = 'PGA';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_PGA_150 = gm;

gmm = 'ASK_14';
imt = 'PGA';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_PGA_150 = gm;

gmm = 'ASK_14';
imt = 'PGA';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_PGA_150 = gm;

gmm = 'ASK_14';
imt = 'PGA';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_PGA_150 = gm;

%% 0.1 Sec

gmm = 'ASK_14';
imt = 'SA0P1';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA0P1_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P1';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA0P1_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P1';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA0P1_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P1';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA0P1_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P1';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA0P1_150 = gm;

%% 0.2 Sec

gmm = 'ASK_14';
imt = 'SA0P2';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA0P2_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P2';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA0P2_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P2';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA0P2_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P2';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA0P2_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P2';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA0P2_150 = gm;

%% 0.3 Sec

gmm = 'ASK_14';
imt = 'SA0P3';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA0P3_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P3';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA0P3_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P3';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA0P3_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P3';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA0P3_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P3';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA0P3_150 = gm;

%% 0.5 Sec

gmm = 'ASK_14';
imt = 'SA0P5';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA0P5_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P5';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA0P5_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P5';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA0P5_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P5';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA0P5_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P5';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA0P5_150 = gm;

%% 0.75 Sec

gmm = 'ASK_14';
imt = 'SA0P75';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA0P75_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P75';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA0P75_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P75';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA0P75_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P75';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA0P75_150 = gm;

gmm = 'ASK_14';
imt = 'SA0P75';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA0P75_150 = gm;

%% 1 Sec

gmm = 'ASK_14';
imt = 'SA1P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA1P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA1P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA1P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA1P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA1P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA1P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA1P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA1P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA1P0_150 = gm;

%% 2 Sec

gmm = 'ASK_14';
imt = 'SA2P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA2P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA2P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA2P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA2P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA2P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA2P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA2P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA2P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA2P0_150 = gm;

%% 3 Sec

gmm = 'ASK_14';
imt = 'SA3P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA3P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA3P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA3P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA3P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA3P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA3P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA3P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA3P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA3P0_150 = gm;

%% 4 Sec

gmm = 'ASK_14';
imt = 'SA4P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA4P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA4P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA4P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA4P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA4P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA4P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA4P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA4P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA4P0_150 = gm;

%% 5 Sec

gmm = 'ASK_14';
imt = 'SA5P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M5_SA5P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA5P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M6_SA5P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA5P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7_SA5P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA5P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M7p5_SA5P0_150 = gm;

gmm = 'ASK_14';
imt = 'SA5P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_ASK_14_M8_SA5P0_150 = gm;

%% BSSA_14

%% PGA

gmm = 'BSSA_14';
imt = 'PGA';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_PGA_150 = gm;

gmm = 'BSSA_14';
imt = 'PGA';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_PGA_150 = gm;

gmm = 'BSSA_14';
imt = 'PGA';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_PGA_150 = gm;

gmm = 'BSSA_14';
imt = 'PGA';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_PGA_150 = gm;

gmm = 'BSSA_14';
imt = 'PGA';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_PGA_150 = gm;

%% 0.1 Sec

gmm = 'BSSA_14';
imt = 'SA0P1';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA0P1_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P1';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA0P1_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P1';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA0P1_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P1';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA0P1_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P1';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA0P1_150 = gm;

%% 0.2 Sec

gmm = 'BSSA_14';
imt = 'SA0P2';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA0P2_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P2';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA0P2_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P2';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA0P2_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P2';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA0P2_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P2';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA0P2_150 = gm;

%% 0.3 Sec

gmm = 'BSSA_14';
imt = 'SA0P3';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA0P3_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P3';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA0P3_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P3';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA0P3_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P3';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA0P3_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P3';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA0P3_150 = gm;

%% 0.5 Sec

gmm = 'BSSA_14';
imt = 'SA0P5';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA0P5_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P5';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA0P5_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P5';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA0P5_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P5';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA0P5_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P5';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA0P5_150 = gm;

%% 0.75 Sec

gmm = 'BSSA_14';
imt = 'SA0P75';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA0P75_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P75';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA0P75_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P75';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA0P75_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P75';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA0P75_150 = gm;

gmm = 'BSSA_14';
imt = 'SA0P75';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA0P75_150 = gm;

%% 1 Sec

gmm = 'BSSA_14';
imt = 'SA1P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA1P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA1P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA1P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA1P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA1P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA1P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA1P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA1P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA1P0_150 = gm;

%% 2 Sec

gmm = 'BSSA_14';
imt = 'SA2P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA2P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA2P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA2P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA2P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA2P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA2P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA2P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA2P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA2P0_150 = gm;

%% 3 Sec

gmm = 'BSSA_14';
imt = 'SA3P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA3P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA3P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA3P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA3P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA3P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA3P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA3P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA3P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA3P0_150 = gm;

%% 4 Sec

gmm = 'BSSA_14';
imt = 'SA4P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA4P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA4P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA4P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA4P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA4P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA4P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA4P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA4P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA4P0_150 = gm;

%% 5 Sec

gmm = 'BSSA_14';
imt = 'SA5P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M5_SA5P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA5P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M6_SA5P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA5P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7_SA5P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA5P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M7p5_SA5P0_150 = gm;

gmm = 'BSSA_14';
imt = 'SA5P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_BSSA_14_M8_SA5P0_150 = gm;

%% CB_14

%% PGA

gmm = 'CB_14';
imt = 'PGA';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_PGA_150 = gm;

gmm = 'CB_14';
imt = 'PGA';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_PGA_150 = gm;

gmm = 'CB_14';
imt = 'PGA';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_PGA_150 = gm;

gmm = 'CB_14';
imt = 'PGA';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_PGA_150 = gm;

gmm = 'CB_14';
imt = 'PGA';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_PGA_150 = gm;

%% 0.1 Sec

gmm = 'CB_14';
imt = 'SA0P1';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA0P1_150 = gm;

gmm = 'CB_14';
imt = 'SA0P1';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA0P1_150 = gm;

gmm = 'CB_14';
imt = 'SA0P1';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA0P1_150 = gm;

gmm = 'CB_14';
imt = 'SA0P1';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA0P1_150 = gm;

gmm = 'CB_14';
imt = 'SA0P1';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA0P1_150 = gm;

%% 0.2 Sec

gmm = 'CB_14';
imt = 'SA0P2';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA0P2_150 = gm;

gmm = 'CB_14';
imt = 'SA0P2';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA0P2_150 = gm;

gmm = 'CB_14';
imt = 'SA0P2';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA0P2_150 = gm;

gmm = 'CB_14';
imt = 'SA0P2';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA0P2_150 = gm;

gmm = 'CB_14';
imt = 'SA0P2';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA0P2_150 = gm;

%% 0.3 Sec

gmm = 'CB_14';
imt = 'SA0P3';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA0P3_150 = gm;

gmm = 'CB_14';
imt = 'SA0P3';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA0P3_150 = gm;

gmm = 'CB_14';
imt = 'SA0P3';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA0P3_150 = gm;

gmm = 'CB_14';
imt = 'SA0P3';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA0P3_150 = gm;

gmm = 'CB_14';
imt = 'SA0P3';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA0P3_150 = gm;

%% 0.5 Sec

gmm = 'CB_14';
imt = 'SA0P5';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA0P5_150 = gm;

gmm = 'CB_14';
imt = 'SA0P5';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA0P5_150 = gm;

gmm = 'CB_14';
imt = 'SA0P5';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA0P5_150 = gm;

gmm = 'CB_14';
imt = 'SA0P5';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA0P5_150 = gm;

gmm = 'CB_14';
imt = 'SA0P5';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA0P5_150 = gm;

%% 0.75 Sec

gmm = 'CB_14';
imt = 'SA0P75';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA0P75_150 = gm;

gmm = 'CB_14';
imt = 'SA0P75';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA0P75_150 = gm;

gmm = 'CB_14';
imt = 'SA0P75';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA0P75_150 = gm;

gmm = 'CB_14';
imt = 'SA0P75';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA0P75_150 = gm;

gmm = 'CB_14';
imt = 'SA0P75';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA0P75_150 = gm;

%% 1 Sec

gmm = 'CB_14';
imt = 'SA1P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA1P0_150 = gm;

gmm = 'CB_14';
imt = 'SA1P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA1P0_150 = gm;

gmm = 'CB_14';
imt = 'SA1P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA1P0_150 = gm;

gmm = 'CB_14';
imt = 'SA1P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA1P0_150 = gm;

gmm = 'CB_14';
imt = 'SA1P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA1P0_150 = gm;

%% 2 Sec

gmm = 'CB_14';
imt = 'SA2P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA2P0_150 = gm;

gmm = 'CB_14';
imt = 'SA2P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA2P0_150 = gm;

gmm = 'CB_14';
imt = 'SA2P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA2P0_150 = gm;

gmm = 'CB_14';
imt = 'SA2P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA2P0_150 = gm;

gmm = 'CB_14';
imt = 'SA2P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA2P0_150 = gm;

%% 3 Sec

gmm = 'CB_14';
imt = 'SA3P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA3P0_150 = gm;

gmm = 'CB_14';
imt = 'SA3P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA3P0_150 = gm;

gmm = 'CB_14';
imt = 'SA3P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA3P0_150 = gm;

gmm = 'CB_14';
imt = 'SA3P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA3P0_150 = gm;

gmm = 'CB_14';
imt = 'SA3P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA3P0_150 = gm;

%% 4 Sec

gmm = 'CB_14';
imt = 'SA4P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA4P0_150 = gm;

gmm = 'CB_14';
imt = 'SA4P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA4P0_150 = gm;

gmm = 'CB_14';
imt = 'SA4P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA4P0_150 = gm;

gmm = 'CB_14';
imt = 'SA4P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA4P0_150 = gm;

gmm = 'CB_14';
imt = 'SA4P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA4P0_150 = gm;

%% 5 Sec

gmm = 'CB_14';
imt = 'SA5P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M5_SA5P0_150 = gm;

gmm = 'CB_14';
imt = 'SA5P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M6_SA5P0_150 = gm;

gmm = 'CB_14';
imt = 'SA5P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7_SA5P0_150 = gm;

gmm = 'CB_14';
imt = 'SA5P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M7p5_SA5P0_150 = gm;

gmm = 'CB_14';
imt = 'SA5P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CB_14_M8_SA5P0_150 = gm;

%% CY_14

%% PGA

gmm = 'CY_14';
imt = 'PGA';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_PGA_150 = gm;

gmm = 'CY_14';
imt = 'PGA';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_PGA_150 = gm;

gmm = 'CY_14';
imt = 'PGA';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_PGA_150 = gm;

gmm = 'CY_14';
imt = 'PGA';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_PGA_150 = gm;

gmm = 'CY_14';
imt = 'PGA';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_PGA_150 = gm;

%% 0.1 Sec

gmm = 'CY_14';
imt = 'SA0P1';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA0P1_150 = gm;

gmm = 'CY_14';
imt = 'SA0P1';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA0P1_150 = gm;

gmm = 'CY_14';
imt = 'SA0P1';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA0P1_150 = gm;

gmm = 'CY_14';
imt = 'SA0P1';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA0P1_150 = gm;

gmm = 'CY_14';
imt = 'SA0P1';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA0P1_150 = gm;

%% 0.2 Sec

gmm = 'CY_14';
imt = 'SA0P2';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA0P2_150 = gm;

gmm = 'CY_14';
imt = 'SA0P2';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA0P2_150 = gm;

gmm = 'CY_14';
imt = 'SA0P2';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA0P2_150 = gm;

gmm = 'CY_14';
imt = 'SA0P2';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA0P2_150 = gm;

gmm = 'CY_14';
imt = 'SA0P2';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA0P2_150 = gm;

%% 0.3 Sec

gmm = 'CY_14';
imt = 'SA0P3';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA0P3_150 = gm;

gmm = 'CY_14';
imt = 'SA0P3';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA0P3_150 = gm;

gmm = 'CY_14';
imt = 'SA0P3';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA0P3_150 = gm;

gmm = 'CY_14';
imt = 'SA0P3';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA0P3_150 = gm;

gmm = 'CY_14';
imt = 'SA0P3';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA0P3_150 = gm;

%% 0.5 Sec

gmm = 'CY_14';
imt = 'SA0P5';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA0P5_150 = gm;

gmm = 'CY_14';
imt = 'SA0P5';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA0P5_150 = gm;

gmm = 'CY_14';
imt = 'SA0P5';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA0P5_150 = gm;

gmm = 'CY_14';
imt = 'SA0P5';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA0P5_150 = gm;

gmm = 'CY_14';
imt = 'SA0P5';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA0P5_150 = gm;

%% 0.75 Sec

gmm = 'CY_14';
imt = 'SA0P75';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA0P75_150 = gm;

gmm = 'CY_14';
imt = 'SA0P75';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA0P75_150 = gm;

gmm = 'CY_14';
imt = 'SA0P75';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA0P75_150 = gm;

gmm = 'CY_14';
imt = 'SA0P75';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA0P75_150 = gm;

gmm = 'CY_14';
imt = 'SA0P75';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA0P75_150 = gm;

%% 1 Sec

gmm = 'CY_14';
imt = 'SA1P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA1P0_150 = gm;

gmm = 'CY_14';
imt = 'SA1P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA1P0_150 = gm;

gmm = 'CY_14';
imt = 'SA1P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA1P0_150 = gm;

gmm = 'CY_14';
imt = 'SA1P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA1P0_150 = gm;

gmm = 'CY_14';
imt = 'SA1P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA1P0_150 = gm;

%% 2 Sec

gmm = 'CY_14';
imt = 'SA2P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA2P0_150 = gm;

gmm = 'CY_14';
imt = 'SA2P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA2P0_150 = gm;

gmm = 'CY_14';
imt = 'SA2P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA2P0_150 = gm;

gmm = 'CY_14';
imt = 'SA2P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA2P0_150 = gm;

gmm = 'CY_14';
imt = 'SA2P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA2P0_150 = gm;

%% 3 Sec

gmm = 'CY_14';
imt = 'SA3P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA3P0_150 = gm;

gmm = 'CY_14';
imt = 'SA3P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA3P0_150 = gm;

gmm = 'CY_14';
imt = 'SA3P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA3P0_150 = gm;

gmm = 'CY_14';
imt = 'SA3P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA3P0_150 = gm;

gmm = 'CY_14';
imt = 'SA3P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA3P0_150 = gm;

%% 4 Sec

gmm = 'CY_14';
imt = 'SA4P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA4P0_150 = gm;

gmm = 'CY_14';
imt = 'SA4P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA4P0_150 = gm;

gmm = 'CY_14';
imt = 'SA4P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA4P0_150 = gm;

gmm = 'CY_14';
imt = 'SA4P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA4P0_150 = gm;

gmm = 'CY_14';
imt = 'SA4P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA4P0_150 = gm;

%% 5 Sec

gmm = 'CY_14';
imt = 'SA5P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M5_SA5P0_150 = gm;

gmm = 'CY_14';
imt = 'SA5P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M6_SA5P0_150 = gm;

gmm = 'CY_14';
imt = 'SA5P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7_SA5P0_150 = gm;

gmm = 'CY_14';
imt = 'SA5P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M7p5_SA5P0_150 = gm;

gmm = 'CY_14';
imt = 'SA5P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_CY_14_M8_SA5P0_150 = gm;

%% IDRISS_14

%% PGA

gmm = 'IDRISS_14';
imt = 'PGA';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_PGA_150 = gm;

gmm = 'IDRISS_14';
imt = 'PGA';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_PGA_150 = gm;

gmm = 'IDRISS_14';
imt = 'PGA';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_PGA_150 = gm;

gmm = 'IDRISS_14';
imt = 'PGA';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_PGA_150 = gm;

gmm = 'IDRISS_14';
imt = 'PGA';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_PGA_150 = gm;

%% 0.1 Sec

gmm = 'IDRISS_14';
imt = 'SA0P1';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA0P1_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P1';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA0P1_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P1';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA0P1_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P1';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA0P1_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P1';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA0P1_150 = gm;

%% 0.2 Sec

gmm = 'IDRISS_14';
imt = 'SA0P2';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA0P2_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P2';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA0P2_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P2';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA0P2_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P2';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA0P2_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P2';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA0P2_150 = gm;

%% 0.3 Sec

gmm = 'IDRISS_14';
imt = 'SA0P3';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA0P3_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P3';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA0P3_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P3';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA0P3_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P3';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA0P3_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P3';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA0P3_150 = gm;

%% 0.5 Sec

gmm = 'IDRISS_14';
imt = 'SA0P5';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA0P5_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P5';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA0P5_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P5';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA0P5_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P5';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA0P5_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P5';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA0P5_150 = gm;

%% 0.75 Sec

gmm = 'IDRISS_14';
imt = 'SA0P75';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA0P75_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P75';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA0P75_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P75';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA0P75_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P75';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA0P75_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA0P75';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA0P75_150 = gm;

%% 1 Sec

gmm = 'IDRISS_14';
imt = 'SA1P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA1P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA1P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA1P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA1P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA1P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA1P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA1P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA1P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA1P0_150 = gm;

%% 2 Sec

gmm = 'IDRISS_14';
imt = 'SA2P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA2P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA2P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA2P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA2P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA2P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA2P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA2P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA2P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA2P0_150 = gm;

%% 3 Sec

gmm = 'IDRISS_14';
imt = 'SA3P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA3P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA3P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA3P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA3P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA3P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA3P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA3P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA3P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA3P0_150 = gm;

%% 4 Sec

gmm = 'IDRISS_14';
imt = 'SA4P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA4P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA4P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA4P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA4P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA4P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA4P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA4P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA4P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA4P0_150 = gm;

%% 5 Sec

gmm = 'IDRISS_14';
imt = 'SA5P0';
M = 5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M5_SA5P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA5P0';
M = 6;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M6_SA5P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA5P0';
M = 7;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7_SA5P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA5P0';
M = 7.5;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M7p5_SA5P0_150 = gm;

gmm = 'IDRISS_14';
imt = 'SA5P0';
M = 8;
allR = [5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500];
gm = ones(length(allR),1);

for i = 1:length(allR)
    R = allR(i);

% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 150.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model

% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]

    result = hazMat.gmmMean(gmm, imt, input);
    gm(i) = exp(result(1));
end

gm_IDRISS_14_M8_SA5P0_150 = gm;

%% Output tables for Excel

clear csvwrite

ASK_14_Median_GMMs = [allR', gm_ASK_14_M5_PGA_150, gm_ASK_14_M6_PGA_150, gm_ASK_14_M7_PGA_150, gm_ASK_14_M7p5_PGA_150, gm_ASK_14_M8_PGA_150,...
                     gm_ASK_14_M5_SA0P1_150, gm_ASK_14_M6_SA0P1_150, gm_ASK_14_M7_SA0P1_150, gm_ASK_14_M7p5_SA0P1_150, gm_ASK_14_M8_SA0P1_150,...
                     gm_ASK_14_M5_SA0P2_150, gm_ASK_14_M6_SA0P2_150, gm_ASK_14_M7_SA0P2_150, gm_ASK_14_M7p5_SA0P2_150, gm_ASK_14_M8_SA0P2_150,...
                     gm_ASK_14_M5_SA0P3_150, gm_ASK_14_M6_SA0P3_150, gm_ASK_14_M7_SA0P3_150, gm_ASK_14_M7p5_SA0P3_150, gm_ASK_14_M8_SA0P3_150,...
                     gm_ASK_14_M5_SA0P5_150, gm_ASK_14_M6_SA0P5_150, gm_ASK_14_M7_SA0P5_150, gm_ASK_14_M7p5_SA0P5_150, gm_ASK_14_M8_SA0P5_150,...
                     gm_ASK_14_M5_SA0P75_150, gm_ASK_14_M6_SA0P75_150, gm_ASK_14_M7_SA0P75_150, gm_ASK_14_M7p5_SA0P75_150, gm_ASK_14_M8_SA0P75_150,...
                     gm_ASK_14_M5_SA1P0_150, gm_ASK_14_M6_SA1P0_150, gm_ASK_14_M7_SA1P0_150, gm_ASK_14_M7p5_SA1P0_150, gm_ASK_14_M8_SA1P0_150,...
                     gm_ASK_14_M5_SA2P0_150, gm_ASK_14_M6_SA2P0_150, gm_ASK_14_M7_SA2P0_150, gm_ASK_14_M7p5_SA2P0_150, gm_ASK_14_M8_SA2P0_150,...
                     gm_ASK_14_M5_SA3P0_150, gm_ASK_14_M6_SA3P0_150, gm_ASK_14_M7_SA3P0_150, gm_ASK_14_M7p5_SA3P0_150, gm_ASK_14_M8_SA3P0_150,...
                     gm_ASK_14_M5_SA4P0_150, gm_ASK_14_M6_SA4P0_150, gm_ASK_14_M7_SA4P0_150, gm_ASK_14_M7p5_SA4P0_150, gm_ASK_14_M8_SA4P0_150,...
                     gm_ASK_14_M5_SA5P0_150, gm_ASK_14_M6_SA5P0_150, gm_ASK_14_M7_SA5P0_150, gm_ASK_14_M7p5_SA5P0_150, gm_ASK_14_M8_SA5P0_150];
                    
csvwrite('ASK_14_Median_GMMs_SS_150.csv',ASK_14_Median_GMMs)   

BSSA_14_Median_GMMs = [allR', gm_BSSA_14_M5_PGA_150, gm_BSSA_14_M6_PGA_150, gm_BSSA_14_M7_PGA_150, gm_BSSA_14_M7p5_PGA_150, gm_BSSA_14_M8_PGA_150,...
                     gm_BSSA_14_M5_SA0P1_150, gm_BSSA_14_M6_SA0P1_150, gm_BSSA_14_M7_SA0P1_150, gm_BSSA_14_M7p5_SA0P1_150, gm_BSSA_14_M8_SA0P1_150,...
                     gm_BSSA_14_M5_SA0P2_150, gm_BSSA_14_M6_SA0P2_150, gm_BSSA_14_M7_SA0P2_150, gm_BSSA_14_M7p5_SA0P2_150, gm_BSSA_14_M8_SA0P2_150,...
                     gm_BSSA_14_M5_SA0P3_150, gm_BSSA_14_M6_SA0P3_150, gm_BSSA_14_M7_SA0P3_150, gm_BSSA_14_M7p5_SA0P3_150, gm_BSSA_14_M8_SA0P3_150,...
                     gm_BSSA_14_M5_SA0P5_150, gm_BSSA_14_M6_SA0P5_150, gm_BSSA_14_M7_SA0P5_150, gm_BSSA_14_M7p5_SA0P5_150, gm_BSSA_14_M8_SA0P5_150,...
                     gm_BSSA_14_M5_SA0P75_150, gm_BSSA_14_M6_SA0P75_150, gm_BSSA_14_M7_SA0P75_150, gm_BSSA_14_M7p5_SA0P75_150, gm_BSSA_14_M8_SA0P75_150,...
                     gm_BSSA_14_M5_SA1P0_150, gm_BSSA_14_M6_SA1P0_150, gm_BSSA_14_M7_SA1P0_150, gm_BSSA_14_M7p5_SA1P0_150, gm_BSSA_14_M8_SA1P0_150,...
                     gm_BSSA_14_M5_SA2P0_150, gm_BSSA_14_M6_SA2P0_150, gm_BSSA_14_M7_SA2P0_150, gm_BSSA_14_M7p5_SA2P0_150, gm_BSSA_14_M8_SA2P0_150,...
                     gm_BSSA_14_M5_SA3P0_150, gm_BSSA_14_M6_SA3P0_150, gm_BSSA_14_M7_SA3P0_150, gm_BSSA_14_M7p5_SA3P0_150, gm_BSSA_14_M8_SA3P0_150,...
                     gm_BSSA_14_M5_SA4P0_150, gm_BSSA_14_M6_SA4P0_150, gm_BSSA_14_M7_SA4P0_150, gm_BSSA_14_M7p5_SA4P0_150, gm_BSSA_14_M8_SA4P0_150,...
                     gm_BSSA_14_M5_SA5P0_150, gm_BSSA_14_M6_SA5P0_150, gm_BSSA_14_M7_SA5P0_150, gm_BSSA_14_M7p5_SA5P0_150, gm_BSSA_14_M8_SA5P0_150];
                    
csvwrite('BSSA_14_Median_GMMs_SS_150.csv',BSSA_14_Median_GMMs)

CB_14_Median_GMMs = [allR', gm_CB_14_M5_PGA_150, gm_CB_14_M6_PGA_150, gm_CB_14_M7_PGA_150, gm_CB_14_M7p5_PGA_150, gm_CB_14_M8_PGA_150,...
                     gm_CB_14_M5_SA0P1_150, gm_CB_14_M6_SA0P1_150, gm_CB_14_M7_SA0P1_150, gm_CB_14_M7p5_SA0P1_150, gm_CB_14_M8_SA0P1_150,...
                     gm_CB_14_M5_SA0P2_150, gm_CB_14_M6_SA0P2_150, gm_CB_14_M7_SA0P2_150, gm_CB_14_M7p5_SA0P2_150, gm_CB_14_M8_SA0P2_150,...
                     gm_CB_14_M5_SA0P3_150, gm_CB_14_M6_SA0P3_150, gm_CB_14_M7_SA0P3_150, gm_CB_14_M7p5_SA0P3_150, gm_CB_14_M8_SA0P3_150,...
                     gm_CB_14_M5_SA0P5_150, gm_CB_14_M6_SA0P5_150, gm_CB_14_M7_SA0P5_150, gm_CB_14_M7p5_SA0P5_150, gm_CB_14_M8_SA0P5_150,...
                     gm_CB_14_M5_SA0P75_150, gm_CB_14_M6_SA0P75_150, gm_CB_14_M7_SA0P75_150, gm_CB_14_M7p5_SA0P75_150, gm_CB_14_M8_SA0P75_150,...
                     gm_CB_14_M5_SA1P0_150, gm_CB_14_M6_SA1P0_150, gm_CB_14_M7_SA1P0_150, gm_CB_14_M7p5_SA1P0_150, gm_CB_14_M8_SA1P0_150,...
                     gm_CB_14_M5_SA2P0_150, gm_CB_14_M6_SA2P0_150, gm_CB_14_M7_SA2P0_150, gm_CB_14_M7p5_SA2P0_150, gm_CB_14_M8_SA2P0_150,...
                     gm_CB_14_M5_SA3P0_150, gm_CB_14_M6_SA3P0_150, gm_CB_14_M7_SA3P0_150, gm_CB_14_M7p5_SA3P0_150, gm_CB_14_M8_SA3P0_150,...
                     gm_CB_14_M5_SA4P0_150, gm_CB_14_M6_SA4P0_150, gm_CB_14_M7_SA4P0_150, gm_CB_14_M7p5_SA4P0_150, gm_CB_14_M8_SA4P0_150,...
                     gm_CB_14_M5_SA5P0_150, gm_CB_14_M6_SA5P0_150, gm_CB_14_M7_SA5P0_150, gm_CB_14_M7p5_SA5P0_150, gm_CB_14_M8_SA5P0_150];
                    
csvwrite('CB_14_Median_GMMs_SS_150.csv',CB_14_Median_GMMs) 

CY_14_Median_GMMs = [allR', gm_CY_14_M5_PGA_150, gm_CY_14_M6_PGA_150, gm_CY_14_M7_PGA_150, gm_CY_14_M7p5_PGA_150, gm_CY_14_M8_PGA_150,...
                     gm_CY_14_M5_SA0P1_150, gm_CY_14_M6_SA0P1_150, gm_CY_14_M7_SA0P1_150, gm_CY_14_M7p5_SA0P1_150, gm_CY_14_M8_SA0P1_150,...
                     gm_CY_14_M5_SA0P2_150, gm_CY_14_M6_SA0P2_150, gm_CY_14_M7_SA0P2_150, gm_CY_14_M7p5_SA0P2_150, gm_CY_14_M8_SA0P2_150,...
                     gm_CY_14_M5_SA0P3_150, gm_CY_14_M6_SA0P3_150, gm_CY_14_M7_SA0P3_150, gm_CY_14_M7p5_SA0P3_150, gm_CY_14_M8_SA0P3_150,...
                     gm_CY_14_M5_SA0P5_150, gm_CY_14_M6_SA0P5_150, gm_CY_14_M7_SA0P5_150, gm_CY_14_M7p5_SA0P5_150, gm_CY_14_M8_SA0P5_150,...
                     gm_CY_14_M5_SA0P75_150, gm_CY_14_M6_SA0P75_150, gm_CY_14_M7_SA0P75_150, gm_CY_14_M7p5_SA0P75_150, gm_CY_14_M8_SA0P75_150,...
                     gm_CY_14_M5_SA1P0_150, gm_CY_14_M6_SA1P0_150, gm_CY_14_M7_SA1P0_150, gm_CY_14_M7p5_SA1P0_150, gm_CY_14_M8_SA1P0_150,...
                     gm_CY_14_M5_SA2P0_150, gm_CY_14_M6_SA2P0_150, gm_CY_14_M7_SA2P0_150, gm_CY_14_M7p5_SA2P0_150, gm_CY_14_M8_SA2P0_150,...
                     gm_CY_14_M5_SA3P0_150, gm_CY_14_M6_SA3P0_150, gm_CY_14_M7_SA3P0_150, gm_CY_14_M7p5_SA3P0_150, gm_CY_14_M8_SA3P0_150,...
                     gm_CY_14_M5_SA4P0_150, gm_CY_14_M6_SA4P0_150, gm_CY_14_M7_SA4P0_150, gm_CY_14_M7p5_SA4P0_150, gm_CY_14_M8_SA4P0_150,...
                     gm_CY_14_M5_SA5P0_150, gm_CY_14_M6_SA5P0_150, gm_CY_14_M7_SA5P0_150, gm_CY_14_M7p5_SA5P0_150, gm_CY_14_M8_SA5P0_150];
                    
csvwrite('CY_14_Median_GMMs_SS_150.csv',CY_14_Median_GMMs)

IDRISS_14_Median_GMMs = [allR', gm_IDRISS_14_M5_PGA_150, gm_IDRISS_14_M6_PGA_150, gm_IDRISS_14_M7_PGA_150, gm_IDRISS_14_M7p5_PGA_150, gm_IDRISS_14_M8_PGA_150,...
                     gm_IDRISS_14_M5_SA0P1_150, gm_IDRISS_14_M6_SA0P1_150, gm_IDRISS_14_M7_SA0P1_150, gm_IDRISS_14_M7p5_SA0P1_150, gm_IDRISS_14_M8_SA0P1_150,...
                     gm_IDRISS_14_M5_SA0P2_150, gm_IDRISS_14_M6_SA0P2_150, gm_IDRISS_14_M7_SA0P2_150, gm_IDRISS_14_M7p5_SA0P2_150, gm_IDRISS_14_M8_SA0P2_150,...
                     gm_IDRISS_14_M5_SA0P3_150, gm_IDRISS_14_M6_SA0P3_150, gm_IDRISS_14_M7_SA0P3_150, gm_IDRISS_14_M7p5_SA0P3_150, gm_IDRISS_14_M8_SA0P3_150,...
                     gm_IDRISS_14_M5_SA0P5_150, gm_IDRISS_14_M6_SA0P5_150, gm_IDRISS_14_M7_SA0P5_150, gm_IDRISS_14_M7p5_SA0P5_150, gm_IDRISS_14_M8_SA0P5_150,...
                     gm_IDRISS_14_M5_SA0P75_150, gm_IDRISS_14_M6_SA0P75_150, gm_IDRISS_14_M7_SA0P75_150, gm_IDRISS_14_M7p5_SA0P75_150, gm_IDRISS_14_M8_SA0P75_150,...
                     gm_IDRISS_14_M5_SA1P0_150, gm_IDRISS_14_M6_SA1P0_150, gm_IDRISS_14_M7_SA1P0_150, gm_IDRISS_14_M7p5_SA1P0_150, gm_IDRISS_14_M8_SA1P0_150,...
                     gm_IDRISS_14_M5_SA2P0_150, gm_IDRISS_14_M6_SA2P0_150, gm_IDRISS_14_M7_SA2P0_150, gm_IDRISS_14_M7p5_SA2P0_150, gm_IDRISS_14_M8_SA2P0_150,...
                     gm_IDRISS_14_M5_SA3P0_150, gm_IDRISS_14_M6_SA3P0_150, gm_IDRISS_14_M7_SA3P0_150, gm_IDRISS_14_M7p5_SA3P0_150, gm_IDRISS_14_M8_SA3P0_150,...
                     gm_IDRISS_14_M5_SA4P0_150, gm_IDRISS_14_M6_SA4P0_150, gm_IDRISS_14_M7_SA4P0_150, gm_IDRISS_14_M7p5_SA4P0_150, gm_IDRISS_14_M8_SA4P0_150,...
                     gm_IDRISS_14_M5_SA5P0_150, gm_IDRISS_14_M6_SA5P0_150, gm_IDRISS_14_M7_SA5P0_150, gm_IDRISS_14_M7p5_SA5P0_150, gm_IDRISS_14_M8_SA5P0_150];
                    
csvwrite('IDRISS_14_Median_GMMs_SS_150.csv',IDRISS_14_Median_GMMs)



