function [gm_log,gm_lin]=calc_BSSA_R60_M6p6(vs30,Z1,period)

%
%nshmpHaz = '/Users/ashumway/USGS_Computer_Codes/nshmp-haz/nshmp-haz/dist/nshmp-haz.jar';
nshmpHaz = '/Users/mmoschetti/progs/NSHMP_dir/GIT/nshmp-haz/build/libs/nshmp-haz.jar';

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

%import org.opensha2.etc.*
import gov.usgs.earthquake.nshmp.etc.*
% =========================================================================
% Single model ground motion calculation:

% Initialize calculator:
hazMat = HazMat.init(nshmpHaz);

% Note that hazMat is stateless and reusable and should therefore be
% initialized external to this script if doing many calculations.

gmm = 'BSSA_14';
%gmm = 'CY_14';
%M = 5;
%R=30;
M = 6.6;
R=60;

% assign period
if abs(period-1)<0.001
  imt = 'SA1P0';
elseif abs(period-2)<0.001
  imt = 'SA2P0';
elseif abs(period-3)<0.001
  imt = 'SA3P0';
elseif abs(period-5)<0.001
  imt = 'SA5P0';
elseif abs(period-10)<0.001
  imt = 'SA10P0';
elseif abs(period-0.1)<0.001
  imt = 'SA0P1';
elseif abs(period-0)<0.001
  imt = 'PGA';
else
  error('No calculation available for this period')
end

% GMM for input vs30 and z1
% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
%input.vs30  = 760.0; % in m/s
input.vs30  = vs30; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
%input.z1p0  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   Z1; % in km; NaN triggers default basin depth model
%input.z1p0  =   NaN; % in km; NaN triggers default basin depth model
% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]
result = hazMat.gmmMean(gmm, imt, input);
gm_vs30_z1 = result(1);

% GMM for input vs30 and z1
% Set up a GMM input parameter object. These data are a source and site
% parameterization that will satisfy all currently implemented Gmms. Note
% that not all models will necessarily use all parameters.
input = GmmParams();
input.Mw    =   M; % moment magnitude
input.rJB   =   R; % Joyner-Boore distance
input.rRup  =   R; % distance to closest point on rupture surface
input.rX    =   R; % distance from source trace; hanging (+); foot (-) wall
input.dip   =  90.0; % in degrees
input.width =  10.0; % in km
input.zTop  =   5.0; % in km
input.zHyp  =   8.0; % in km
input.rake  =   0.0; % in degrees
input.vs30  = 760.0; % in m/s
input.vsInf =  true; % boolean
input.z2p5  =   NaN; % in km; NaN triggers default basin depth model
input.z1p0  =   NaN; % in km; NaN triggers default basin depth model
% Do a calculation. The MatUtil.calc(gmm, imt, gmmInput) method returns an
% array of [ln(median ground motion), sigma]
result = hazMat.gmmMean(gmm, imt, input);
gm_760_def = result(1);
gm_log=gm_760_def;
gm_lin=exp(gm_log);

% ratio of ground motions (linear)
%gm_ratio_vs30_z1_to_default=exp(gm_vs30_z1-gm_760_def);



end
%----------------------------------------------------
