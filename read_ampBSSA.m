%-----------------------------------------------
function [lon,lat,z1,z2p5,vs30_wills,ln_Amp760]=read_ampBSSA(fileIn)

%
% lon,lat,Vs30_Wills,Z2.5,Z1,ln_Amp760
ff=readtable(fileIn);
lon=ff.lon;
lat=ff.lat;
vs30_wills=ff.Vs30_Wills;
%z2p5=ff.Z2p5;
z2p5=ff.Z2_5;
z1=ff.Z1;
ln_Amp760=ff.ln_Amp760;

end
%-----------------------------------------------

