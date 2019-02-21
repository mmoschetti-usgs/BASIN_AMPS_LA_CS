%-----------------------------------------------
function [lon,lat,z1,z2p5,vs30_wills,B_adj760,B_varVs30]=read_adjusted_Bmaps_BSSA(fileIn)

%
% lon,lat,B_vs30_var,vs30_wills,Z2p5,Z1,B_vs30_760
ff=readtable(fileIn);
%lon=ff.lon;
%lat=ff.lat;
%z1=ff.Z1;
%z2p5=ff.Z2p5;
%vs30_wills=ff.vs30_wills;
%B_varVs30=ff.B_vs30_var;
%B_adj760=ff.B_vs30_760;
lon=ff.Var1;
lat=ff.Var2;
z1=ff.Var6;
z2p5=ff.Var5;
vs30_wills=ff.Var4;
B_varVs30=ff.Var3;
B_adj760=ff.Var7;

end
%-----------------------------------------------
