function []=compare_amps()

%
close all

% read files
[lon,lat,z1,z2p5,vs30_wills,B_adj760,B_varVs30]=read_adjusted_Bmaps_BSSA('B_3.00_BA14_Vs30_Zx_adj.csv');
[lon_gm,lat_gm,z1_gm,z2p5_gm,vs30_wills_gm,lnAmp_BSSA_760]=read_ampBSSA('ampBA_varVs30_760_3s.csv');
[lon1_cs,lat_cs,B_CS]=read_Bmaps_CS('CS-LA15.4_3.00_Bs.dat');

sum_BSSA_B=sum(B_varVs30)
sum_BSSA_B760=sum(B_adj760)
sum_CS_B=sum(B_CS)

%figure
%plot(vs30_wills,'bs')

% loop through points and extract those with valid Vs30 values
%length(lon)
cnt=1;
fid=fopen('ampCS_varVs30_760_3s.csv','w');
fprintf(fid,'lon,lat,Vs30_Wills,Z2p5,Z1,ln_Amp760\n');
for ii=1:length(lon)
  if ~isnan(vs30_wills(ii))
%    disp(sprintf('%s',vs30_wills(ii)))
%    disp(sprintf('%.1f',vs30_wills(ii)))
    lon_arr(cnt)=lon(ii);
    lat_arr(cnt)=lat(ii);
    z1_arr(cnt)=z1(ii);
    amp_CS_BA760(cnt)=B_CS(ii)-B_adj760(ii);
%    amp_CS_BSSA(cnt)=amp_CS_BA760(cnt)/lnAmp_BSSA_760(ii);
    fprintf(fid,'%.4f,%.4f,%.1f,%.1f,%.1f,%.4f\n',lon(ii),lat(ii),vs30_wills(ii),z2p5(ii),z1(ii),amp_CS_BA760(cnt));
    cnt=cnt+1;
  end
end
lon_arr=lon_arr';
lat_arr=lat_arr';
amp_CS_BA760=amp_CS_BA760';
disp(cnt)
%whos

%
sval=25;
figure
subplot(2,2,1)
scatter(lon_gm,lat_gm,sval,lnAmp_BSSA_760,'filled')
title('BSSA,amp,relative to BSSA_{760}')
colorbar
subplot(2,2,2)
scatter(lon_arr,lat_arr,sval,amp_CS_BA760,'filled')
title('CS,amp,relative to BSSA_{760}')
colorbar
subplot(2,2,3)
scatter(lon_arr,lat_arr,sval,(amp_CS_BA760-lnAmp_BSSA_760),'filled')
title('CS,amp,relative to BSSA,amp')
colorbar
subplot(2,2,4)
plot(z1_arr,(amp_CS_BA760-lnAmp_BSSA_760),'bs')
ylabel('CS,amp,relative to BSSA,amp')
xlabel('Z1 (km)')
%scatter(lon_arr,lat_arr,sval,amp_CS_BSSA,'filled')

fig2=0
if fig2
figure
subplot(1,2,1)
plot(lon_gm,lon_arr,'bs')
subplot(1,2,2)
plot(lat_gm,lat_arr,'bs')
end

end
%-----------------------------------------------

%-----------------------------------------------
function [lon,lat,Bval]=read_Bmaps_CS(fileIn)

%
ff=load(fileIn);
lon=ff(:,1);
lat=ff(:,2);
Bval=ff(:,3);

end
%-----------------------------------------------

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
