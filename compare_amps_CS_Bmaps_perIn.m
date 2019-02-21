function []=compare_amps_CS_Bmaps_perIn(perIn)

%
close all

% write_file_names
fileNm_adjB=sprintf('B_%d.00_BA14_Vs30_Zx_adj.csv',perIn)
fileNm_ampBSSA=sprintf('ampBA_varVs30_760_%ds.csv',perIn)
fileNm_BmapsCS=sprintf('CS-LA15.4_%d.00_Bs.dat',perIn)
fileNm_ampCS=sprintf('ampCS_varVs30_760_%ds.csv',perIn)
plotNm=sprintf('pl_amps_CS_Bmaps_%ds.png', perIn);

% read files
[lon,lat,z1,z2p5,vs30_wills,B_adj760,B_varVs30]=read_adjusted_Bmaps_BSSA(fileNm_adjB);
[lon_gm,lat_gm,z1_gm,z2p5_gm,vs30_wills_gm,lnAmp_BSSA_760]=read_ampBSSA(fileNm_ampBSSA);
[lon1_cs,lat_cs,B_CS]=read_Bmaps_CS(fileNm_BmapsCS);

sum_BSSA_B=sum(B_varVs30)
sum_BSSA_B760=sum(B_adj760)
sum_CS_B=sum(B_CS)

%figure
%plot(vs30_wills,'bs')

% loop through points and extract those with valid Vs30 values
%length(lon)
cnt=1;
fid=fopen(fileNm_ampCS,'w');
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
title(sprintf('ln(BSSA,amp/BSSA_{760}), %d s',perIn))
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
%
print('-f','-dpng','-r300',plotNm)

%
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




