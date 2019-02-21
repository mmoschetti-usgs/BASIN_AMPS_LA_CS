function []=compare_amps_allPer()

%
close all

%
figure(1)
figure(2)
sval=25;
cntv=0;


% loop over periods
for perIn=[2 3 5 10]
% set file names
  fileNm_adjB=sprintf('B_%d.00_BA14_Vs30_Zx_adj.csv',perIn);
  fileNm_ampBSSA=sprintf('ampBA_varVs30_760_%ds.csv',perIn);
  fileNm_BmapsCS=sprintf('CS-LA15.4_%d.00_Bs.dat',perIn);
% read file data
  [lon,lat,z1,z2p5,vs30_wills,B_adj760,B_varVs30]=read_adjusted_Bmaps_BSSA(fileNm_adjB);
  [lon_gm,lat_gm,z1_gm,z2p5_gm,vs30_wills_gm,lnAmp_BSSA_760]=read_ampBSSA(fileNm_ampBSSA);
  [lon1_cs,lat_cs,B_CS]=read_Bmaps_CS(fileNm_BmapsCS);
% compute amplifications and 
  cnt=1;
  cnt1=1;
  for ii=1:length(lon)
    if ~isnan(vs30_wills(ii))
      lon_arr(cnt)=lon(ii);
      lat_arr(cnt)=lat(ii);
      z1_arr(cnt)=z1(ii);
      z2p5_arr(cnt)=z2p5(ii);
      amp_CS_BA760(cnt)=B_CS(ii)-B_adj760(ii);
      amp_CS_BAvarVs30(cnt)=B_CS(ii)-B_varVs30(ii);
%      fprintf(fid,'%.4f,%.4f,%.1f,%.1f,%.1f,%.4f\n',lon(ii),lat(ii),vs30_wills(ii),z2p5(ii),z1(ii),amp_CS_BA760(cnt));
      if z2p5(ii)>2
        lon_arr_basin(cnt1)=lon_arr(cnt);
        lat_arr_basin(cnt1)=lat_arr(cnt);
        z1_arr_basin(cnt1)=z1_arr(cnt);
        z2p5_arr_basin(cnt1)=z2p5_arr(cnt);
        amp_CS_BA760_basin(cnt1)=amp_CS_BA760(cnt);
        amp_CS_BAvarVs30_basin(cnt1)=amp_CS_BAvarVs30(cnt);
        cnt1=cnt1+1;
      end
      cnt=cnt+1;
    end
  end
%
  lon_arr=lon_arr';
  lat_arr=lat_arr';
  z1_arr=z1_arr';
  z2p5_arr=z2p5_arr';
  amp_CS_BA760=amp_CS_BA760';
  amp_CS_BAvarVs30=amp_CS_BAvarVs30';
%
  lon_arr_basin=lon_arr_basin';
  lat_arr_basin=lat_arr_basin';
  z1_arr_basin=z1_arr_basin';
  z2p5_arr_basin=z2p5_arr_basin';
  amp_CS_BA760_basin=amp_CS_BA760_basin';
  amp_CS_BAvarVs30_basin=amp_CS_BAvarVs30_basin';
%
  figure(1)
  subplot(3,4,1+cntv)
  scatter(lon_arr,lat_arr,sval,amp_CS_BA760,'filled');
  title(sprintf('T=%d s',perIn))
  ylabel('ln(CS/BSSA_{Vs30})')
  subplot(3,4,5+cntv)
  scatter(lon_arr,lat_arr,sval,amp_CS_BAvarVs30,'filled');
  ylabel('ln(CS/BSSA_{Vs30=760})')
  subplot(3,4,9+cntv)
  plot(z1_arr,amp_CS_BAvarVs30,'bs');
  ylabel('ln(CS/BSSA_{Vs30})')
  xlabel('Z1 (km)')
%
  figure(2)
  subplot(3,4,1+cntv)
  scatter(lon_arr_basin,lat_arr_basin,sval,amp_CS_BA760_basin,'filled');
  title(sprintf('T=%d s',perIn))
  ylabel('ln(CS/BSSA_{Vs30})')
  subplot(3,4,5+cntv)
  scatter(lon_arr_basin,lat_arr_basin,sval,amp_CS_BAvarVs30_basin,'filled');
  ylabel('ln(CS/BSSA_{Vs30=760})')
  subplot(3,4,9+cntv)
  plot(z1_arr_basin,amp_CS_BAvarVs30_basin,'bs');
  ylabel('ln(CS/BSSA_{Vs30})')
  xlabel('Z1 (km)')

%
  cntv=cntv+1;
end
