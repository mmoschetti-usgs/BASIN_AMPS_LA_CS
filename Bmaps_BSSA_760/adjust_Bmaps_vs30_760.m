function []=adjust_Bmaps_vs30_760(inputF_Bmap,period)

% 
outFile=strrep(inputF_Bmap,'Zx','Zx_adj');

% read all data 
ff=csvread(inputF_Bmap);
lon=ff(:,1);
lat=ff(:,2);
Bvalue=ff(:,3);
Vs30_wills2015=ff(:,4);
Vs30_wills2006=ff(:,5);
Z2p5=ff(:,6);
Z1=ff(:,7);
%size(ff)

% loop over all B-values and adjust from site-specific values to Vs30=760 m/s
% default Z1=0.041 (BSSA, Vs30=760)
fid=fopen(outFile,'w');
fprintf(fid,'lon,lat,B(vs30,variable),Vs30-Wills,Z2.5,Z1,B(vs30=760m/s)\n');
for cnt=1:length(lon)
%for cnt=1:50
  vs30_site=Vs30_wills2015(cnt);
  Z1_site=Z1(cnt);
  if isnan(vs30_site)
    ln_amp=0;
  else
    [lnGM_vs30_z1,lnGM_760_def,gm_ratio]=calc_BSSA_Vs30_Z1_input(vs30_site,Z1_site,period);
    ln_amp=lnGM_760_def-lnGM_vs30_z1;
  end
  modBvalue=Bvalue(cnt)+ln_amp;
  if isnan(vs30_site)
    fprintf(fid,'%.4f,%.4f,%.5f,%s,%.1f,%.1f,%.5f\n',lon(cnt),lat(cnt),Bvalue(cnt),Vs30_wills2015(cnt),Z2p5(cnt),Z1(cnt),modBvalue);
  else
    fprintf(fid,'%.4f,%.4f,%.5f,%.1f,%.1f,%.1f,%.5f\n',lon(cnt),lat(cnt),Bvalue(cnt),Vs30_wills2015(cnt),Z2p5(cnt),Z1(cnt),modBvalue);
  end
  
end

%
fclose(fid);
disp(sprintf('Wrote to file, %s', outFile));

% quick check
ff=csvread(outFile,1,0);
lon=ff(:,1);
lat=ff(:,2);
B_orig=ff(:,3);
vs30=ff(:,4);
z1=ff(:,5);
z2p5=ff(:,6);
B_adj_760=ff(:,7);

%
figure
subplot(1,2,1)
plot(B_orig,B_adj_760,'bs')
xlabel('B (variable Vs30)')
ylabel('B (Vs30=760 m/s)')
title(sprintf('Adjusted-B, %d s', period))
%
subplot(1,2,2)
plot(vs30,B_adj_760-B_orig,'bs')
xlabel('Vs30 (m/s)')
ylabel('\Delta{B} (B,760-B,vs30)')
%
prName=sprintf('pl_adj_B_%d.png',period);
print('-f','-dpng','-r300',prName)




end
%--------------------------------------------------
