function []= comp_Bmaps_per()

%
close all

% 
%fileOrig=sprintf('B_%.2f_BA14_Vs30_Zx.csv',PerIn);

% read and plot
% T=2 s
PerIn=2;
% files adjusted, original
fileMod=sprintf('B_%.2f_BA14_Vs30_Zx_adj.csv',PerIn);
[lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileMod,1);
% files sent to JB
file1=sprintf('B_2.00_BA14_Vs30_760.csv');
[lon,lat,Bval_760]=read_format2_file(file1);
% files returned from JB
file2=sprintf('BSSA_760_siteAmp.csv');
[lon,lat,Bval_760]=read_format3_file(file2);

%
%figure(1)
%plot_periods()
%print('-f','-dpng','-r300','pl_diff_Bmaps_Vs30.png')
%figure(2)
%plot_periods()
%print('-f','-dpng','-r300','pl_diff_Bmaps_Z1.png')
%figure(3)
%plot_periods()
%print('-f','-dpng','-r300','pl_diff_Bmaps_maps.png')


end
%-----------------------------------------------------


%-----------------------------------------------------
function []=plot_periods()

%
subplot(2,2,1)
title('T=2 s')
subplot(2,2,2)
title('T=3 s')
subplot(2,2,3)
title('T=5 s')
subplot(2,2,4)
title('T=10 s')

end
%-----------------------------------------------------


%-----------------------------------------------------
function [lon,lat,Bval]=read_format3_file(fileIn)

%
disp(sprintf('Reading from modified Vs30=760 file, %s', fileIn))

%
ff=csvread(fileIn,1,2);
lon=ff(:,1);
lat=ff(:,2);
Bval=ff(:,9);

Bval(1:10)

figure(1)
sval=40;
subplot(1,3,3)
scatter(lon,lat,sval,Bval,'filled')
caxis([-2.2 0.8])
colorbar
%
figure(2)
hold on
sval=30
scatter(lon,lat,sval,Bval,'filled')
scatter(lon,lat,sval,-2.5*ones(length(lon),1))


end
%-----------------------------------------------------



%-----------------------------------------------------
function [lon,lat,Bval]=read_format2_file(fileIn)

%
disp(sprintf('Reading from modified Vs30=760 file, %s', fileIn))

%
ff=csvread(fileIn,1,0);
lon=ff(:,1);
lat=ff(:,2);
Bval=ff(:,3);
min(Bval)
max(Bval)


figure(1)
sval=40;
subplot(1,3,2)
scatter(lon,lat,sval,Bval,'filled')
caxis([-2.2 0.8])
colorbar
%
figure(2)
sval=50;
scatter(lon,lat,sval,Bval,'filled')
caxis([-2.2 0.8])
colorbar


end
%-----------------------------------------------------

%-----------------------------------------------------
function [lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileIn,Nsubplot)

%
disp(sprintf('Reading from modified file, %s', fileIn))

%
ff=csvread(fileIn,1,0);
lon=ff(:,1);
lat=ff(:,2);
Bval_orig=ff(:,3);
%vs30_wills2015=ff(:,4);
vs30=ff(:,4);
z2p5=ff(:,5);
z1=ff(:,6);
Bval_760=ff(:,7);
min(Bval_760)
max(Bval_760)

%
Bdiff=Bval_760-Bval_orig;

% plotting
% Vs30
sval=40
figure(1)
subplot(1,3,Nsubplot)
scatter(lon,lat,sval,Bval_760,'filled')
caxis([-2.2 0.8])
colorbar
%axis([150 800 -2 0.1])
% Z1
%figure(2)
%subplot(2,2,Nsubplot)
%scatter(z1,Bdiff,sval,vs30,'filled')
%xlabel('Z1 (km)')
%ylabel('ln(760/Vs30)')
%colorbar
%axis([0 1.5 -2 0.1])
%% maps
%figure(3)
%subplot(2,2,Nsubplot)
%scatter(lon,lat,sval,Bdiff,'filled')
%xlabel('ln(760/Vs30)')
%colorbar


end
%-----------------------------------------------------
