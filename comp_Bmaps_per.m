function []= comp_Bmaps_per()

%
close all

% 
%fileOrig=sprintf('B_%.2f_BA14_Vs30_Zx.csv',PerIn);

% read and plot
% T=2 s
PerIn=2;
fileMod=sprintf('B_%.2f_BA14_Vs30_Zx_adj.csv',PerIn);
[lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileMod,1);
% T=3 s
PerIn=3;
fileMod=sprintf('B_%.2f_BA14_Vs30_Zx_adj.csv',PerIn);
[lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileMod,2);
% T=5 s
PerIn=5;
fileMod=sprintf('B_%.2f_BA14_Vs30_Zx_adj.csv',PerIn);
[lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileMod,3);
% T=10 s
PerIn=10;
fileMod=sprintf('B_%.2f_BA14_Vs30_Zx_adj.csv',PerIn);
[lon,lat,vs30,z2p5,z1,Bval_orig,Bval_760]=read_mod_file(fileMod,4);

%
figure(1)
plot_periods()
print('-f','-dpng','-r300','pl_diff_Bmaps_Vs30.png')
figure(2)
plot_periods()
print('-f','-dpng','-r300','pl_diff_Bmaps_Z1.png')
figure(3)
plot_periods()
print('-f','-dpng','-r300','pl_diff_Bmaps_maps.png')


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
function []=read_orig_file(fileIn)

%
disp(sprintf('Reading from original file, %s', fileIn))

%
ff=csvread(fileIn,0,0);
lon=ff(:,1);
lat=ff(:,2);
Bval=ff(:,3);
vs30_wills2015=ff(:,4);
vs30_wills2006=ff(:,5);
z2p5=ff(:,6);
z1=ff(:,7);

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

%
Bdiff=Bval_760-Bval_orig;

% plotting
% Vs30
sval=40
figure(1)
subplot(2,2,Nsubplot)
scatter(vs30,Bdiff,sval,z1,'filled')
xlabel('Vs30 (m/s)')
ylabel('ln(760/Vs30)')
colorbar
axis([150 800 -2 0.1])
% Z1
figure(2)
subplot(2,2,Nsubplot)
scatter(z1,Bdiff,sval,vs30,'filled')
xlabel('Z1 (km)')
ylabel('ln(760/Vs30)')
colorbar
axis([0 1.5 -2 0.1])
% maps
figure(3)
subplot(2,2,Nsubplot)
scatter(lon,lat,sval,Bdiff,'filled')
xlabel('ln(760/Vs30)')
colorbar


end
%-----------------------------------------------------
