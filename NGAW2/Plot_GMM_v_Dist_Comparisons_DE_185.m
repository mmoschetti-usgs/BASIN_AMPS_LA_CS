%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Name: Plot_GMM_v_Distance_Comparisons_DE_185.m
%
% Author: Allison M. Shumway (ashumway@usgs.gov)
%
% Date Created: 03/03/16
%
% Last Modified: 05/15/17
%
% Purpose: This function is used to plot hazard curves for a single site, or can be used
% to compare multiple hazard curves for comparison. MODIFIED FOR GMM QA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import Data from Excel Spreadsheets
    
    addpath('/Users/ashumway/Documents/Projects/Multi-Point_Spectra/GMM_QA_2014_NSHM_nshmp-haz/NGAW2/nshmp-haz/')
    
    nshmp_haz_data_ASK_14 = csvread('ASK_14_Median_GMMs_SS_185.csv');
    R = nshmp_haz_data_ASK_14(1:19,1);
    nshmp_haz_data_BSSA_14 = csvread('BSSA_14_Median_GMMs_SS_185.csv');
    nshmp_haz_data_CB_14 = csvread('CB_14_Median_GMMs_SS_185.csv');
    nshmp_haz_data_CY_14 = csvread('CY_14_Median_GMMs_SS_185.csv');
    nshmp_haz_data_IDRISS_14 = csvread('IDRISS_14_Median_GMMs_SS_185.csv');

    
%% Plot GMM vs. Distance comparisons of WUS GMMs implemented in nshmp-haz for 5 magnitudes and 11 periods

% NOTE: Idriss14 is not valid at Site Class DE!

%% M5

% PGA
    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,2), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,2), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,2), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,2), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,2), 'k-','LineWidth',1); 
    ylabel('Median Peak Ground Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_PGA_185_GM_vs_Distance_Plot_M5','png')

% 0.1 Sec (10 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,7), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,7), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,7), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,7), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,7), 'k-','LineWidth',1); 
    ylabel('Median 0.1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P1_185_GM_vs_Distance_Plot_M5','png')    

% 0.2 Sec (5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,12), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,12), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,12), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,12), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,12), 'k-','LineWidth',1); 
    ylabel('Median 0.2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P2_185_GM_vs_Distance_Plot_M5','png')      
    
% 0.3 Sec (3.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,17), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,17), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,17), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,17), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,17), 'k-','LineWidth',1); 
    ylabel('Median 0.3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P3_185_GM_vs_Distance_Plot_M5','png')          
    
% 0.5 Sec (2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,22), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,22), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,22), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,22), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,22), 'k-','LineWidth',1); 
    ylabel('Median 0.5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P5_185_GM_vs_Distance_Plot_M5','png')     

% 0.75 Sec (1.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,27), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,27), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,27), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,27), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,27), 'k-','LineWidth',1); 
    ylabel('Median 0.75 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P75_185_GM_vs_Distance_Plot_M5','png')       
    
% 1 Sec (1 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,32), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,32), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,32), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,32), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,32), 'k-','LineWidth',1); 
    ylabel('Median 1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA1P0_185_GM_vs_Distance_Plot_M5','png')      
    
% 2 Sec (0.5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,37), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,37), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,37), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,37), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,37), 'k-','LineWidth',1); 
    ylabel('Median 2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA2P0_185_GM_vs_Distance_Plot_M5','png')     
    
% 3 Sec (0.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,42), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,42), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,42), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,42), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,42), 'k-','LineWidth',1); 
    ylabel('Median 3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA3P0_185_GM_vs_Distance_Plot_M5','png')     
        
% 4 Sec (0.25 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,47), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,47), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,47), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,47), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,47), 'k-','LineWidth',1); 
    ylabel('Median 4 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA4P0_185_GM_vs_Distance_Plot_M5','png')    
    
% 5 Sec (0.2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,52), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,52), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,52), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,52), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,52), 'k-','LineWidth',1); 
    ylabel('Median 5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA5P0_185_GM_vs_Distance_Plot_M5','png')   
    
%% M6

% PGA
    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,3), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,3), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,3), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,3), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,3), 'k-','LineWidth',1); 
    ylabel('Median Peak Ground Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_PGA_185_GM_vs_Distance_Plot_M6','png')

% 0.1 Sec (10 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,8), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,8), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,8), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,8), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,8), 'k-','LineWidth',1); 
    ylabel('Median 0.1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P1_185_GM_vs_Distance_Plot_M6','png')    

% 0.2 Sec (5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,13), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,13), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,13), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,13), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,13), 'k-','LineWidth',1); 
    ylabel('Median 0.2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P2_185_GM_vs_Distance_Plot_M6','png')      
    
% 0.3 Sec (3.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,18), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,18), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,18), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,18), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,18), 'k-','LineWidth',1); 
    ylabel('Median 0.3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P3_185_GM_vs_Distance_Plot_M6','png')          
    
% 0.5 Sec (2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,23), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,23), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,23), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,23), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,23), 'k-','LineWidth',1); 
    ylabel('Median 0.5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P5_185_GM_vs_Distance_Plot_M6','png')     

% 0.75 Sec (1.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,28), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,28), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,28), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,28), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,28), 'k-','LineWidth',1); 
    ylabel('Median 0.75 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P75_185_GM_vs_Distance_Plot_M6','png')       
    
% 1 Sec (1 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,33), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,33), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,33), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,33), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,33), 'k-','LineWidth',1); 
    ylabel('Median 1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA1P0_185_GM_vs_Distance_Plot_M6','png')      
    
% 2 Sec (0.5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,38), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,38), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,38), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,38), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,38), 'k-','LineWidth',1); 
    ylabel('Median 2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA2P0_185_GM_vs_Distance_Plot_M6','png')     
    
% 3 Sec (0.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,43), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,43), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,43), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,43), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,43), 'k-','LineWidth',1); 
    ylabel('Median 3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA3P0_185_GM_vs_Distance_Plot_M6','png')     
        
% 4 Sec (0.25 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,48), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,48), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,48), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,48), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,48), 'k-','LineWidth',1); 
    ylabel('Median 4 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA4P0_185_GM_vs_Distance_Plot_M6','png')    
    
% 5 Sec (0.2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,53), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,53), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,53), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,53), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,53), 'k-','LineWidth',1); 
    ylabel('Median 5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M6 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA5P0_185_GM_vs_Distance_Plot_M6','png') 
    
%% M7

% PGA
    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,4), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,4), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,4), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,4), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,4), 'k-','LineWidth',1); 
    ylabel('Median Peak Ground Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_PGA_185_GM_vs_Distance_Plot_M7','png')

% 0.1 Sec (10 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,9), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,9), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,9), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,9), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,9), 'k-','LineWidth',1); 
    ylabel('Median 0.1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P1_185_GM_vs_Distance_Plot_M7','png')    

% 0.2 Sec (5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,14), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,14), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,14), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,14), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,14), 'k-','LineWidth',1); 
    ylabel('Median 0.2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P2_185_GM_vs_Distance_Plot_M7','png')      
    
% 0.3 Sec (3.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,19), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,19), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,19), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,19), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,19), 'k-','LineWidth',1); 
    ylabel('Median 0.3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P3_185_GM_vs_Distance_Plot_M7','png')          
    
% 0.5 Sec (2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,24), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,24), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,24), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,24), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,24), 'k-','LineWidth',1); 
    ylabel('Median 0.5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P5_185_GM_vs_Distance_Plot_M7','png')     

% 0.75 Sec (1.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,29), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,29), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,29), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,29), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,29), 'k-','LineWidth',1); 
    ylabel('Median 0.75 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P75_185_GM_vs_Distance_Plot_M7','png')       
    
% 1 Sec (1 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,34), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,34), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,34), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,34), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,34), 'k-','LineWidth',1); 
    ylabel('Median 1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA1P0_185_GM_vs_Distance_Plot_M7','png')      
    
% 2 Sec (0.5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,39), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,39), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,39), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,39), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,39), 'k-','LineWidth',1); 
    ylabel('Median 2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA2P0_185_GM_vs_Distance_Plot_M7','png')     
    
% 3 Sec (0.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,44), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,44), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,44), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,44), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,44), 'k-','LineWidth',1); 
    ylabel('Median 3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA3P0_185_GM_vs_Distance_Plot_M7','png')     
        
% 4 Sec (0.25 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,49), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,49), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,49), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,49), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,49), 'k-','LineWidth',1); 
    ylabel('Median 4 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA4P0_185_GM_vs_Distance_Plot_M7','png')    
    
% 5 Sec (0.2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,54), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,54), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,54), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,54), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,54), 'k-','LineWidth',1); 
    ylabel('Median 5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA5P0_185_GM_vs_Distance_Plot_M7','png')  
    
%% M7.5

% PGA
    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,5), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,5), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,5), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,5), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,5), 'k-','LineWidth',1); 
    ylabel('Median Peak Ground Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_PGA_185_GM_vs_Distance_Plot_M7p5','png')

% 0.1 Sec (10 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,10), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,10), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,10), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,10), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,10), 'k-','LineWidth',1); 
    ylabel('Median 0.1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P1_185_GM_vs_Distance_Plot_M7p5','png')    

% 0.2 Sec (5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,15), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,15), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,15), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,15), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,15), 'k-','LineWidth',1); 
    ylabel('Median 0.2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P2_185_GM_vs_Distance_Plot_M7p5','png')      
    
% 0.3 Sec (3.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,20), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,20), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,20), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,20), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,20), 'k-','LineWidth',1); 
    ylabel('Median 0.3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P3_185_GM_vs_Distance_Plot_M7p5','png')          
    
% 0.5 Sec (2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,25), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,25), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,25), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,25), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,25), 'k-','LineWidth',1); 
    ylabel('Median 0.5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P5_185_GM_vs_Distance_Plot_M7p5','png')     

% 0.75 Sec (1.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,30), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,30), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,30), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,30), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,30), 'k-','LineWidth',1); 
    ylabel('Median 0.75 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P75_185_GM_vs_Distance_Plot_M7p5','png')       
    
% 1 Sec (1 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,35), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,35), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,35), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,35), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,35), 'k-','LineWidth',1); 
    ylabel('Median 1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7p5';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA1P0_185_GM_vs_Distance_Plot_M7p5','png')      
    
% 2 Sec (0.5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,40), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,40), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,40), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,40), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,40), 'k-','LineWidth',1); 
    ylabel('Median 2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA2P0_185_GM_vs_Distance_Plot_M7p5','png')     
    
% 3 Sec (0.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,45), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,45), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,45), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,45), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,45), 'k-','LineWidth',1); 
    ylabel('Median 3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA3P0_185_GM_vs_Distance_Plot_M7p5','png')     
        
% 4 Sec (0.25 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,50), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,50), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,50), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,50), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,50), 'k-','LineWidth',1); 
    ylabel('Median 4 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA4P0_185_GM_vs_Distance_Plot_M7p5','png')    
    
% 5 Sec (0.2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,55), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,55), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,55), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,55), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,55), 'k-','LineWidth',1); 
    ylabel('Median 5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M7.5 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA5P0_185_GM_vs_Distance_Plot_M7p5','png') 
    
%% M8

% PGA
    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,6), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,6), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,6), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,6), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,6), 'k-','LineWidth',1); 
    ylabel('Median Peak Ground Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_PGA_185_GM_vs_Distance_Plot_M8','png')

% 0.1 Sec (10 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,11), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,11), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,11), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,11), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,11), 'k-','LineWidth',1); 
    ylabel('Median 0.1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P1_185_GM_vs_Distance_Plot_M8','png')    

% 0.2 Sec (5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,16), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,16), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,16), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,16), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,16), 'k-','LineWidth',1); 
    ylabel('Median 0.2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P2_185_GM_vs_Distance_Plot_M8','png')      
    
% 0.3 Sec (3.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,21), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,21), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,21), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,21), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,21), 'k-','LineWidth',1); 
    ylabel('Median 0.3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P3_185_GM_vs_Distance_Plot_M8','png')          
    
% 0.5 Sec (2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,26), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,26), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,26), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,26), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,26), 'k-','LineWidth',1); 
    ylabel('Median 0.5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P5_185_GM_vs_Distance_Plot_M8','png')     

% 0.75 Sec (1.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,31), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,31), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,31), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,31), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,31), 'k-','LineWidth',1); 
    ylabel('Median 0.75 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA0P75_185_GM_vs_Distance_Plot_M8','png')       
    
% 1 Sec (1 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,36), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,36), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,36), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,36), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,36), 'k-','LineWidth',1); 
    ylabel('Median 1 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA1P0_185_GM_vs_Distance_Plot_M8','png')      
    
% 2 Sec (0.5 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,41), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,41), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,41), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,41), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,41), 'k-','LineWidth',1); 
    ylabel('Median 2 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA2P0_185_GM_vs_Distance_Plot_M8','png')     
    
% 3 Sec (0.33 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,46), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,46), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,46), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,46), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,46), 'k-','LineWidth',1); 
    ylabel('Median 3 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA3P0_185_GM_vs_Distance_Plot_M8','png')     
        
% 4 Sec (0.25 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,51), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,51), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,51), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,51), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,51), 'k-','LineWidth',1); 
    ylabel('Median 4 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA4P0_185_GM_vs_Distance_Plot_M8','png')    
    
% 5 Sec (0.2 Hz)

    figure('visible','off')
    % GM vs. Distance Plots
    loglog(R, nshmp_haz_data_ASK_14(1:19,56), 'b-','LineWidth',1); 
    hold on
    grid on 
    loglog(R, nshmp_haz_data_BSSA_14(1:19,56), 'r-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CB_14(1:19,56), 'g-','LineWidth',1); 
    loglog(R, nshmp_haz_data_CY_14(1:19,56), 'c-','LineWidth',1); 
    %loglog(R, nshmp_haz_data_IDRISS_14(1:19,56), 'k-','LineWidth',1); 
    ylabel('Median 5 Second Spectral Acceleration (g)')
    ylim([10^-4 10])
    xlim([5 300])
    xlabel('Source Distance (km)');
    h = legend('ASK14','BSSA14','CB14','CY14');
    set(h,'Interpreter','none','FontSize',8);
    h = title('Parameters = strike-slip fault, dip = 90, zTop = 5 km, width = 10 km, zHyp = 8 km','FontSize',10); % Subtitle
    set(h,'Interpreter','none');
    set(h,'Position',[40 10]); % Position of subtilte, [center value, y-axis]
    title1 = 'NGA-W2 Median GMM Comparison - M8 - Vs30 = 185 m/s';
    text(12,20,title1,'FontWeight','bold'); % Position of Title [x-axis, y-axis]
    saveas(gcf,'NGAW2_SA5P0_185_GM_vs_Distance_Plot_M8','png') 