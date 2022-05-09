%% Linear Array
close all; clear; clc; % clean the work environment
addpath(genpath('.')); % add folder to MATLAB path

%% Parameters
c = 343; % speed of sound [m/s]
freqs = 0:10:8e3; % array of frequencies [Hz] to compute at
desired_bw = 15; % desired beamwidth [deg]
beta_res = 1e-3; % resolution of Kaiser window shape factor. Example: 1e-3
L_support_option = 'supports'; % choose what the active sensors are. Options: 'single', 'supports', 'custom'
use_continuous_kaiser = true; % if to sample the continuous Kaiser window. Otherwise, uses the discrete Kaiser window
use_trapezoidal_integration = true; % if to use the trapezoidal integration technique

% figures
linewd = 1.5;
linewd_figure_3 = 0.8;
linewidth_markers = 0.8;
markerSize = 5;
hcfontsize = 9;

%% Get Linear Array coordinates
x_LA = [0.038    0.079    0.143    0.292    0.748]; % kaiser_mod
x_LA = [-x_LA(end:-1:1), 0, x_LA]';

% x_LA = nonIterativeAlgorithm_sensorPositions(c, freqs, 11, desired_bw, 0.1e-2, 1.36, 0.034);

%% constant beamwidth calculation
[weights_LA, Zopt_LA, directivityFactor_LA, WNG_LA, BW_LA, phi, fmin_LA, beta_opt_LA, L_opt_LA] = Algorithm1_AttainingTheWeights...
    (c, freqs, x_LA, desired_bw, beta_res, L_support_option, use_continuous_kaiser, use_trapezoidal_integration);

%% plot the arrays
M = length(x_LA);
N = (M-1)/2; % amount of elements that have positive x coordinates
x_lim_max = max(x_LA)*100+5;
y_diff = 1;
figure()
tiledlayout(2,1,'TileSpacing','tight','Padding','compact')
nexttile
hAX = gca;
plot(x_LA*100, zeros(size(x_LA)), 'o','linewidth',linewidth_markers,'MarkerSize',markerSize)
title('Sensor Positions','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$x~(\mathrm{cm})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlim([-x_lim_max,x_lim_max])
set(gca,'ytick',[])
ylim([-0.9, 0.9])
daspect([6 1 1])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(hAX, 'LineWidth', linewd);

%% plot weights
figure()
T_fig3 = tiledlayout(1,2,'TileSpacing','compact','Padding','tight');
T_weights = tiledlayout(T_fig3,1,1,'TileSpacing','tight','Padding','tight');
T_weights.Layout.Tile = 1;
T_weights.Layout.TileSpan = [1 1];

clim_m = max(weights_LA,[],'all');
w1 = nexttile(T_weights);
clim_weights = [0 clim_m];
imagesc(0:size(weights_LA,1)-1,freqs/1e3,weights_LA',clim_weights)
title('Weights','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$l~(\mathrm{index})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$f~(\mathrm{kHz})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
colormap(w1,jet);

set(gca,'YDir','normal')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd_figure_3);
set(gca,'XTick',0:N);

colorbar_weights = colorbar;
colorbar_weights.Layout.Tile = 'east';

% plot beampattern
T_pattern = tiledlayout(T_fig3,1,1,'TileSpacing','tight','Padding','tight');
T_pattern.Layout.Tile = 2;
T_pattern.Layout.TileSpan = [1 1];

b1 = nexttile(T_pattern);
clim_beampatter = [-30 0];
imagesc(phi*180/pi,freqs/1e3,db(Zopt_LA)',clim_beampatter)
title('Beampattern','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$\varphi~(\mathrm{deg})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
colormap(b1,jet)

set(gca,'YDir','normal')
set(gca,'TitleFontWeight','normal');
% add a contour
hold on
contour(phi*180/pi, freqs/1e3, db(Zopt_LA)', [db(sqrt(0.5)) db(sqrt(0.5))], '--k','LineWidth',linewd_figure_3); % Add 3dB contour line
hold off
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd_figure_3);
set(gca,'XTick',0:45:180);

colorbar_beampattern = colorbar;
colorbar_beampattern.Title.String = '(dB)';
colorbar_beampattern.Title.FontSize = hcfontsize;
colorbar_beampattern.Title.Interpreter = 'latex';
colorbar_beampattern.Title.FontName = 'Times New Roman';
colorbar_beampattern.Layout.Tile = 'east';

%% plot performance measures
% plot BW vs. f
figure()
tiledlayout(3,1,'TileSpacing','tight','Padding','tight')
nexttile
plot(freqs/1e3, BW_LA,'LineWidth',linewd)
title('Beamwidth','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$b_{\varphi}~(\mathrm{deg})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylim([0 180])
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

% plot WNG
nexttile
plot(freqs/1e3, 10*log10(WNG_LA),'LineWidth',linewd)
title('White Noise Gain','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('${\cal W}~(\mathrm{dB})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
yl = ylim;
ylim([0 yl(2)])
legend(['${\cal W}_{[',num2str(freqs(1)/1e3),',~',num2str(freqs(end)/1e3),'~\mathrm{kHz}]} = ',num2str(calcWidebandDI(WNG_LA),'%01.1f'),'~\mathrm{dB}$'],'Location','best','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

% plot Directivity
nexttile
plot(freqs/1e3, 10*log10(directivityFactor_LA),'LineWidth',linewd)
title('Narrowband Directivity Factor','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$f~(\mathrm{kHz})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('${\cal D}~(\mathrm{dB})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
legend(['${\cal DI}_{[',num2str(freqs(1)/1e3),',~',num2str(freqs(end)/1e3),'~\mathrm{kHz}]} = ',num2str(calcWidebandDI(directivityFactor_LA),'%01.1f'),'~\mathrm{dB}$'],'Location','best','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

%% plot the window parameters
figure()
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')
nexttile
plot(freqs/1e3, beta_opt_LA,'LineWidth',linewd)
title('Kaiser Window Shape Factor - $\beta$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$\beta$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylim([0 10])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

nexttile
plot(freqs/1e3, L_opt_LA,'LineWidth',linewd)

title('Index of the Window Support','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$f~(\mathrm{kHz})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$i~(\mathrm{index})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylim([1-0.5 N+0.5])
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
