%% Planar Array
close all; clear; clc; % clean the work environment
addpath(genpath('.')); % add folder to MATLAB path

%% Parameters
c = 343; % speed of sound [m/s]

%constant beamwidth parameters
desired_bw_xz = 15; % desired beamwidth in the XZ plain [deg]
desired_bw_yz = 30; % desired beamwidth in the YZ plain [deg]

freqs = 0:10:8e3; % array of frequencies [Hz] to compute at
f_plot = [2e3,5e3,8e3]; % array of frequencies [Hz] to plot beampatterns at

planar_pattern = 'rectangle_star'; % Planar array geometry. Examples: 'full' or 'rectangle_star'. For more examples, see getSupportOfPA.m
do_kronecker = false; % use Kronecker product. Otherwise, uses trade-off beamformer
alpha = 0.5; % trade-off parameter

% figures
linewd = 1.5;
linewd_figure_3 = 0.8;
linewidth_markers = 0.8;
markerSize = 5;
hcfontsize = 9;

%% Design Linear Array
% The first one we design to maintain a beamwidth of desired_bw_xz
LA_A = nonIterativeAlgorithm_sensorPositions(c, freqs, 11, desired_bw_xz, 0.1e-2, 1.36, 0.034);
% The second one we design to maintain a beamwidth of desired_bw_yz
LA_B = nonIterativeAlgorithm_sensorPositions(c, freqs,  9, desired_bw_yz, 0.1e-2, 1.00, 0.034);

%% Build Planar Array
[weights_PA,Zopt_PA,DF_PA,WNG_PA,BW_xz,BW_yz,PA_x,PA_y,phi,theta] = planar_function(LA_A, LA_B, desired_bw_xz, desired_bw_yz, freqs, f_plot, c, planar_pattern, 0.01, do_kronecker);

%% plot the arrays
xlim_max = max(PA_x)*100+5;
ylim_max = max(PA_y)*100+5;
lim_max = max(xlim_max, ylim_max);

figure()
tiledlayout(1,1,'TileSpacing','tight','Padding','compact')
nexttile
hAX = gca;
plot(PA_x*100, PA_y*100, 'o','linewidth',linewidth_markers,'MarkerSize',markerSize)
title('Sensor Positions','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$x~(\mathrm{cm})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$y~(\mathrm{cm})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlim([-xlim_max,xlim_max])
ylim([-ylim_max,ylim_max])

set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(hAX, 'LineWidth', linewd);

%% plot beampattern 3D
try
    close(9)
catch
end
figure()
global title_text

T_fig_3D = tiledlayout(1, length(f_plot), 'TileSpacing','compact','Padding','loose');

clim = -20;
tile_count = 0;
for s_freq = 1:length(f_plot)
    tile_count = tile_count + 1;
    s_tile = nexttile(tile_count);

    Zopt_PA_curr = Zopt_PA(:,:,freqs==f_plot(s_freq));

    title_text = ['$f$ = ',num2str(f_plot(s_freq)/1e3),' kHz'];
    patternCustomEditted(max(db(Zopt_PA_curr).',clim), theta*180/pi, phi*180/pi,'CoordinateSystem','polar');

    colormap(s_tile, jet)
    hcb = colorbar;
    set(hcb,'visible','off')

    set(gca,'YDir','normal')
    set(gca,'TitleFontWeight','normal');
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', hcfontsize);
    set(gca, 'LineWidth', linewd_figure_3);
end

colorbar_beampattern = colorbar;
colorbar_beampattern.Title.String = '(dB)';
colorbar_beampattern.Title.FontSize = hcfontsize;
colorbar_beampattern.Title.Interpreter = 'latex';
colorbar_beampattern.Title.FontName = 'Times New Roman';
colorbar_beampattern.Layout.Tile = 'east';

%% plot BW
figure()
tiledlayout(1,1,'TileSpacing','tight','Padding','tight')

% BW
nexttile
plot(freqs/1e3, BW_xz,'LineWidth',linewd)
hold on
plot(freqs/1e3, BW_yz,'--','LineWidth',linewd)
hold off
title('Beamwidth','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$f~(\mathrm{kHz})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('$b_{\theta}~(\mathrm{deg})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylim([0 180])
legend('XZ Plain','YZ Plain','Location','best','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

%% plot WNG and Directivity
figure()
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')

nexttile
plot(freqs/1e3, 10*log10(WNG_PA),'LineWidth',linewd)
title('White Noise Gain','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('${\cal W}~(\mathrm{dB})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
yl = ylim;
ylim([0 yl(2)])
legend(['${\cal W}_{[',num2str(freqs(1)/1e3),',~',num2str(freqs(end)/1e3),'~\mathrm{kHz}]} = ',num2str(calcWidebandDI(WNG_PA),'%01.1f'),'~\mathrm{dB}$'],'Location','best','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);

% plot Directivity
nexttile
plot(freqs/1e3, 10*log10(DF_PA),'LineWidth',linewd)
title('Narrowband Directivity Factor','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
xlabel('$f~(\mathrm{kHz})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
ylabel('${\cal D}~(\mathrm{dB})$','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
yl = ylim;
ylim([0 yl(2)])
legend(['${\cal DI}_{[',num2str(freqs(1)/1e3),',~',num2str(freqs(end)/1e3),'~\mathrm{kHz}]} = ',num2str(calcWidebandDI(DF_PA),'%01.1f'),'~\mathrm{dB}$'],'Location','best','FontSize',hcfontsize,'Interpreter','latex','FontName','Times New Roman')
set(gca,'TitleFontWeight','normal');
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd);
