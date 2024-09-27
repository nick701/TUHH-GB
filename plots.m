%% Publications plots
clear all;clc;close all

NRMSE_Seq1 = load("NRMSE_Seq1_s1_100_win5000_r5.mat");
NRMSE_Seq2 = load("NRMSE_Seq2_s1_100_win5000_r5.mat");

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "NRMS_AUG_s1_100_win_5000.eps";
fullFilePath = fullfile(folderPath, filename);

figure;
lin_width = 1.2;
lin_width_tick = 1;
font_size = 8;

t = tiledlayout("vertical");
t.TileSpacing = 'compact';
t.Padding = 'compact';



nexttile
hold on
plot(NRMSE_Seq1.s_save,NRMSE_Seq1.G1_ERR,'LineWidth',lin_width,'Color','#66c2a5','Marker','diamond','MarkerFaceColor','#66c2a5')
plot(NRMSE_Seq1.s_save,NRMSE_Seq1.G2_ERR,'LineWidth',lin_width,'Color','#fc8d62','Marker','diamond','MarkerFaceColor','#fc8d62')
plot(NRMSE_Seq1.s_save,NRMSE_Seq1.G3_ERR,'LineWidth',lin_width,'Color','#8da0cb','Marker','diamond','MarkerFaceColor','#8da0cb')
plot(NRMSE_Seq1.s_save,NRMSE_Seq1.G4_ERR,'LineWidth',lin_width,'Color','#e78ac3','Marker','diamond','MarkerFaceColor','#e78ac3')
plot(NRMSE_Seq1.s_save,NRMSE_Seq1.G5_ERR,'LineWidth',lin_width,'Color','#a6d854','Marker','diamond','MarkerFaceColor','#a6d854')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
title('EA','Interpreter','latex')
ylim([0 0.5])
xlim([0 101])
grid minor
hold off
box on

nexttile
hold on
plot(NRMSE_Seq1.s_save,NRMSE_Seq2.G1_ERR,'LineWidth',lin_width,'Color','#66c2a5','Marker','diamond','MarkerFaceColor','#66c2a5')
plot(NRMSE_Seq1.s_save,NRMSE_Seq2.G2_ERR,'LineWidth',lin_width,'Color','#fc8d62','Marker','diamond','MarkerFaceColor','#fc8d62')
plot(NRMSE_Seq1.s_save,NRMSE_Seq2.G3_ERR,'LineWidth',lin_width,'Color','#8da0cb','Marker','diamond','MarkerFaceColor','#8da0cb')
plot(NRMSE_Seq1.s_save,NRMSE_Seq2.G4_ERR,'LineWidth',lin_width,'Color','#e78ac3','Marker','diamond','MarkerFaceColor','#e78ac3')
plot(NRMSE_Seq1.s_save,NRMSE_Seq2.G5_ERR,'LineWidth',lin_width,'Color','#a6d854','Marker','diamond','MarkerFaceColor','#a6d854')
ylim([0 0.5])
xlim([0 101])
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
title('EB','Interpreter','latex')
grid minor
hold off
box on


leg = legend('$G_1$','$G_2$','$G_3$','$G_4$','$G_5$','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Interpreter = 'latex';
leg.FontSize = 10;

ylabel(t, 'NRMSE', 'Interpreter', 'latex','FontSize',10);
xlabel(t, 'No. delayed coordinates', 'Interpreter', 'latex','FontSize',10);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 8.152])

%print(gcf,fullFilePath,'-depsc','-vector');
exportgraphics(gcf,'NRMS_AUG_s1_100_win_5000.png','Resolution',300)
%%
clear all; clc; close all;

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "DMD-AUGDMD-Seq1.eps";
fullFilePath = fullfile(folderPath, filename);

F_DMD_AUG_Seq1 = load("f_dmd_aug_win500_r5_seq1_predictionWin250.mat");
%F_DMD_AUG_Seq2 = load("f_dmd_aug_win500_r5_s50_seq2.mat");
%F_DMD_AUG_Seq3 = load("f_dmd_aug_win500_r5_s50_seq3.mat");

F_DMD_STAN_Seq1 = load('f_dmd_standard_win500_r5_seq1_predictionWin250.mat');
%F_DMD_STAN_Seq2 = load('f_dmd_standard_win500_r5_seq2.mat');
%F_DMD_STAN_Seq3 = load('f_dmd_standard_win500_r5_seq3.mat');
%time = load('time.mat');

f = figure;
lin_width = 1.2;
lin_width_tick = 1;
font_size = 8;

t = tiledlayout("vertical");
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
hold on
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_DMD_1*10,'Color',"#984ea3",'LineWidth',lin_width);
plot(F_DMD_AUG_Seq1.t,F_DMD_AUG_Seq1.F_DMD_1*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_AUG_1*10,'Color',"#4daf4a",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
    'Label','$\rightarrow\hat{x}_{i}$','Interpreter','latex','LabelOrientation','horizontal',...
    'LabelVerticalAlignment','bottom');
hold off
box on
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
%xlim([0 max(time.t)])
title('$\mathbf{G}_1$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
%grid minor

nexttile
hold on

plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_DMD_2*10,'Color',"#984ea3",'LineWidth',lin_width);
plot(F_DMD_AUG_Seq1.t,F_DMD_AUG_Seq1.F_DMD_2*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_AUG_2*10,'Color',"#4daf4a",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
    'Label','$\rightarrow\hat{x}_{i}$','Interpreter','latex','LabelOrientation','horizontal',...
    'LabelVerticalAlignment','bottom');
hold off
box on
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
%xlim([0 max(time.t)])
title('$\mathbf{G}_2$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size);
grid on
%grid minor 

nexttile
hold on

plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_DMD_3*10,'Color',"#984ea3",'LineWidth',lin_width);
plot(F_DMD_AUG_Seq1.t,F_DMD_AUG_Seq1.F_DMD_3*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_AUG_3*10,'Color',"#4daf4a",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
    'Label','$\rightarrow\hat{x}_{i}$','Interpreter','latex','LabelOrientation','horizontal',...
    'LabelVerticalAlignment','bottom');
hold off
box on
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
%xlim([0 max(F_DMD_AUG_Seq1.t)])
title('$\mathbf{G}_3$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
%grid minor

nexttile
hold on

plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_DMD_4*10,'Color',"#984ea3",'LineWidth',lin_width);
plot(F_DMD_AUG_Seq1.t,F_DMD_AUG_Seq1.F_DMD_4*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_AUG_4*10,'Color',"#4daf4a",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
    'Label','$\rightarrow\hat{x}_{i}$','Interpreter','latex','LabelOrientation','horizontal',...
    'LabelVerticalAlignment','bottom');
hold off
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
%xlim([0 max(F_DMD_AUG_Seq1.t)])
title('$\mathbf{G}_4$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
yticks([-0.2 0 0.2])
grid on
%grid minor
box on
nexttile
hold on

plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_DMD_5*10,'Color',"#984ea3",'LineWidth',lin_width);
plot(F_DMD_AUG_Seq1.t,F_DMD_AUG_Seq1.F_DMD_5*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(F_DMD_STAN_Seq1.t,F_DMD_STAN_Seq1.F_AUG_5*10,'Color',"#4daf4a",'LineWidth',lin_width);
%xlim([0 max(time.t)])
xline(F_DMD_STAN_Seq1.t(end)-F_DMD_STAN_Seq1.extraTime,'LineWidth',lin_width,...
    'Label','$\rightarrow\hat{x}_{i}$','Interpreter','latex','LabelOrientation','horizontal',...
    'LabelVerticalAlignment','bottom');
hold off
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex','Margin',5)
%xlabel('Time [s]','Interpreter','latex')
title('$\mathbf{G}_5$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
%yticks([-0.002 0 0.002])
box on
grid on

leg = legend('DMD','AugDMD','Measured','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Interpreter = 'latex';
leg.FontSize = 10;
leg.Direction = "reverse";
ylabel(t, '$\dot{z}$(t) [mm/s]', 'Interpreter', 'latex','FontSize',10);
xlabel(t, 'Time [s]', 'Interpreter', 'latex','FontSize',10);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 12.88])

%print(gcf,fullFilePath,'-depsc','-vector');
exportgraphics(gcf,'DMD_AUGDMD_SEQ1.png','Resolution',300)

%%

%%
clear all; clc; close all;

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "NRMSE_prediction.eps";
fullFilePath = fullfile(folderPath, filename);

G_ERR_Seq1 = load("NRMSE_Seq1_prediction4.mat");
G_ERR_Seq2 = load("NRMSE_Seq2_prediction4.mat");


figure;
lin_width = 1.2;
lin_width_tick = 1;
font_size = 8;

t = tiledlayout("vertical");
t.TileSpacing = 'compact';
t.Padding = 'compact';
nexttile
hold on
plot(G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq1.G1_ERR,'LineWidth',lin_width,'Color','#66c2a5')
plot(G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq1.G2_ERR,'LineWidth',lin_width,'Color','#fc8d62')
plot(G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq1.G3_ERR,'LineWidth',lin_width,'Color','#8da0cb')
plot(G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq1.G4_ERR,'LineWidth',lin_width,'Color','#e78ac3')
plot(G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq1.G5_ERR,'LineWidth',lin_width,'Color','#a6d854')
hold off
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
ylim([0 0.3])
xlim([0 max(G_ERR_Seq1.t)])
%xlim([G_ERR_Seq1.t(G_ERR_Seq1.t_extra_0+1) G_ERR_Seq1.t(end)])
title('EA','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
grid minor
box on
nexttile
hold on
plot(G_ERR_Seq2.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq2.G1_ERR,'LineWidth',lin_width,'Color','#66c2a5')
plot(G_ERR_Seq2.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq2.G2_ERR,'LineWidth',lin_width,'Color','#fc8d62')
plot(G_ERR_Seq2.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq2.G3_ERR,'LineWidth',lin_width,'Color','#8da0cb')
plot(G_ERR_Seq2.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq2.G4_ERR,'LineWidth',lin_width,'Color','#e78ac3')
plot(G_ERR_Seq2.t(G_ERR_Seq1.t_extra_0+1:end),G_ERR_Seq2.G5_ERR,'LineWidth',lin_width,'Color','#a6d854')
hold off
ylim([0 0.3])
xlim([0 max(G_ERR_Seq2.t)])
%ylabel('$\dot{z}$(t) [m/s$^2$]','Interpreter','latex')
%xlabel('Time [s]','Interpreter','latex')
title('EB','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
grid minor
box on
leg = legend('$G_1$','$G_2$','$G_3$','$G_4$','$G_5$','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Interpreter = 'latex';
leg.FontSize = 10;

ylabel(t, 'NRMSE', 'Interpreter', 'latex','FontSize',10);
xlabel(t, 'Time [s]', 'Interpreter', 'latex','FontSize',10);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 8.152])


%print(gcf,fullFilePath,'-depsc','-vector');
exportgraphics(gcf,'NRMSE_prediction.png','Resolution',500)
%%
clear all; clc; close all;

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "DMD-AUGDMD-Seq2_fullWindow.eps";
fullFilePath = fullfile(folderPath, filename);

F_DMD_STAN_Seq2 = load("f_dmd_standard_win_500_r5_seq2_prediction.mat");
F_DMD_AUG_Seq2 = load("f_dmd_aug_win_500_r5_seq2_prediction.mat");


%F_DMD_STAN_Seq2 = load("f_dmd_standard_win_full_r5_seq2_predictionWin250.mat");
%F_DMD_AUG_Seq2 = load("f_dmd_aug_win_full_r5_seq2_predictionWin250.mat");


% F_DMD_STAN_Seq1 = load('f_dmd_standard_win500_r5_seq1.mat');
% F_DMD_STAN_Seq2 = load('f_dmd_standard_win500_r5_seq2.mat');

time = F_DMD_STAN_Seq2.t;

figure;
lin_width = 0.8;
lin_width_tick = 1;
font_size = 8;
label_prediction_offset = 1.5;
t = tiledlayout("vertical");
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
hold on
plot(time,F_DMD_STAN_Seq2.F_AUG_1*10,'Color',"#4daf4a",'LineWidth',lin_width);
plot(time,F_DMD_AUG_Seq2.F_DMD_1*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(time,F_DMD_STAN_Seq2.F_DMD_1*10,'Color',"#984ea3",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,'LineWidth',lin_width*1.5);
text(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,label_prediction_offset*max(F_DMD_STAN_Seq2.F_AUG_1*10),...
    '$\rightarrow\hat{x}_{i}$','Interpreter','latex')
hold off
xlim([0 max(time)])
ylim([-abs(max(F_DMD_AUG_Seq2.F_DMD_1*10)) abs(max(F_DMD_AUG_Seq2.F_DMD_1*10))])
title('$\mathbf{G}_1$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
box on

nexttile
hold on
plot(time,F_DMD_STAN_Seq2.F_AUG_2*10,'Color',"#4daf4a",'LineWidth',lin_width);
plot(time,F_DMD_AUG_Seq2.F_DMD_2*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(time,F_DMD_STAN_Seq2.F_DMD_2*10,'Color',"#984ea3",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,'LineWidth',lin_width*1.5);
text(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,label_prediction_offset*max(F_DMD_STAN_Seq2.F_AUG_2*10),...
    '$\rightarrow\hat{x}_{i}$','Interpreter','latex')
hold off
xlim([0 max(time)])
ylim([-abs(max(F_DMD_AUG_Seq2.F_DMD_2*10)) abs(max(F_DMD_AUG_Seq2.F_DMD_2*10))])
title('$\mathbf{G}_2$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on

box on  
nexttile
hold on
plot(time,F_DMD_STAN_Seq2.F_AUG_3*10,'Color',"#4daf4a",'LineWidth',lin_width);
plot(time,F_DMD_AUG_Seq2.F_DMD_3*10,'Color',"#ff7f00",'LineWidth',lin_width);

plot(time,F_DMD_STAN_Seq2.F_DMD_3*10,'Color',"#984ea3",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,'LineWidth',lin_width*1.5);
text(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,label_prediction_offset*max(F_DMD_STAN_Seq2.F_AUG_3*10),...
    '$\rightarrow\hat{x}_{i}$','Interpreter','latex')
hold off
xlim([0 max(time)])
ylim([-abs(max(F_DMD_AUG_Seq2.F_DMD_3*10)) abs(max(F_DMD_AUG_Seq2.F_DMD_3*10))])
%ylim([-1.5 1.5])
title('$\mathbf{G}_3$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
grid on
box on

nexttile
hold on
plot(time,F_DMD_STAN_Seq2.F_AUG_4*10,'Color',"#4daf4a",'LineWidth',lin_width);
plot(time,F_DMD_AUG_Seq2.F_DMD_4*10,'Color',"#ff7f00",'LineWidth',lin_width);
plot(time,F_DMD_STAN_Seq2.F_DMD_4*10,'Color',"#984ea3",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,'LineWidth',lin_width*1.5);
text(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,label_prediction_offset*max(F_DMD_STAN_Seq2.F_AUG_4*10),...
    '$\rightarrow\hat{x}_{i}$','Interpreter','latex')
hold off
xlim([0 max(time)])
ylim([-abs(max(F_DMD_AUG_Seq2.F_DMD_4*10)) abs(max(F_DMD_AUG_Seq2.F_DMD_4*10))])
%ylim([-0.25 0.25])
title('$\mathbf{G}_4$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
%yticks([-0.02 0 0.02])
grid on
box on

nexttile
hold on
plot(time,F_DMD_STAN_Seq2.F_AUG_5*10,'Color',"#4daf4a",'LineWidth',lin_width);
plot(time,F_DMD_AUG_Seq2.F_DMD_5*10,'Color',"#ff7f00",'LineWidth',lin_width);

plot(time,F_DMD_STAN_Seq2.F_DMD_5*10,'Color',"#984ea3",'LineWidth',lin_width);
xline(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,'LineWidth',lin_width*1.5);
text(F_DMD_STAN_Seq2.t(end)-F_DMD_STAN_Seq2.extraTime,label_prediction_offset*max(F_DMD_STAN_Seq2.F_AUG_5*10),...
    '$\rightarrow\hat{x}_{i}$','Interpreter','latex')
xlim([0 max(time)])
ylim([-0.025 0.025])
hold off
title('$\mathbf{G}_5$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
%yticks([-0.002 0 0.002])
grid on
box on

leg = legend('Measured','AugDMD','DMD','Orientation', 'Horizontal');
leg.Layout.Tile = 'north';
leg.Interpreter = 'latex';
leg.FontSize = 10;
ylabel(t, '$\dot{z}$(t) [mm/s]', 'Interpreter', 'latex','FontSize',10);
xlabel(t, 'Time [s]', 'Interpreter', 'latex','FontSize',10);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 45, 20.88])

%print(gcf,fullFilePath,'-depsc','-vector');

exportgraphics(gcf,'DMD-AUGDMD-Seq2_fullWindow.png','Resolution',300)

%% FFT
clear all
clc
close all

folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
filename = "FF1_Seq1.eps";
fullFilePath = fullfile(folderPath, filename);

name = "Seq1";
load(name+".mat")

n = 100;
G1 = downsample(G1V3,n);
G2 = downsample(G2V3,n);
G3 = downsample(G3V3,n);
G4 = downsample(G4V3,n);
G5 = downsample(G5V3,n);
time = downsample(Time,n);

% G1 = G1V3;
% G2 = G2V3;
% G3 = G3V3;
% G4 = G4V3;
% G5 = G5V3;
% time = Time;
lin_width = 1.5;
lin_width_tick = 1;
font_size = 8;
% we calculate the delta of time
T = time(2)-time(1);
Fs = 1/T;
L = length(G1);
Y_G1 = fft(G1);
Y_G2 = fft(G2);
Y_G3 = fft(G3);
Y_G4 = fft(G4);
Y_G5 = fft(G5);
grid on
hold on
P2 = abs(Y_G1/L);
P1_G1 = P2(1:L/2+1);
P1_G1(2:end-1) = 2*P1_G1(2:end-1);
f = Fs/L*(0:(L/2));
plot(f,P1_G1,"LineWidth",lin_width,'Color',"#66c2a5")
%yline(max(P1_G1),'Color',"#0072BD",'LineStyle','--')
[max_P1_G1,indx_P1_G1] = max(P1_G1);
plot(f(indx_P1_G1),max_P1_G1,'Marker','o','MarkerEdgeColor',"#66c2a5",'MarkerFaceColor',"#66c2a5")

P2 = abs(Y_G2/L);
P1_G2 = P2(1:L/2+1);
P1_G2(2:end-1) = 2*P1_G2(2:end-1);
f = Fs/L*(0:(L/2));
plot(f,P1_G2,"LineWidth",lin_width,'Color',"#fc8d62")
%yline(max(P1_G2),'Color',"#D95319",'LineStyle','--')
[max_P1_G2,indx_P1_G2] = max(P1_G2);
plot(f(indx_P1_G2),max_P1_G2,'Marker','o','MarkerEdgeColor',"#fc8d62",'MarkerFaceColor',"#fc8d62")

P2 = abs(Y_G3/L);
P1_G3 = P2(1:L/2+1);
P1_G3(2:end-1) = 2*P1_G3(2:end-1);
f = Fs/L*(0:(L/2));
plot(f,P1_G3,"LineWidth",lin_width,'Color',"#8da0cb")
%yline(max(P1_G3),'Color',"#EDB120",'LineStyle','--')
[max_P1_G3,indx_P1_G3] = max(P1_G3);
plot(f(indx_P1_G3),max_P1_G3,'Marker','o','MarkerEdgeColor',"#8da0cb",'MarkerFaceColor',"#8da0cb")

P2 = abs(Y_G4/L);
P1_G4 = P2(1:L/2+1);
P1_G4(2:end-1) = 2*P1_G4(2:end-1);
f = Fs/L*(0:(L/2));
plot(f,P1_G4,"LineWidth",lin_width,'Color',"#e78ac3")
%yline(max(P1_G4),'Color',"#7E2F8E",'LineStyle','--')
[max_P1_G4,indx_P1_G4] = max(P1_G4);
plot(f(indx_P1_G4),max_P1_G4,'Marker','o','MarkerEdgeColor',"#e78ac3",'MarkerFaceColor',"#e78ac3")

P2 = abs(Y_G5/L);
P1_G5 = P2(1:L/2+1);
P1_G5(2:end-1) = 2*P1_G5(2:end-1);
f = Fs/L*(0:(L/2));
plot(f,P1_G5,"LineWidth",lin_width,'Color',"#a6d854")
%yline(max(P1_G5),'Color',"#77AC30",'LineStyle','--')
[max_P1_G1,indx_P1_G5] = max(P1_G5);
plot(f(indx_P1_G5),max_P1_G1,'Marker','o','MarkerEdgeColor',"#a6d854",'MarkerFaceColor',"#a6d854")
hold off


% create a new pair of axes inside current figure

ylim([0 0.8])
text(30, -0.06,'39.18','Interpreter','latex','FontSize',font_size,'FontWeight','bold')
set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
%set(gca,'xtick', 39.18, 'xticklabel', {'39.18 Hz'})

leg = legend('$G_1$','','$G_2$','','$G_3$','','$G_4$','', '$G_5$','','Orientation','horizontal');
%leg.Layout.Tile = 'north';
leg.Interpreter = 'latex';
leg.FontSize = 10;
leg.Location = "north";
ylabel('$|G_n(f)|$', 'Interpreter', 'latex','FontSize',10);
xlabel('f [Hz]', 'Interpreter', 'latex','FontSize',10);
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 20, 8.152/2])
% 
% axes('position',[.55 .255 .25 .55])
% box on % put box around new pair of axes
% indexOfInterest = (f < 200) & (f > 0); % range of t near perturbation
% hold on
% plot(f(indexOfInterest),P1_G1(indexOfInterest),'LineWidth',lin_width) % plot on new axes
% plot(f(indexOfInterest),P1_G2(indexOfInterest),'LineWidth',lin_width) % plot on new axes
% plot(f(indexOfInterest),P1_G3(indexOfInterest),'LineWidth',lin_width) % plot on new axes
% plot(f(indexOfInterest),P1_G4(indexOfInterest),'LineWidth',lin_width) % plot on new axes
% plot(f(indexOfInterest),P1_G5(indexOfInterest),'LineWidth',lin_width) % plot on new axes
% set(gca,'TickLabelInterpreter','latex','LineWidth',lin_width_tick,'FontSize',font_size)
% 
% ylabel('$|G(f)|$', 'Interpreter', 'latex','FontSize',10);
% xlabel('f [Hz]', 'Interpreter', 'latex','FontSize',10);
% grid on

print(gcf,fullFilePath,'-depsc','-vector');

%axis tight
