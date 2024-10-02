%% DMD

clear all
close all
tic

%%
% Load the raw data from the installation of the 3rd pile. For the 5 first
% geophones. Only vertical velocity is considered.
%
% Seq1 -> 1.55 - 2.00 min; penetration 0.9 - 1.0 m (approx)
% Seq2 -> 4.40 - 4.45 min; penetration 1.75 - 2.0 m (approx)
% Seq3 -> 10.55 - 11.0 min; penetration 3.50 - 3.60 (approx)

name = "Seq1";
load(name+".mat")

% resample the data to a maximum of n = 100 (after that you have alliasing)
% (neglecting the V3)

n = 100;
G1 = downsample(G1V3,n);
G2 = downsample(G2V3,n);
G3 = downsample(G3V3,n);
G4 = downsample(G4V3,n);
G5 = downsample(G5V3,n);
time = downsample(Time,n);

% we calculate the delta of time
dt = time(2)-time(1);

% Arrange all of our data in a matrix
f_all = [G1'; G2'; G3'; G4'; G5']';
% -------
% Number of Geophones to analyze
num_geo=5;

% --------
% Control what to plot
% live plot of geophones
plot_geophones = 1;

% plot the NRMSE
plot_NRMSE = 1;
plot_NRMSE_prediction = 1;
% -------
% Control what to print (save to file)
% print the NRMSE
print_NRMSE = 0;
print_prediction_error = 0;

% Folder to save the plots
folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
% ------
% DMD Parameters
% rank truncation
r = 5;

% number of samples considered for the DMD
window = 500;%length(G1)-1;      

% total amount of rows of the aug matric
l = ceil(window/2);     

% num of steps to move the window
nstep = 1;             

% number of delay coordinates. s=1 corresponds to standard DMD
s = 50;               

% size of the step. If step < window it will be an overlap in the data
step = window;

% time considered for the extrapolation (prediction)
extraTime = 0;

% index to indicate the begining of the analysis
initstep = 1;%1;

% creation of augmented matrix (to reduce computational time
f_aug = zeros(l,s*num_geo);

% time vector including extrapolation time
t = 0:dt:(l-1)*dt+extraTime;

% data matrices accoring to DMD
X1 = zeros(s*num_geo,l-1);
X2 = zeros(s*num_geo,l-1);

% length of the extra time
t_extra_0 = length(0:dt:(l-1)*dt);

% index to create store the error data
jj = 0;

% ------
% analysis 
% standard DMD
standard_DMD = 0;

prediction_error = 0;
% augmented DMD
aug_DMD = 0;
s_aug = 50;

% number of delay coordinates
s_analysis = 1; % YOU NEED TO SET THIS TO 1 TO PERFORM VARIATON OF S
s_0 = 1; 
s_f = 101;
s_step = 10;

% RMS 
G1_ERR = [];
G2_ERR = [];
G3_ERR = [];
G4_ERR = [];
G5_ERR = [];

if (standard_DMD == 1 && aug_DMD == 1)
    error('need to set either standard DMD or Aug DMD')
end
if aug_DMD == 1
    s_0 = s;
    s_f = s;
    s_step = 1;
end

if standard_DMD == 1
    s_0 = 1;
    s_f = 1;
    s_step = 1;
    disp('Using a constant delay coordinate (or standard DMD)')
end

if (s_analysis == 1 && aug_DMD ==1)
    error('need to set either s analysis or AugDMD')
end

for s = s_0:s_step:s_f
    for j = initstep:step:nstep*step+initstep-1
    
        jj = jj+1;
        for i = 1:s
            f = f_all(j:j+window,:);
            f_aug(:,[i s+i 2*s+i 3*s+i 4*s+i]) = [f(i:l+i-1,1) f(i:l+i-1,2) f(i:l+i-1,3) f(i:l+i-1,4) f(i:l+i-1,5)];
            %f_aug = f;
        end

        X1 = f_aug(1:l-1,:)';
        X2 = f_aug(2:l,:)';
        % Perform Singular Value Decomposition (SVD) on X1
        [U, S, V] = svd(X1, "econ");
        % Truncate the SVD to the desired rank 'r'
        Ur = U(:, 1:r);
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
        % Build the low-rank approximation of the state transition matrix Atilde
        Atilde = Ur' * X2 * Vr *Sr^(-1);
        % Compute the eigenvalues (D) and eigenvectors (W) of Atilde
        [W, D] = eig(Atilde);
        % Calculate the DMD modes (Phi) using X2, Vr, Sr, and W
        Phi_aug = X2 * Vr / Sr * W;
        % Calculate the continuous-time frequencies (omega) from the eigenvalues
        lambda = diag(D);
        omega_aug = log(lambda) / dt;
        % Determine the initial condition (b) in DMD space
        b_aug = Phi_aug \ X1(:, 1);
        f_dmd_extra = Phi_aug*(b_aug.*exp(omega_aug*t));
    
        if plot_geophones == 1
            subplot(3,2,1)
            plot(t,real(f_dmd_extra(1,:)),'--')
            hold on
            plot(t,f_all(j:length(f_dmd_extra(1,:))+j-1,1))
            xline(l*dt,'--',{'Predicted','Signal'});
            legend('augDMD','Measurement')
            title('Geophone 1')
            ylabel('Velocity [m/s]')
        
            subplot(3,2,2)
            plot(t,real(f_dmd_extra(s+1,:)),'--')
            hold on
            plot(t,f_all(j:length(f_dmd_extra(1,:))+j-1,2))
            xline(l*dt,'--',{'Predicted','Signal'});
                legend('augDMD','Measurement')
            title('Geophone 2')
        
            subplot(3,2,3)
            plot(t,real(f_dmd_extra(2*s+1,:)),'--')
            hold on
            plot(t,f_all(j:length(f_dmd_extra(1,:))+j-1,3))
            xline(l*dt,'--',{'Predicted','Signal'});
                legend('augDMD','Measurement')
            title('Geophone 3')
            ylabel('Velocity [m/s]')
        
            subplot(3,2,4)
            plot(t,real(f_dmd_extra(3*s+1,:)),'--')
            hold on
            plot(t,f_all(j:length(f_dmd_extra(1,:))+j-1,4))
            xline(l*dt,'--',{'Predicted','Signal'});
                legend('augDMD','Measurement')
            title('Geophone 4')
        
            subplot(3,2,5)
            plot(t,real(f_dmd_extra(4*s+1,:)),'--')
            hold on
            plot(t,f_all(j:length(f_dmd_extra(1,:))+j-1,5))
            xline(l*dt,'--',{'Predicted','Signal'});
                legend('augDMD','Measurement')
            title('Geophone 5')
              xlabel('Time [s]')
              ylabel('Velocity [m/s]')
        
            subplot(3,2,6)
           
            plot(t,abs(real(f_dmd_extra(0*s+1,:))'-f_all(j:length(f_dmd_extra(1,:))+j-1,1)),'r--')
            hold on
            plot(t,abs(real(f_dmd_extra(1*s+1,:))'-f_all(j:length(f_dmd_extra(1,:))+j-1,2)),'b--')
            hold on
            plot(t,abs(real(f_dmd_extra(2*s+1,:))'-f_all(j:length(f_dmd_extra(1,:))+j-1,3)),'g--')
            hold on
            plot(t,abs(real(f_dmd_extra(3*s+1,:))'-f_all(j:length(f_dmd_extra(1,:))+j-1,4)),'y--')
            hold on
            plot(t,abs(real(f_dmd_extra(4*s+1,:))'-f_all(j:length(f_dmd_extra(1,:))+j-1,5)),'b--')
            xline(l*dt,'--',{'Predicted','Error'});
            ylim([0 1])
            legend('G1','G2','G3','G4','G5')
            title('Error')
            xlabel('Time [s]')
            drawnow
            clf
        end
            F_DMD_1 = real(f_dmd_extra(0*s+1,:)');
            F_DMD_2 = real(f_dmd_extra(1*s+1,:)');
            F_DMD_3 = real(f_dmd_extra(2*s+1,:)');
            F_DMD_4 = real(f_dmd_extra(3*s+1,:)');
            F_DMD_5 = real(f_dmd_extra(4*s+1,:)');
    
            F_AUG_1 = f_all(j:length(f_dmd_extra(1,:))+j-1,1);
            F_AUG_2 = f_all(j:length(f_dmd_extra(2,:))+j-1,2);
            F_AUG_3 = f_all(j:length(f_dmd_extra(3,:))+j-1,3);
            F_AUG_4 = f_all(j:length(f_dmd_extra(4,:))+j-1,4);
            F_AUG_5 = f_all(j:length(f_dmd_extra(5,:))+j-1,5);
    
          
            ERR_1_c = 0;
            ERR_2_c = 0;
            ERR_3_c = 0;
            ERR_4_c = 0;
            ERR_5_c = 0;

            if prediction_error == 1
                for i = 2:length(F_AUG_1(t_extra_0:end))
                    % 
                    % ERR_1 = rmse(F_DMD_1(t_extra_0:i+t_extra_0-1),F_AUG_1(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_1(t_extra_0:i+t_extra_0-1)))-min(abs(F_AUG_1(t_extra_0:i+t_extra_0-1))));
                    % ERR_2 = rmse(F_DMD_2(t_extra_0:i+t_extra_0-1),F_AUG_2(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_2(t_extra_0:i+t_extra_0-1)))-min(abs(F_AUG_2(t_extra_0:i+t_extra_0-1))));
                    % ERR_3 = rmse(F_DMD_3(t_extra_0:i+t_extra_0-1),F_AUG_3(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_3(t_extra_0:i+t_extra_0-1)))-min(abs(F_AUG_3(t_extra_0:i+t_extra_0-1))));
                    % ERR_4 = rmse(F_DMD_4(t_extra_0:i+t_extra_0-1),F_AUG_4(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_4(t_extra_0:i+t_extra_0-1)))-min(abs(F_AUG_4(t_extra_0:i+t_extra_0-1))));
                    % ERR_5 = rmse(F_DMD_5(t_extra_0:i+t_extra_0-1),F_AUG_5(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_5(t_extra_0:i+t_extra_0-1)))-min(abs(F_AUG_5(t_extra_0:i+t_extra_0-1))));


                    ERR_1 = rmse(F_DMD_1(t_extra_0:i+t_extra_0-1),F_AUG_1(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_1))-min(abs(F_AUG_1)));
                    ERR_2 = rmse(F_DMD_2(t_extra_0:i+t_extra_0-1),F_AUG_2(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_2))-min(abs(F_AUG_2)));
                    ERR_3 = rmse(F_DMD_3(t_extra_0:i+t_extra_0-1),F_AUG_3(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_3))-min(abs(F_AUG_3)));
                    ERR_4 = rmse(F_DMD_4(t_extra_0:i+t_extra_0-1),F_AUG_4(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_4))-min(abs(F_AUG_4)));
                    ERR_5 = rmse(F_DMD_5(t_extra_0:i+t_extra_0-1),F_AUG_5(t_extra_0:i+t_extra_0-1))/(max(abs(F_AUG_5))-min(abs(F_AUG_5)));

                    % ERR_1_c = ERR_1_c + ERR_1;
                    % ERR_2_c = ERR_2_c + ERR_2;
                    % ERR_3_c = ERR_3_c + ERR_3;
                    % ERR_4_c = ERR_4_c + ERR_4;
                    % ERR_5_c = ERR_5_c + ERR_5;
                    G1_ERR = [G1_ERR;ERR_1];
                    G2_ERR = [G2_ERR;ERR_2];
                    G3_ERR = [G3_ERR;ERR_3];
                    G4_ERR = [G4_ERR;ERR_4];
                    G5_ERR = [G5_ERR;ERR_5];

                    % G1_ERR = [G1_ERR;ERR_1_c];
                    % G2_ERR = [G2_ERR;ERR_2_c];
                    % G3_ERR = [G3_ERR;ERR_3_c];
                    % G4_ERR = [G4_ERR;ERR_4_c];
                    % G5_ERR = [G5_ERR;ERR_5_c];
                    
                end
            else
                ERR(1,jj) = rmse(F_DMD_1,F_AUG_1)/(max(abs(F_AUG_1))-min(abs(F_AUG_1)));
                ERR(2,jj) = rmse(F_DMD_2,F_AUG_2)/(max(abs(F_AUG_2))-min(abs(F_AUG_2)));
                ERR(3,jj) = rmse(F_DMD_3,F_AUG_3)/(max(abs(F_AUG_3))-min(abs(F_AUG_3)));
                ERR(4,jj) = rmse(F_DMD_4,F_AUG_4)/(max(abs(F_AUG_4))-min(abs(F_AUG_4)));
                ERR(5,jj) = rmse(F_DMD_5,F_AUG_5)/(max(abs(F_AUG_5))-min(abs(F_AUG_5)));
            end
    
    end
    
    if prediction_error ~= 1
        G1_ERR = [G1_ERR;ERR(1,jj)];
        G2_ERR = [G2_ERR;ERR(2,jj)];
        G3_ERR = [G3_ERR;ERR(3,jj)];
        G4_ERR = [G4_ERR;ERR(4,jj)];
        G5_ERR = [G5_ERR;ERR(5,jj)];
    end

end

toc
%%
% skip_first_no_steps = 10;
% if plot_NRMSE_prediction == 1
%     f1 = figure;
%     hold on
%     plot(t(t_extra_0+skip_first_no_steps:end),G1_ERR(skip_first_no_steps:end),'LineWidth',1.5,'Color','#66c2a5')
%     plot(t(t_extra_0+skip_first_no_steps:end),G2_ERR(skip_first_no_steps:end),'LineWidth',1.5,'Color','#fc8d62')
%     plot(t(t_extra_0+skip_first_no_steps:end),G3_ERR(skip_first_no_steps:end),'LineWidth',1.5,'Color','#8da0cb')
%     plot(t(t_extra_0+skip_first_no_steps:end),G4_ERR(skip_first_no_steps:end),'LineWidth',1.5,'Color','#e78ac3')
%     plot(t(t_extra_0+skip_first_no_steps:end),G5_ERR(skip_first_no_steps:end),'LineWidth',1.5,'Color','#a6d854')
%     hold off
%     set(gca,'TickLabelInterpreter','latex')
%     %grid on
%     leg = legend({'$G_1$','$G_2$','$G_3$','$G_4$','$G_5$'},'Interpreter','latex');
%     leg.Location = "southeast";
%     %xlabel('No. delay coordinates','Interpreter','latex')
%     %ylabel('$\displaystyle \left| \frac{RMS(G_{DMD})-RMS(G_n)}{RMS(G_n)}\right|\times 100\,\%$','Interpreter','latex')
% 
%     xlabel('Time [s]','Interpreter','latex')
%     if print_prediction_error == 1
%         filename = "NRMS_"+"win_"+int2str(window)+"_"+"s_"+int2str(s_0)+"predic"+int2str(extraTime)+"_"+name+".eps";
%         fullFilePath = fullfile(folderPath, filename);
%         print(gcf,fullFilePath,'-depsc','-vector');
%     end
%     ylabel('NRMSE','Interpreter','latex')
    
%     if print_NRMSE == 1
%         filename = "NRMS_"+"win_"+int2str(window)+"_"+"s_"+int2str(s_0)+"_"+int2str(s_f)+"_"+name+".eps";
% 
%         fullFilePath = fullfile(folderPath, filename);
%         exportgraphics(f1, fullFilePath, 'ContentType', 'vector', 'Resolution', 1500);
% 
%     end
% 
%     t_save = t(t_extra_0+skip_first_no_steps:end);
%     G1_ERR_save = G1_ERR(skip_first_no_steps:end);
%     G2_ERR_save = G2_ERR(skip_first_no_steps:end);
%     G3_ERR_save = G3_ERR(skip_first_no_steps:end);
%     G4_ERR_save = G4_ERR(skip_first_no_steps:end);
%     G5_ERR_save = G5_ERR(skip_first_no_steps:end);
% 
%     save("NRMSE_prediction_"+name+".mat",'t_save',...
%         'G1_ERR_save','G2_ERR_save','G3_ERR_save','G4_ERR_save','G5_ERR_save');
% end

%%
figure;

% for i=1:5:length(F_DMD_1)
%     hold on
%     plot(F_AUG_1(i),F_DMD_1(i),'*')
%     pause(0.5)
% end
%hold on
subplot(2,3,1)
plot(downsample(F_AUG_1,10),downsample(F_AUG_1,10),'Color','black')
hold on
plot(downsample(F_AUG_1,10),downsample(F_DMD_1,10),'diamond','Color','#66c2a5')
ylim([min(F_AUG_1) max(F_AUG_1)])
xlim([min(F_AUG_1) max(F_AUG_1)])
xlabel('Measured','Interpreter','latex')
ylabel('AugDMD','Interpreter','latex')
title('$\mathbf{G_1}$','Interpreter','latex','FontSize',11)
set(gca, 'TickLabelInterpreter', 'latex');
grid on

subplot(2,3,2)
plot(downsample(F_AUG_2,10),downsample(F_AUG_2,10),'Color','black')
hold on
plot(downsample(F_AUG_2,10),downsample(F_DMD_2,10),'diamond','Color','#fc8d62')
ylim([min(F_AUG_2) max(F_AUG_2)])
xlim([min(F_AUG_2) max(F_AUG_2)])
xlabel('Measured','Interpreter','latex')
ylabel('AugDMD','Interpreter','latex')
title('$\mathbf{G_2}$','Interpreter','latex','FontSize',11)
set(gca, 'TickLabelInterpreter', 'latex');
grid on

subplot(2,3,3)
plot(downsample(F_AUG_3,10),downsample(F_AUG_3,10),'Color','black')
hold on
plot(downsample(F_AUG_3,10),downsample(F_DMD_3,10),'diamond','Color','#8da0cb')
ylim([min(F_AUG_3) max(F_AUG_3)])
xlim([min(F_AUG_3) max(F_AUG_3)])
xlabel('Measured','Interpreter','latex')
ylabel('AugDMD','Interpreter','latex')
title('$\mathbf{G_3}$','Interpreter','latex','FontSize',11)
set(gca, 'TickLabelInterpreter', 'latex');
grid on

subplot(2,3,4)
plot(downsample(F_AUG_4,10),downsample(F_AUG_4,10),'Color','black')
hold on
plot(downsample(F_AUG_4,10),downsample(F_DMD_4,10),'diamond','Color','#e78ac3')
ylim([min(F_AUG_4) max(F_AUG_4)])
xlim([min(F_AUG_4) max(F_AUG_4)])
xlabel('Measured','Interpreter','latex')
ylabel('AugDMD','Interpreter','latex')
title('$\mathbf{G_4}$','Interpreter','latex','FontSize',11)
set(gca, 'TickLabelInterpreter', 'latex');
grid on

subplot(2,3,5)
plot(downsample(F_AUG_5,10),downsample(F_AUG_5,10),'Color','black')
hold on
plot(downsample(F_AUG_5,10),downsample(F_DMD_5,10),'diamond','Color','#a6d854')
ylim([min(F_AUG_5) max(F_AUG_5)])
xlim([min(F_AUG_5) max(F_AUG_5)])
xlabel('Measured','Interpreter','latex')
ylabel('AugDMD','Interpreter','latex')
title('$\mathbf{G_5}$','Interpreter','latex','FontSize',11)
set(gca, 'TickLabelInterpreter', 'latex');
grid on

% subplot(2,3,6)
% b = bar(ERR);
% b.FaceColor = 'flat';
% b.EdgeColor = 'none';
% b.CData(1,:) = [102,194,165]/255;
% b.CData(2,:) = [252,141,98]/255;
% b.CData(3,:) = [141,160,203]/255;
% b.CData(4,:) = [231,138,195]/255;
% b.CData(5,:) = [166,216,84]/255;
% xticklabels({'G1', 'G2', 'G3', 'G4', 'G5'});
% set(gca, 'TickLabelInterpreter', 'latex');
% title('\textbf{NRMSE}','Interpreter','latex');
% grid on
%%
if s_analysis == 1
hold on
    % plot(s_0:s_step:s_f,G1_ERR)
    % plot(s_0:s_step:s_f,G2_ERR)
    % plot(s_0:s_step:s_f,G3_ERR)
    % plot(s_0:s_step:s_f,G4_ERR)
    % plot(s_0:s_step:s_f,G5_ERR)
s_save = s_0:s_step:s_f;
save("NRMSE_"+name+"_s1_100_win5000_r5"+".mat",...
    'G1_ERR','G2_ERR','G3_ERR','G4_ERR','G5_ERR','s_save')
end
% 
if prediction_error == 1
hold on
    % plot(s_0:s_step:s_f,G1_ERR)
    % plot(s_0:s_step:s_f,G2_ERR)
    % plot(s_0:s_step:s_f,G3_ERR)
    % plot(s_0:s_step:s_f,G4_ERR)
    % plot(s_0:s_step:s_f,G5_ERR)
  a=   "NRMSE_"+name+"_prediction"+int2str(extraTime)+".mat";
%save("NRMSE_"+name+"_prediction"+int2str(extraTime)+".mat",...
%    'G1_ERR','G2_ERR','G3_ERR','G4_ERR','G5_ERR','t','t_extra_0')
end