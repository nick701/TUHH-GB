% Window size analysis (Approximation error)

%% DMD
clear all
%close all
tic
%%
% Load the raw data from the installation of the 3rd pile. For the 5 first
% geophones. Only vertical velocity is considered.
%
% Seq1 -> 1.55 - 2.00 min; penetration 0.9 - 1.0 m (approx)
% Seq2 -> 4.40 - 4.45 min; penetration 1.75 - 2.0 m (approx)
% Seq3 -> 10.55 - 11.0 min; penetration 3.50 - 3.60 (approx)
name = "mat-files/Seq1";
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

% Folder to save the plots
folderPath = 'E:\TUHH\Publications\in_preparation\2024_RASD\figures'; 
% ------
% DMD Parameters
% rank truncation
r = 5;

% ploting
plot_NRMSE_approx_error = 1;
print_NRMSE_approx_error = 0;

plot_geophones = 0;
% num of steps to move the window
nstep = 1;             

% number of delay coordinates. s=1 corresponds to standard DMD
s = 50;               

% size of the step. If step < window it will be an overlap in the data


% time considered for the extrapolation (prediction)
extraTime = 0;%0.1;

% index to indicate the begining of the analysis
initstep = 1;%1;

% creation of augmented matrix (to reduce computational time
plot_geophones = 0;

% time vector including extrapolation time


% data matrices accoring to DMD


% length of the extra time
%t_extra_0 = length(0:dt:(l-1)*dt);

% index to create store the error data
jj = 0;

% ------
% analysis 
% standard DMD
standard_DMD = 0;

% augmented DMD
aug_DMD = 1;
s_aug = 50;

% RMS 
G1_ERR = [];
G2_ERR = [];
G3_ERR = [];
G4_ERR = [];
G5_ERR = [];

win_0 = 100;
win_step = 200;
win_f = 9999;

for window = win_0:win_step:win_f
    l=ceil(window/2);
    f_aug = zeros(l,s*num_geo);
    step = window;
    t = 0:dt:(l-1)*dt+extraTime;
    X1 = zeros(s*num_geo,l-1);
    X2 = zeros(s*num_geo,l-1);
    for j = initstep:step:nstep*step+initstep-1
    
        jj = jj+1;
        for i = 1:s
            f = f_all(j:j+window,:);
            f_aug(:,[i s+i 2*s+i 3*s+i 4*s+i]) = [f(i:l+i-1,1) f(i:l+i-1,2) f(i:l+i-1,3) f(i:l+i-1,4) f(i:l+i-1,5)];
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
    
            ERR(1,jj) = rmse(F_DMD_1,F_AUG_1)/(max(abs(F_AUG_1))-min(abs(F_AUG_1)));
            ERR(2,jj) = rmse(F_DMD_2,F_AUG_2)/(max(abs(F_AUG_2))-min(abs(F_AUG_2)));
            ERR(3,jj) = rmse(F_DMD_3,F_AUG_3)/(max(abs(F_AUG_3))-min(abs(F_AUG_3)));
            ERR(4,jj) = rmse(F_DMD_4,F_AUG_4)/(max(abs(F_AUG_4))-min(abs(F_AUG_4)));
            ERR(5,jj) = rmse(F_DMD_5,F_AUG_5)/(max(abs(F_AUG_5))-min(abs(F_AUG_5)));
    
    end
    
    G1_ERR = [G1_ERR;ERR(1,jj)];
    G2_ERR = [G2_ERR;ERR(2,jj)];
    G3_ERR = [G3_ERR;ERR(3,jj)];
    G4_ERR = [G4_ERR;ERR(4,jj)];
    G5_ERR = [G5_ERR;ERR(5,jj)];


end

toc
%%
if plot_NRMSE_approx_error == 1
    f1 = figure;
    hold on
    plot(G1_ERR,'LineWidth',1.5,'Color','#66c2a5')
    plot(G2_ERR,'LineWidth',1.5,'Color','#fc8d62')
    plot(G3_ERR,'LineWidth',1.5,'Color','#8da0cb')
    plot(G4_ERR,'LineWidth',1.5,'Color','#e78ac3')
    plot(G5_ERR,'LineWidth',1.5,'Color','#a6d854')
    hold off
    set(gca,'TickLabelInterpreter','latex')
    grid on
    legend({'$G_1$','$G_2$','$G_3$','$G_4$','$G_5$'},'Interpreter','latex')
    xlabel('No. of samples','Interpreter','latex')
    %ylabel('$\displaystyle \left| \frac{RMS(G_{DMD})-RMS(G_n)}{RMS(G_n)}\right|\times 100\,\%$','Interpreter','latex')
    ylabel('NRMSE','Interpreter','latex')
    
    % if print_NRMSE == 1
    %     filename = "NRMSE_Approx_error_"+"win_"+int2str(win_0)+int2str(winf)+"_"+name+".eps";
    % 
    %     fullFilePath = fullfile(folderPath, filename);
    %     exportgraphics(f1, fullFilePath, 'ContentType', 'vector', 'Resolution', 1500);
    % 
    % end
end

