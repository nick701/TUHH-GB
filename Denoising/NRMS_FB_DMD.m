%% Forward/Backward DMD with Denoising and Continuous NRMSE
clear all
close all
tic

%% Load the dataset
name = "mat-files/Seq1";
load(name + ".mat")

% Parameters
n = 100; % Resampling rate
num_geo = 5; % Number of geophones
r = 5; % Rank truncation
window = 1000; % DMD window size
nstep = 1; % Step size for window
initstep = 1; % Initial analysis step
dt = 0.01; % Time step

G1 = downsample(G1V3, n);
G2 = downsample(G2V3, n);
G3 = downsample(G3V3, n);
G4 = downsample(G4V3, n);
G5 = downsample(G5V3, n);
time = downsample(Time, n);

f_all = [G1'; G2'; G3'; G4'; G5']';

X1 = zeros(num_geo, window - 1);
X2 = zeros(num_geo, window - 1);

lambda = 0.1;

%% DMD Analysis
figure;
for j = initstep:nstep*window:length(f_all) - window
    f_window = f_all(j:j+window-1, :);
    
    % Forward DMD
    X1 = f_window(1:end-1, :)';
    X2 = f_window(2:end, :)';
    
    % Apply SVD
    [U, S, V] = svd(X1, 'econ');
    
    % Truncate SVD
    Ur = U(:, 1:r);
    Sr = S(1:r, 1:r);
    Vr = V(:, 1:r);
    
    Atilde = Ur' * X2 * Vr / Sr;
    
    [W, D] = eig(Atilde);
    Phi = X2 * Vr / Sr * W;
    lambda = diag(D);
    omega = log(lambda) / dt;
    
    b = Phi \ X1(:, 1);
    
    time_dmd = (0:window-1) * dt;
    f_dmd = real(Phi * (b .* exp(omega * time_dmd)));
    
    f_denoised = max(0, abs(f_dmd) - lambda) .* sign(f_dmd);
    
    %% NRMSE Calculation
    ERR = zeros(num_geo, window);
    for g = 1:num_geo
        F_DMD = real(f_dmd(g, :)');
        F_AUG = f_window(:, g);
        % Compute RMSE
        rmse_val = sqrt(mean((F_DMD - F_AUG).^2, 2));
        % Normalize RMSE
        ERR(g, :) = rmse_val ./ (max(abs(F_AUG)) - min(abs(F_AUG)));
    end
    nrmse_over_time = mean(ERR, 1); % Average NRMSE for all geophones over time
    
    %% Visualization
    
    subplot(2, 1, 1)
    plot(time(j:j+window-1), f_window(:, 1), 'DisplayName', 'Original')
    hold on
    plot(time(j:j+window-1), f_dmd(1, :), '--', 'DisplayName', 'DMD Reconstruction')
    legend
    title('Geophone 1 Signal Reconstruction')
    
    subplot(2, 1, 2)
    plot(time(j:j+window-1), nrmse_over_time, 'r', 'DisplayName', 'NRMSE')
    legend
    title('NRMS Error')
    xlabel('Time')
    ylabel('NRMSE')
    
    drawnow
    pause(1)
    clf
end

toc
