%% Combined DMD and Forward/Backward DMD Comparison
clear all;
close all;
tic;

% Load the dataset
name = "mat-files/Seq1"; % Adjust the path if necessary
load(name + ".mat");

% Parameters
n = 100; % Resampling rate
r = 5; % Rank truncation
window = 1000; % DMD window size
dt = 0.01; % Time step
nstep = 1; % Step size for window
initstep = 1; % Initial analysis step

% Downsample data
G1 = downsample(G1V3, n);
G2 = downsample(G2V3, n);
G3 = downsample(G3V3, n);
G4 = downsample(G4V3, n);
G5 = downsample(G5V3, n);
time = downsample(Time, n);

f_all = [G1'; G2'; G3'; G4'; G5']';

% Preallocate data for DMD comparison
f_dmd_results = zeros(size(f_all));
f_fwd_bwd_results = zeros(size(f_all));

%% Regular DMD Analysis
for j = initstep:nstep*window:length(f_all) - window
    f_window = f_all(j:j+window-1, :);
    
    % Regular DMD
    X1 = f_window(1:end-1, :)';
    X2 = f_window(2:end, :)';
    [U, S, V] = svd(X1, 'econ');
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
    
    % Store results for comparison
    f_dmd_results(j:j+window-1, :) = f_dmd.';
end

%% Forward/Backward DMD Analysis
for j = initstep:nstep*window:length(f_all) - window
    f_window = f_all(j:j+window-1, :);
    
    % Forward/Backward DMD with denoising
    X1 = f_window(1:end-1, :)';
    X2 = f_window(2:end, :)';
    [U, S, V] = svd(X1, 'econ');
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
    
    % Apply denoising
    lambda_denoise = 0.1; % Denoising threshold
    f_denoised = max(0, abs(f_dmd) - lambda_denoise) .* sign(f_dmd);
    
    % Store results for comparison
    f_fwd_bwd_results(j:j+window-1, :) = f_denoised.';
end

%% Plot Comparison
figure;
% Plot all geophones combined for each method
plot(time, mean(f_all, 2), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original');
hold on;
plot(time, mean(f_dmd_results, 2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Regular DMD');
plot(time, mean(f_fwd_bwd_results, 2), 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Fwd/Back DMD');
title('Comparison of DMD Methods (All Geophones)');
xlabel('Time [s]');
ylabel('Amplitude');
legend;
grid on;

toc;
