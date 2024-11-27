%% Forward/Backward DMD with Denoising
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
window = 500; % DMD window size
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
    
    time_dmd = (0:size(X1, 2) - 1) * dt;
    f_dmd = real(Phi * (b .* exp(omega * time_dmd)));
    
    f_denoised = max(0, abs(f_dmd) - lambda) .* sign(f_dmd);
    
    
    %% Visualization
    figure(1)
    subplot(2, 1, 1)
    plot(time(j:j+window-1), f_window(:, 1), 'DisplayName', 'Original')
    hold on
    plot(time(j:j+window-1), f_dmd(1, :), '--', 'DisplayName', 'DMD Reconstruction')
    legend
    title('Geophone 1 Signal Reconstruction')
    
    subplot(2, 1, 2)
    plot(time(j:j+window-1), f_window(:, 1) - f_denoised(1, :), 'r', 'DisplayName', 'Error')
    legend
    title('Denoising Error')
    
    drawnow
end

toc
