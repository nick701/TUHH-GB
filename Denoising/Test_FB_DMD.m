%% FB-DMD Implementation

clear all
close all
tic

% Load the raw data
name = "mat-files/Seq1";
load(name + ".mat")

% Resample the data to avoid aliasing
n = 100;
G1 = downsample(G1V3, n);
G2 = downsample(G2V3, n);
G3 = downsample(G3V3, n);
G4 = downsample(G4V3, n);
G5 = downsample(G5V3, n);
time = downsample(Time, n);

% Calculate time delta
dt = time(2) - time(1);

% Arrange data in matrix
f_all = [G1'; G2'; G3'; G4'; G5']';

% Parameters
r = 5; % Truncation rank
plot_live = 1; % Enable live plotting
plot_error = 1; % Enable error plotting

% FB-DMD Implementation
window = 500; % Window length
step = 100; % Step size for overlapping windows

% Initialize NRMSE storage
nrmse_values = [];
start_indices = [];

for start_idx = 1:step:(size(f_all, 1) - window)
    % Extract data windows
    Y = f_all(start_idx:start_idx+window-1, :)';
    Yp = Y(:, 2:end);
    Ym = Y(:, 1:end-1);
    
    % Forward DMD
    [U_f, S_f, V_f] = svd(Ym, 'econ');
    f_Atilde = U_f' * Yp * V_f / S_f;
    
    % Backward DMD
    [U_b, S_b, V_b] = svd(Yp, 'econ');
    b_Atilde = U_b' * Ym * V_b / S_b;
    
    % Combined FB-DMD
    Atilde = (f_Atilde / b_Atilde) ^ 0.5;
    [W, D] = eig(Atilde);
    Phi = Yp * V_f / S_f * W;
    
    % DMD Spectra
    lambda = diag(D);
    omega = log(lambda) / dt;
    
    % Initial condition
    b = Phi \ Ym(:, 1);
    
    % Reconstructed signal
    time_recon = 0:dt:(window-1)*dt;
    Y_recon = real(Phi * (b .* exp(omega * time_recon)));
    
    % Live plotting
    if plot_live
        figure(1);
        clf;
        hold on;
        for geo = 1:size(Y, 1)
            plot(time_recon, Y(geo, :), 'DisplayName', ['Geophone ', num2str(geo), ' Original']);
            plot(time_recon, Y_recon(geo, :), '--', 'DisplayName', ['Geophone ', num2str(geo), ' Reconstructed']);
        end
        xline(window * dt, '--', {'Predicted', 'Signal'});
        title('Geophones: Original vs. Reconstructed Signals');
        xlabel('Time [s]');
        ylabel('Velocity [m/s]');
        legend('show');
        hold off;
        drawnow;
    end
    
    % Compute and store NRMSE
    NRMSE = sqrt(mean((Y(:) - Y_recon(:)).^2)) / (max(Y(:)) - min(Y(:)));
    nrmse_values = [nrmse_values, NRMSE];
    start_indices = [start_indices, start_idx];
end

% Plot NRMSE
if plot_error
    figure(2);
    plot(start_indices, nrmse_values, '-r', 'LineWidth', 1.5);
    title('NRMSE Over Time Windows');
    xlabel('Window Start Index');
    ylabel('NRMSE');
    grid on;
end

toc
