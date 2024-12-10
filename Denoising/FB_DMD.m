%% FB DMD 

clear all
close all
tic

name = "mat-files/Seq1";
load(name + ".mat")

n = 100;
G1 = downsample(G1V3, n);
G2 = downsample(G2V3, n);
G3 = downsample(G3V3, n);
G4 = downsample(G4V3, n);
G5 = downsample(G5V3, n);
time = downsample(Time, n);

% Calculate time delta
dt = time(2) - time(1);

f_all = [G1'; G2'; G3'; G4'; G5']';

% Parameters
r = 5; % Truncation rank
plot_live = 1; % Enable live plotting
plot_error = 1; % Enable error plotting

window = 500; % Window length
step = 100; % Step size for overlapping windows

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
        for geo = 1:size(Y, 1)
            subplot(3, 2, geo);
            plot(time_recon, Y(geo, :), 'b', 'DisplayName', 'Original Signal');
            hold on;
            plot(time_recon, Y_recon(geo, :), 'r--', 'DisplayName', 'Reconstructed Signal');
            title(['Geophone ', num2str(geo)]);
            xlabel('Time [s]');
            ylabel('Velocity [m/s]');
            legend();
            hold off;
        end
        drawnow;
    end
    
    % Compute and plot error
    if plot_error
        NRMSE = sqrt(mean((Y(:) - Y_recon(:)).^2)) / (max(Y(:)) - min(Y(:)));
        subplot(3, 2, 6);
        plot(start_idx, NRMSE, 'ko', 'MarkerFaceColor', 'r');
        hold on;
        xlabel('Time Window Start');
        ylabel('NRMSE');
        title('Error Analysis');
        drawnow;
    end
end

toc
