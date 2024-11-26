%% TLS-DMD
clear all;
close all;
clc;
tic;

name = 'mat-files/Seq1';
load(name + ".mat");

%% Parameters
n = 100; 
G1 = downsample(G1V3, n);
G2 = downsample(G2V3, n);
G3 = downsample(G3V3, n);
G4 = downsample(G4V3, n);
G5 = downsample(G5V3, n);
time = downsample(Time, n);

dt = time(2) - time(1); 

f_all = [G1'; G2'; G3'; G4'; G5']';

num_geo = 5;    % Number of geophones
r = 5;          % Rank truncation
window = 500;   % Window size
step = window;  % Sliding step for window
denoise = true; % Enable denoising
extraTime = 0;  % Time for prediction

%% TLS-DMD Analysis
[rows, cols] = size(f_all);
t = 0:dt:(window-1)*dt+extraTime;
errors = zeros(num_geo, cols-window);

for j = 1:step:(cols-window)
    X = f_all(j:j+window-1, :)';
    X1 = X(:, 1:end-1);
    X2 = X(:, 2:end);

    % SVD for Total Least Squares
    [U, S, V] = svd([X1; X2], 'econ');
    rTLS = min(r, size(S, 2));
    U_r = U(:, 1:rTLS);
    S_r = S(1:rTLS, 1:rTLS);
    V_r = V(:, 1:rTLS);

    X1_denoised = U_r(1:size(X1, 1), :) * S_r * V_r';
    X2_denoised = U_r(size(X1, 1)+1:end, :) * S_r * V_r';

    Atilde = X1_denoised \ X2_denoised;
    [W, D] = eig(Atilde);
    Phi = X2_denoised * V_r / S_r * W;
   
    lambda = diag(D);
    omega = log(lambda) / dt;

   
    b = Phi \ X1_denoised(:, 1);

    % Reconstruct signals
    t_dmd = 0:dt:(size(X1, 2)-1)*dt;
    f_dmd = real(Phi * (b .* exp(omega * t_dmd)));

    % Denoising application
    if denoise
        f_all(j:j+size(f_dmd, 2)-1, :) = f_dmd';
    end

    errors(:, j:j+size(f_dmd, 2)-1) = abs(f_dmd - X1_denoised);

    % Visualization
    if j == 1
        figure;
        for i = 1:num_geo
            subplot(ceil(num_geo/2), 2, i);
            plot(t, f_all(j:j+window-1, i), 'k', 'DisplayName', 'Raw');
            hold on;
            plot(t, f_dmd(i, :), 'r--', 'DisplayName', 'DMD');
            legend('show');
            title(['Geophone ', num2str(i)]);
            xlabel('Time [s]');
            ylabel('Amplitude');
        end
    end
end

%% Plot Errors
figure;
for i = 1:num_geo
    subplot(ceil(num_geo/2), 2, i);
    plot(errors(i, :), 'b');
    title(['Error for Geophone ', num2str(i)]);
    xlabel('Time');
    ylabel('Error');
end

toc;
