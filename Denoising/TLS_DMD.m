% Load data
load('mat-files/Seq1.mat');  % Load the dataset

% Combine signals into a matrix
% Each signal is a row in the matrix
X_original = [G1V3'; G2V3'; G3V3'; G4V3'; G5V3'];  % Transpose to make time snapshots columns

% Extract snapshots X and X' from the dataset
X = X_original(:, 1:end-1);   % Past states
X_prime = X_original(:, 2:end); % Future states

% Denoising function using Singular Value Thresholding (SVT)
function X_denoised = denoise(X)
    [U, S, V] = svd(X, 'econ');
    % Thresholding: keep singular values above a certain threshold
    threshold = 0.1 * max(diag(S)); % Example: 10% of maximum singular value
    S_denoised = diag(max(diag(S) - threshold, 0)); % Apply soft-thresholding
    X_denoised = U * S_denoised * V';
end

% Apply denoising to X and X_prime
X_denoised = denoise(X);
X_prime_denoised = denoise(X_prime);

% Concatenate denoised data for TLS
Z = [X_denoised; X_prime_denoised];

% Perform Singular Value Decomposition
[U, S, V] = svd(Z, 'econ');

% Partition U, S, V for TLS
n = size(X, 1); % Number of rows in X
m = size(X_prime, 1); % Number of rows in X_prime
r = rank(S); % Rank of the S matrix

% Ensure reduced partitions are extracted correctly
U_r = U(1:n, 1:r);    % Top part of U with r columns
S_r = S(1:r, 1:r);    % Top left part of S
V_r = V(:, 1:r);      % First r columns of V

% Compute TLS-DMD operator
A_TLS = X_prime_denoised * V_r * pinv(S_r) * U_r';

% Eigen decomposition of A_TLS
[W, D] = eig(A_TLS);

% Reconstruct signal
timesteps = size(X, 2); % Number of timesteps
initial_state = X_denoised(:, 1); % Initial condition
X_reconstructed = zeros(size(X));
for t = 1:timesteps
    X_reconstructed(:, t) = W * (D^(t-1)) * pinv(W) * initial_state;
end

% Plot the original and reconstructed signals
figure;
subplot(2, 1, 1);
plot(1:timesteps, X(1, :), 'r', 'LineWidth', 1.5); hold on;
plot(1:timesteps, X_reconstructed(1, :), 'b--', 'LineWidth', 1.5);
title('Original vs Reconstructed Signal - Component 1');
xlabel('Time Step');
ylabel('Amplitude');
legend('Original', 'Reconstructed');

subplot(2, 1, 2);
plot(1:timesteps, X(2, :), 'r', 'LineWidth', 1.5); hold on;
plot(1:timesteps, X_reconstructed(2, :), 'b--', 'LineWidth', 1.5);
title('Original vs Reconstructed Signal - Component 2');
xlabel('Time Step');
ylabel('Amplitude');
legend('Original', 'Reconstructed');
grid on;

% Compute and display Signal-to-Noise Ratio (SNR)
reconstruction_error = norm(X - X_reconstructed, 'fro');
snr_value = 20 * log10(norm(X, 'fro') / reconstruction_error);
disp(['Reconstruction SNR: ', num2str(snr_value), ' dB']);
