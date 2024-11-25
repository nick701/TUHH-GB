% Load the dataset
load('mat-files/Seq1.mat'); % Load the .mat file containing the signal data

% Combine the relevant signal variables into a single data matrix (assume G1V3, G2V3, etc.)
data = [G1V3, G2V3, G3V3, G4V3, G5V3];  % Combine signals column-wise

% Prepare the state matrices
X = data(:, 1:end-1);  % Original states
X_prime = data(:, 2:end);  % Advanced states

% Forward DMD
[U, S, V] = svd(X, 'econ'); % Singular Value Decomposition
A_tilde = U' * X_prime * V * diag(1 ./ diag(S)); % Approximation of A
[eigVec, eigVal] = eig(A_tilde);
Phi = X_prime * V * diag(1 ./ diag(S)) * eigVec; % DMD modes
lambda = diag(eigVal); % DMD eigenvalues
omega = log(lambda); % Continuous-time eigenvalues

% Backward DMD
X_rev = fliplr(X);  % Reverse the time dimension
X_prime_rev = fliplr(X_prime);

[U_b, S_b, V_b] = svd(X_rev, 'econ');
A_tilde_b = U_b' * X_prime_rev * V_b * diag(1 ./ diag(S_b));
[eigVec_b, eigVal_b] = eig(A_tilde_b);
Phi_b = X_prime_rev * V_b * diag(1 ./ diag(S_b)) * eigVec_b;
lambda_b = diag(eigVal_b);
omega_b = log(lambda_b);

% Reconstruct forward and backward signals
time_vec = 1:size(data, 2);
X_denoised_fwd = real(Phi * diag(exp(omega * time_vec)) * eigVec');
X_denoised_bwd = real(Phi_b * diag(exp(omega_b * time_vec)) * eigVec_b');

% Combine denoised signals (e.g., simple averaging)
X_denoised = 0.5 * (X_denoised_fwd + fliplr(X_denoised_bwd));

% Plot results
figure;
subplot(3,1,1);
plot(data');
title('Original Signal');

subplot(3,1,2);
plot(X_denoised_fwd');
title('Forward Denoised Signal');

subplot(3,1,3);
plot(X_denoised');
title('Combined Denoised Signal');
