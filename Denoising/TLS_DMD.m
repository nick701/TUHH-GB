% Load data
load('mat-files/Seq1.mat');

% Combine signals into a matrix
X = [G1V3; G2V3; G3V3; G4V3; G5V3];

% Define parameters
dt = mean(diff(Time)); % Time step
noiseLevels = [0.01, 0.05, 0.1]; % Noise levels

% Storage for results
DenoisedData = cell(length(noiseLevels), 1);
Eigenvalues = cell(length(noiseLevels), 1);
Modes = cell(length(noiseLevels), 1);

for i = 1:length(noiseLevels)
    % Add Gaussian noise
    sigma = noiseLevels(i);
    X_noisy = X + sigma * randn(size(X));
    X_prime_noisy = X(:, 2:end) + sigma * randn(size(X(:, 2:end)));

    % Total Least Squares DMD
    Z = [X_noisy; X_prime_noisy];
    [U, S, V] = svd(Z, 'econ');
    r = rank(S); % Truncate based on energy threshold if needed
    U_r = U(:, 1:r);
    S_r = S(1:r, 1:r);
    V_r = V(:, 1:r);
    Z_r = U_r * S_r * V_r';

    % Separate back into X and X'
    d = size(X, 1);
    X_r = Z_r(1:d, :);
    X_prime_r = Z_r(d+1:end, :);

    % Koopman operator
    A_tls = X_prime_r / X_r;

    % Eigen decomposition
    [W, Lambda] = eig(A_tls);

    % DMD Modes
    Modes{i} = W;
    Eigenvalues{i} = diag(Lambda);

    % Reconstruct the signal
    b = W \ X(:, 1); % Initial condition projection
    t = Time - Time(1);
    X_reconstructed = W * (b .* exp(Eigenvalues{i} .* t));

    % Store reconstructed signal
    DenoisedData{i} = X_reconstructed;
end





% Visualization of results
for i = 1:length(noiseLevels)
    figure;
    plot(Time, X(1, :), 'k', 'LineWidth', 1.5); hold on;
    plot(Time, real(DenoisedData{i}(1, :)), '--r', 'LineWidth', 1.2);
    title(['Noise Level: ', num2str(noiseLevels(i))]);
    legend('Original', 'Denoised');
    xlabel('Time');
    ylabel('Amplitude');
end
