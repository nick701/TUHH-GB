% Define parameters
epsilon = 1e-1; % Perturbation parameter
dt = 1e-1; % Time step
tspan = 0:dt:10; % Time range
x0 = [1; 0]; % Initial condition

% Step 1: Define dynamics
dyn = @(t,x) ([0, 1+epsilon*t; -(1+epsilon*t), 0])*x;

% Step 2: Generate data by solving the ODE
[tq, xq] = ode45(dyn, tspan, x0); 
xq = xq'; tq = tq'; % Transpose for consistency

% Extract snapshot pairs
x = xq(:, 1:end-1); % Current state
y = xq(:, 2:end);   % Next state

% Step 3: Set the number of delay coordinates
% User-defined number of delay coordinates
num_delays = 1; % <-- Modify this value to test different numbers of delay coordinates

% Step 4: Define a function to construct Hankel matrix with delay embedding
function [X_delayed] = construct_delay_coordinates(x, num_delays)
    [n, m] = size(x); % n: number of rows (states), m: number of columns (time steps)
    X_delayed = [];
    for i = 1:(m-num_delays)
        delayed_snapshot = reshape(x(:,i:i+num_delays-1), [n*num_delays, 1]);
        X_delayed = [X_delayed, delayed_snapshot];
    end
end

% Construct delay coordinates for X and Y
X_delayed = construct_delay_coordinates(x, num_delays);
Y_delayed = construct_delay_coordinates(y, num_delays);

% Step 5: Perform DMD
% Compute the SVD of X_delayed
[U, S, V] = svd(X_delayed, 'econ');

% Compute the reduced-order system matrix
A_tilde = U' * Y_delayed * V / S;

% Compute eigenvalues and DMD modes
[W, D] = eig(A_tilde);
Phi = Y_delayed * V / S * W; % DMD modes

% Extract eigenvalues
eigenvalues = diag(D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Eigenvalues for %d delays:\n', num_delays);
disp(eigenvalues);


figure;
plot(real(eigenvalues), imag(eigenvalues), 'o');
xlabel('Real Part');
ylabel('Imaginary Part');
title(sprintf('Eigenvalues for %d Delays', num_delays));
grid on;
