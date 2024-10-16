run('DataGenerator.m'); % Load 'x', 'y', 'time', etc. from DataGenerator.m
% run('DataGenerator2.m'); % Load 'x', 'y', 'time', etc. from DataGenerator.m


% Delay Coordinates
num_delays = 150; 


function [X_delayed] = construct_delay_coordinates(x, num_delays)
    [n, m] = size(x); 
    X_delayed = [];
    for i = 1:(m-num_delays)
        delayed_snapshot = reshape(x(:,i:i+num_delays-1), [n*num_delays, 1]);
        X_delayed = [X_delayed, delayed_snapshot];
    end
end


X_delayed = construct_delay_coordinates(x, num_delays);
Y_delayed = construct_delay_coordinates(y, num_delays);


% SVD
[U, S, V] = svd(X_delayed, 'econ');

A_tilde = U' * Y_delayed * V / S;

% Eigenvalues and DMD modes
[W, D] = eig(A_tilde);
Phi = Y_delayed * V / S * W; % DMD modes


eigenvalues = diag(D);


fprintf('Eigenvalues for %d delays:\n', num_delays);
disp(eigenvalues);

% Eigenvalue Plot
figure;
plot(real(eigenvalues), imag(eigenvalues), 'o', ...
    'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
xlabel('Re');
ylabel('Im');
title(sprintf('Eigenvalues for %d Delays', num_delays));
grid on;