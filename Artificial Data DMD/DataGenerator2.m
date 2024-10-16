% define dynamics
epsilon = 1e-1;
dyn = @(t,x) ([0, 1+epsilon*t; -(1+epsilon*t),0])*x;

% generate data

% Increase tspan or reduce dt to generate more data
dt = 1e-2;
tspan = 0:dt:10;

x0 = [1;0];
[tq,xq] = ode45(dyn, tspan, x0);

% add small random noise to the data
noise_level = 1e-3;  % adjust this value as needed for testing
xq = xq + noise_level * randn(size(xq)); % add noise to the signal

% extract snapshot pairs
xq = xq'; tq = tq';
x = xq(:,1:end-1); y = xq(:,2:end); time = tq(2:end);

% true dynamics, eigenvalues
[n, m] = size(x);
A = zeros(n,n,m);
evals = zeros(n,m);

for k = 1:m
    A(:,:,k) = [0, 1+epsilon*time(k); -(1+epsilon*time(k)),0]; % continuous time dynamics
    evals(:,k) = eig(A(:,:,k)); % analytical continuous time eigenvalues
end
