clear all; close all; clc;
%%


% Seq1 -> 1.55 - 2.00 min; penetration 0.9 - 1.0 m (approx)
% Seq2 -> 4.40 - 4.45 min; penetration 1.75 - 2.0 m (approx)
% Seq3 -> 10.55 - 11.0 min; penetration 3.50 - 3.60 (approx)

sequence = 2;

load ("Seq"+int2str(sequence)+".mat");

% Downsampling

down_sampling = 150;
G1 = downsample(G1V3,down_sampling);
G2 = downsample(G2V3,down_sampling);
G3 = downsample(G3V3,down_sampling);
G4 = downsample(G4V3,down_sampling);
G5 = downsample(G5V3,down_sampling);
time = downsample(Time,down_sampling);

dt = time(2)-time(1);

f = [G1'; G2'; G3'; G4'; G5'];

num_geo = 100;



ls =  ceil(length(f)/2);    
s = 1;
for i = 1:s
    f_aug([i s+i 2*s+i 3*s+i 4*s+i],:) = [f(1,i:ls+i-1) ;f(2,i:ls+i-1) ;f(3,i:ls+i-1) ;f(4,i:ls+i-1) ;f(5,i:ls+i-1)];
    %f_aug = f;
end
% f_aug = f;
x = f_aug(:,1:end-1); y = f_aug(:,2:end);
[n,m] = size(x);

%regularization
lambda = 1e-0;

odmd = OnlineDMD(n,0.9);
odmd.initialize(x(:,1:num_geo),y(:,1:num_geo),lambda);

% online DMD
t_prediction_horizon = 0.5;
figure;

for k = num_geo+1:m
    t_prediction = time(k):dt:t_prediction_horizon;


    odmd.update(x(:,k),y(:,k));
    [evals, modes] = odmd.computemodes();


    omega = log(evals) / dt;

    b = modes \ x(:,k);  % Use first snapshot

    f_dmd_update = real(modes * (b .* exp(omega*t_prediction)));  % Predicted state evolution

    clf

    subplot(3,2,1)
    plot(f(1,k:k+length(t_prediction)))
    hold on
    plot(f_dmd_update(1,:))
    legend( 'Measurements','ODMD')
    title('Geo 1')

    subplot(3,2,2)
    plot(f(2,k+1:k+length(t_prediction)))
    hold on
    plot(f_dmd_update(2,:))
    legend('Measurements','ODMD')

    title('Geo 2')

    subplot(3,2,3)
    plot(f(3,k+1:k+length(t_prediction)))
    hold on
    plot(f_dmd_update(3,:))
    legend('Measurements','ODMD')
    title('Geo 3')

    subplot(3,2,4)
    plot(f(4,k+1:k+length(t_prediction)))
    hold on
    plot(f_dmd_update(4,:))
    legend('Measurements','ODMD')
    title('Geo 4')

    subplot(3,2,5)
    plot(f(5,k+1:k+length(t_prediction)))
    hold on
    plot(f_dmd_update(5,:))
    legend('Measurements','ODMD')
    title('Geo 5')
    
    subplot(3,2,6)
    plot(real(omega),imag(omega),"o")
    xlabel('Re')
    ylabel('Img')
    title('Eigenvalues')
    axis equal
    grid on
    hold on

    drawnow 
    pause(.1)
end





