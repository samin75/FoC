%% Part 1 — For-loop, fixed N (precision vs cost)
clear; clc; close all; rng(42);              % reproducible randomness

Ns = round(logspace(2,6,9));                 % N = 1e2 ... 1e6 (9 points)
pi_hat = zeros(size(Ns));
abs_err = zeros(size(Ns));
t_sec  = zeros(size(Ns));

for k = 1:numel(Ns)
    N = Ns(k);
    tic
    x = rand(N,1); 
    y = rand(N,1);
    H = sum(x.^2 + y.^2 <= 1);               % quarter-circle hits
    pi_hat(k) = 4*H/N;
    t_sec(k)  = toc;
    abs_err(k) = abs(pi_hat(k) - pi);        % deviation from true π
end

% 1) π-hat vs N
figure; 
plot(Ns, pi_hat, 'o-'); hold on; 
yline(pi,'--');
set(gca,'XScale','log'); grid on
xlabel('N (samples)','Interpreter','latex');
ylabel('$\hat{\pi}$ estimate','Interpreter','latex');
title('Monte Carlo $\pi$ vs N','Interpreter','latex');
legend({'$\hat{\pi}$','true $\pi$'},'Location','best','Interpreter','latex');

% 2) |error| vs N
figure; 
loglog(Ns, abs_err, 'o-'); grid on
xlabel('N (samples)','Interpreter','latex');
ylabel('Absolute error $|\hat{\pi}-\pi|$','Interpreter','latex');
title('Absolute error vs N','Interpreter','latex');

% 3) |error| vs runtime (precision–cost curve)
figure; 
loglog(t_sec, abs_err, 'o-'); grid on
xlabel('Runtime (s)','Interpreter','latex');
ylabel('Absolute error $|\hat{\pi}-\pi|$','Interpreter','latex');
title('Precision vs computational cost','Interpreter','latex');
