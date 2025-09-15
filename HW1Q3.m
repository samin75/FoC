% Q3 — Demo of the mcPi_sigfigs function
% Calls the function with user-defined precision

[pi_out, N_used, ci95] = mcPi_sigfigs(3, 'Batch', 2000, 'Plot', true);

fprintf('Returned pi = %s, N = %d, CI = [%.6f, %.6f]\n', ...
        num2str(pi_out), N_used, ci95(1), ci95(2));

mcPi_sigfigs(2, 'Batch', 2000, 'Plot', true, 'Figure', 1);
mcPi_sigfigs(3, 'Batch', 2000, 'Plot', true, 'Figure', 2);
