%% init vectors
R_SF = zeros(251, 1);

%% optimize parameters
[mu_init, nu_init] = init_parameters(1);
% mu_init = 0.716
% nu_init = 0.008
[mu, nu, R_SF_opt] = optimize_parameters(mu_init, nu_init);

%% from 1km to 251km, calculate key rate with statistical fluctutation
for L = 1:1:251
    R_SF(L) = keyrate(0.5, 0.05, L);
end

%% plotting
L = 1:1:251;

figure
plot(L, mu')
hold on
plot(L, nu')
hold off
legend('mu', 'nu')

figure
semilogy(L, R_SF_opt')
hold on 
semilogy(L, R_SF', '.')
hold off
legend('optimized', 'original')