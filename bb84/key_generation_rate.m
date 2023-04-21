function R = key_generation_rate(Y_1, e_1, Q_mu, E_mu)
    q = 0.5;
    % (P_1)^mu
    Pmu_1 = mu*exp(-mu);
    Q_1 = Pmu_1 * Y_1;
    R = q * (Q_1 * (1-binary_shannon_entropy(e_1)) - Q_mu * error_correction_rate(E_mu) * binary_shannon_entropy(E_mu) ) 
end