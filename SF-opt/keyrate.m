function R = keyrate(mu, nu, cnt)

    %% initialize

    %  protocol-related efficiency 
    q = 0.5;
    %  attenuation coefficient (in 1550nm fiber optics) (dB/km)
    alpha = 0.2;
    %  misalignment error rate 
    e_d = 0.02;

    %  symbol error rate of vacuum pulse
    e_0 = 0.5;
    %  response rate of vacuum pulse
    Y_0 = 10^-6;

    %  detect efficiency at bob's side
    eta_bob = 0.5;

    % %  average photon number
    % mu = 0.5;
    % %  decoy state strength
    % nu = 0.05;

    % %%--------------------------------------------------------------------------------------------

    %% 3-decoy state with statistical fluctutations

    % parameter for base pulse number
    N = 10^10;

    % parameter response rate in single photon, statistic fluctutations included  
    gamma = 5.3;

    % Probabilitites
    % photon strength selection
    Pr_mu = 0.6;
    Pr_nu = 0.3;
    Pr_o  = 0.1;
    %% total photon number of different decoy states
    N_mu = Pr_mu * N;
    % N_mu = 6 * 10^9;
    N_nu = Pr_nu * N;
    % N_nu = 3 * 10^9;
    N_o  = Pr_o * N;
    % N_o  = 10^9;
    % base selection
    Pr_Z = 0.5;
    Pr_X = 0.5;
    %% the number of specific base and decoy state
    NZ_mu = N_mu * Pr_Z * Pr_Z;
    % NZ_mu = 1.5 * 10^9;
    NZ_nu = N_nu * Pr_Z * Pr_Z;
    % NZ_nu = 7.5 * 10^8;
    NX_nu = N_nu * Pr_X * Pr_X;
    % NX_nu = 7.5 * 10^8;
    %% the upper and lower limit of Y_0, the background rate
    YL_0 = Y_0 * (1 - gamma/sqrt(N_o*Y_0));
    YU_0 = Y_0 * (1 + gamma/sqrt(N_o*Y_0));

    %% probability of i-photon(i=0,1,2) with strength mu or nu
    Pmu_0 = poisspdf(0, mu); %% (P_0)^mu
    Pmu_1 = poisspdf(1, mu); %% (P_1)^mu
    Pmu_2 = poisspdf(2, mu); %% (P_2)^mu
    Pnu_0 = poisspdf(0, nu); %% (P_0)^nu
    Pnu_1 = poisspdf(1, nu); %% (P_1)^nu
    Pnu_2 = poisspdf(2, nu); %% (P_2)^nu

    %% efficiency
    eta_channel = power(10, -alpha*cnt/10);
    eta = eta_channel * eta_bob;

    %% decoy state
    %  gain 
    Q_nu = Y_0 + 1 - exp(-eta*nu);
    %  qubit error rate
    T_nu = e_0 * Y_0 + e_d * (1-exp(-eta*nu));

    %% signal state
    %  gain 
    Q_mu = Y_0 + 1 - exp(-eta*mu);
    %  qubit error rate
    T_mu = e_0 * Y_0 + e_d * (1-exp(-eta*mu));
    E_mu = T_mu / Q_mu;

    DELTA_Z_nu = gamma / sqrt(NZ_nu * Q_nu);
    DELTA_Z_mu = gamma / sqrt(NZ_mu * Q_mu);
    DELTA_PRIME_X_nu = gamma / sqrt(NZ_nu * T_nu);

    %% single-photon rate
    SF_YZL_1 = (Pmu_2*Q_nu*(1 - DELTA_Z_nu) - Pnu_2*Q_mu*(1 - DELTA_Z_mu) + (Pnu_2*Pmu_0 - Pmu_2*Pnu_0)*YL_0) / (Pnu_1*Pmu_2 - Pmu_1*Pnu_2);

    %% single-photon error rate
    SF_eXU_1 = (T_nu*(1 + DELTA_PRIME_X_nu) - Pnu_0*e_0*YL_0) / (Pnu_1 * SF_YZL_1);

    %% Binary Shannon Entropy
    H_21 = -log2(SF_eXU_1)*SF_eXU_1 - (1-SF_eXU_1)*log2(1-SF_eXU_1);
    H_22 = -log2(E_mu)*E_mu - (1-E_mu)*log2(1-E_mu);

    %% single-photon key generation rate
    SF_R = q * (Pmu_1 * SF_YZL_1 * (1 - H_21) - Q_mu * 1.16 * H_22);

    R = SF_R;

% %%--------------------------------------------------------------------------------------------
end