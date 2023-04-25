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

%  average photon number
mu = 0.5;
%  decoy state strength
nu = 0.05;

%%--------------------------------------------------------------------------------------------

%% infinite decoy state bb84 simulate

%% init empty vector
Y_1_inf = zeros(1,250);
e_1_inf = zeros(1,250);
R_inf = zeros(1,250);

for cnt = 1:1:250

    %% efficiency
    eta_channel = power(10, -alpha*cnt/10);
    eta = eta_channel * eta_bob;
    % efficiency for 1-photon state
    eta_1 = 1 - power(1-eta, 1);

    %% decoy state
    %  gain 
    Q_mu = Y_0 + 1 - exp(-eta*mu);
    %  qubit error rate
    T_mu = e_0*Y_0 + e_d*(1-exp(-eta*mu));
    E_mu = T_mu / Q_mu;

    %% single-photon rate
    Y_1_inf(1,cnt) = Y_0 + eta_1 - Y_0*eta_1;

    %% single-photon error rate
    e_1_inf(1,cnt) = (e_0*Y_0 + e_d*eta_1) / Y_1_inf(1,cnt);

    %% single-photon key generation rate
    % 1-photon gain
    Q_1 = Y_1_inf(1,cnt)*poisspdf(1,mu);
    R_inf(1,cnt) = q * (-Q_mu*error_correction_rate(E_mu)*binary_shannon_entropy(E_mu) + Q_1*(1-binary_shannon_entropy(e_1_inf(1,cnt))));

end

% %%--------------------------------------------------------------------------------------------

%% 3-decoy-state bb84 simulate

%% vaccum state 
%  gain
Q_0 = Y_0;
%  qubit error rate
T_0 = Y_0 * e_0;

%% probability of i-photon(i=0,1,2) with strength mu or nu
Pmu_0 = poisspdf(0, mu); %% (P_0)^mu
Pmu_1 = poisspdf(1, mu); %% (P_1)^mu
Pmu_2 = poisspdf(2, mu); %% (P_2)^mu
Pnu_0 = poisspdf(0, nu); %% (P_0)^nu
Pnu_1 = poisspdf(1, nu); %% (P_1)^nu
Pnu_2 = poisspdf(2, nu); %% (P_2)^nu

%% init empty vector
Y_1 = zeros(1,250);
e_1 = zeros(1,250);
R = zeros(1,250);

for cnt = 1:1:250

    %% efficiency
    eta_channel = power(10, -alpha*cnt/10);
    eta = eta_channel * eta_bob;

    %% decoy state
    %  gain 
    Q_nu = overall_qubit_gain(nu, Y_0, eta);
    %  qubit error rate
    T_nu = overall_qubit_error_rate(nu, Y_0, eta);

    %% signal state
    %  gain 
    Q_mu = overall_qubit_gain(mu, Y_0, eta);
    %  qubit error rate
    T_mu = overall_qubit_error_rate(mu, Y_0, eta);
    E_mu = average_qubit_error_rate(Q_mu, T_mu);

    %% single-photon rate
    Y_1(1,cnt) = (Pmu_2*(Q_nu - Pnu_0*Y_0) - Pnu_2*(Q_mu - Pmu_0*Y_0)) / (Pnu_1*Pmu_2 - Pmu_1*Pnu_2);
    % Y_1_alt(1,cnt) = mu/(mu*nu - nu*nu) * (Q_nu*exp(nu) - Q_mu*exp(mu)*(nu*nu/(mu*mu)) - ((mu*mu - nu*nu)/(mu*mu))*Y_0)

    %% single-photon error rate
    e_1(1,cnt) = (T_nu - Pnu_0*Y_0*e_0) / (Pnu_1 * Y_1(1,cnt));
    % e_1_alt(1,cnt) = (average_qubit_error_rate(Q_nu, T_nu)*Q_nu*exp(nu) - e_0*Y_0)/(Y_1_alt(1,cnt)*nu)

    %% single-photon key generation rate
    %  single-photon gain
    Q_1 = Y_1(1,cnt)*poisspdf(1,mu);
    
    R(1,cnt) = q * (Q_1*(1-binary_shannon_entropy(e_1(1,cnt))) - Q_mu*error_correction_rate(E_mu)*binary_shannon_entropy(E_mu));

end

% %%--------------------------------------------------------------------------------------------

%% 3-decoy state with statistical fluctutations

% parameter for base pulse number
N = 10^10;

% parameter response rate in single photon, statistic fluctutations included  
Gamma = 5.3;

% Probabilitites
% photon strength selection
Pr_mu = 0.6;
Pr_nu = 0.3;
Pr_o  = 0.1;
%% total photon number of different decoy states
N_mu = Pr_mu * N;
N_nu = Pr_nu * N;
N_o  = Pr_o * N;
% base selection
Pr_Z = 0.5;
Pr_X = 0.5;
%% the number of specific base and decoy state
NZ_mu = N_mu * Pr_Z * Pr_Z;
NZ_nu = N_nu * Pr_Z * Pr_Z;
NX_nu = N_nu * Pr_X * Pr_X;
%% the upper and lower limit of Y_0, the background rate
YL_0 = Y_0 * (1 - Gamma/sqrt(N_o*Y_0));
YU_0 = Y_0 * (1 + Gamma/sqrt(N_o*Y_0));

%% probability of i-photon(i=0,1,2) with strength mu or nu
Pmu_0 = poisspdf(0, mu); %% (P_0)^mu
Pmu_1 = poisspdf(1, mu); %% (P_1)^mu
Pmu_2 = poisspdf(2, mu); %% (P_2)^mu
Pnu_0 = poisspdf(0, nu); %% (P_0)^nu
Pnu_1 = poisspdf(1, nu); %% (P_1)^nu
Pnu_2 = poisspdf(2, nu); %% (P_2)^nu

%% init empty vector
SF_YZL_1 = zeros(1,250);
SF_eXU_1 = zeros(1,250);
SF_R = zeros(1,250);

for cnt = 1:1:250

    %% efficiency
    eta_channel = power(10, -alpha*cnt/10);
    eta = eta_channel * eta_bob;

    %% decoy state
    %  gain 
    Q_nu = overall_qubit_gain(nu, Y_0, eta);
    %  qubit error rate
    T_nu = overall_qubit_error_rate(nu, Y_0, eta);

    %% signal state
    %  gain 
    Q_mu = overall_qubit_gain(mu, Y_0, eta);
    %  qubit error rate
    T_mu = overall_qubit_error_rate(mu, Y_0, eta);
    E_mu = average_qubit_error_rate(Q_mu, T_mu);

    DELTA_Z_nu = Gamma / sqrt(NZ_nu * Q_nu);
    DELTA_Z_mu = Gamma / sqrt(NZ_mu * Q_mu);
    DELTA_PRIME_X_nu = Gamma / sqrt(NZ_nu * T_nu);

    %% single-photon rate
    SF_YZL_1(1,cnt) = (Pmu_2*Q_nu*(1 - DELTA_Z_nu) - Pnu_2*Q_mu*(1 - DELTA_Z_mu) + (Pnu_2*Pmu_0 - Pmu_2*Pnu_0)*YL_0) / (Pnu_1*Pmu_2 - Pmu_1*Pnu_2);

    %% single-photon error rate
    SF_eXU_1(1,cnt) = (T_nu*(1 + DELTA_PRIME_X_nu) - Pnu_0*e_0*YL_0) / (Pnu_1 * SF_YZL_1(1,cnt));

    %% single-photon key generation rate
    SF_R(1,cnt) = q * (Pmu_1*SF_YZL_1(1,cnt)*(1-binary_shannon_entropy(SF_eXU_1(1,cnt))) - Q_mu*error_correction_rate(E_mu)*binary_shannon_entropy(E_mu));

end
% %%--------------------------------------------------------------------------------------------

%% figures

L1 = 1:1:250;

% Y_1 figure
figure
L11 = plot(L1,Y_1,'.r');
hold on;
L12 = plot(L1,Y_1_inf,'b');
L13 = plot(L1,SF_YZL_1,'g');
xlabel('Distance (km)');
ylabel('Y_1');
% title('Single Photon response rate')

% e_1 figure
figure
L21 = plot(L1,e_1,'.r');
hold on;
L22 = plot(L1,e_1_inf,'b');
L23 = plot(L1,SF_eXU_1,'g');
xlabel('Distance (km)');
ylabel('e_1');
% title('Single Photon error rate')

% R figure
figure
L31 = semilogy(L1,R,'.r');
hold on;
L32 = semilogy(L1,R_inf,'b');
L33 = semilogy(L1,SF_R,'g');
xlabel('Distance (km)');
ylabel('Key generation rate');

% %%--------------------------------------------------------------------------------------------

