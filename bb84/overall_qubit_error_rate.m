function T_x = overall_qubit_error_rate(lambda, Y_0, eta)
    %% misalignment error rate 
    e_d = 0.02;
    %% symbol error rate of vacuum pulse
    e_0 = 0.5;

    T_x = e_0 * Y_0 + e_d * (1-exp(-eta*lambda));

    % T_x = 0;
    % for n=0:1:100
    %     P_n = poisspdf(n, lambda);
    %     Y_n = yield(n, Y_0, eta);
    %     e_n = photon_error_rate(n, Y_0, Y_n);
    %     T_x = T_x + P_n * Y_n * e_n;
    % end
end
