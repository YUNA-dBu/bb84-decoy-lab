function Q_x = overall_qubit_gain(lambda, Y_0, eta)
    Q_x = Y_0 + 1 - exp(-eta*lambda);
    % Q_x = 0;
    % for n=0:1:100
    %     P_n = poisspdf(n, lambda);
    %     Y_n = yield(n, Y_0, eta);
    %     Q_x = Q_x + P_n * Y_n;
    % end
end
