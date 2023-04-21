function Q_x = overall_qubit_gain(lambda, Y_0, eta)
    Q_x = 0;
    for n=0:1:10
        P_n = poisspdf(n, lambda);
        Y_n = yield(n, Y_0, eta);
        T_x = T_x + P_n * Y_n;
    end
end
