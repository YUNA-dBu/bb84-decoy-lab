function e_n = photon_error_rate(n, Y_0, Y_n)

    %% misalignment error rate 
    e_d = 0.02 
    %% symbol error rate of vacuum pulse
    e_0 = 0.5

    e_n = (e_0*Y_0 + e_d*(Y_n-Y_0)) / Y_n;

end
