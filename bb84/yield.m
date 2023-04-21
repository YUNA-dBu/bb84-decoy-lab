function Y_n = yield(n, Y_0, eta)
    Y_n = 1 - (1-Y_0)*power(1-eta, n);
end	
