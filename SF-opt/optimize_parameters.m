function [mu, nu, val] = optimize_parameter(mu_init, nu_init)
    %% init optimal value
    mu = zeros(251, 1);
    nu = zeros(251, 1);
    val= zeros(251, 1);     % optimized R VALue
    mu(1) = mu_init;        % mu(1) = 0.716;
    nu(1) = nu_init;        % nu(1) = 0.008;

    %% loop over the distance
    for L = 2:1:251         % distance
            
        piv    = 99999;     % pivot, cache the R of previous optimization loop
        val(L) = 0;         % init the maximum with 0
        mu(L)  = mu(L-1);   % inherit the continuous parameter from previous loop
        nu(L)  = nu(L-1);   % same

        % when piv is almost the same as val, 
        % the program should have reached the optimized (mu, nu) pair for the distance. 
        while(abs(piv - val(L)) > 10^(-4))
            piv = val(L);       % cache the previous val
            ptr = 99999;        % pointer of mu's increment vector(d_i) No. Set > 0 so that the loop can start
            % fix nu, optimize mu
            while(ptr > 0) 
                z_i = zeros(4, 1);      % cache the value of 4 increment candidate pairs
                d_i = [-1  1  0  0];    % increment vector
                
                % flag the temporary maximum value
                % init var with current value as default
                % use real() to get the real part, to make MATLAB happy in some cases like L > 100
                ptr = 0;                % suppose not found
                val(L) = real(keyrate(mu(L), nu(L), L)); 

                % iterate the vectors, find a larger value
                % the increment vector is 1000x scaled up, so divide it by 1000
                for k = 1:1:2
                    z_i(k) = real(keyrate(mu(L) + d_i(k)/1000, nu(L), L));
                    if(z_i(k) > val(L))
                        ptr    = k;
                        val(L) = z_i(k);
                    end
                end

                % (ptr > 0) <--> found
                if (ptr > 0)
                    mu(L) = mu(L) + d_i(ptr)/1000;
                end
            end
            
            ptr = 99999;        % pointer of nu's increment vector(d_j) No. Set > 0 so that the loop can start
            % fix mu, optimize nu
            while(ptr > 0) 
                z_j = zeros(4, 1);      % cache the value of 4 increment candidate pair
                d_j = [ 0  0 -1  1];    % increment vector
                
                % flag the temporary maximum value
                % init var with current value as default
                % use real() to get the real part, to make MATLAB happy in some cases like L > 100
                ptr = 0;                % suppose not found
                val(L) = real(keyrate(mu(L), nu(L), L)); 

                % iterate the vectors, find a larger value
                % the increment vector is 1000x scaled up, so divide it by 1000
                for k=3:1:4
                    z_j(k) = real(keyrate(mu(L), nu(L) + d_j(k)/1000, L));
                    if(z_j(k) > val(L))
                        ptr    = k;
                        val(L) = z_j(k);
                    end
                end

                % (ptr > 0) <--> found
                if (ptr > 0)
                    nu(L) = nu(L) + d_j(ptr)/1000;
                end
            end 
        end
    end
end