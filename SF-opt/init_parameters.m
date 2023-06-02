function [mu_init, nu_init] = init_parameters(L)

    z = zeros(1000);

    for i = 1:1000
        for j = 1:(i-1) % nu < mu
            mu = i / 1000
            nu = j / 1000;

            z(i,j) = real(keyrate(mu, nu, L));
            
            if (z(i,j) < 0)
                z(i,j) = 0;
            end
        end
    end

    [maxZ, maxIndex] = max(z(:));
    [maxRowIndex, maxColIndex] = ind2sub(size(z), maxIndex);

    mu_init = maxRowIndex / 1000;
    nu_init = maxColIndex / 1000;

end

%% Result

% maxZ =

%     0.0473


% maxIndex =

%         7716


% maxRowIndex =

%    716


% maxColIndex =

%      8