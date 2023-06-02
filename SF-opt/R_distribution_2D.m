maxRowIndex = zeros(1, 251);
maxColIndex = zeros(1, 251);
maxZ = zeros(1, 251);

figure

for L = 1:10:251
    z = zeros(100, 100);
    L
    for i = 1:100
        for j = 1:i-1
                mu = i / 100;
                nu = j / 100;
                z(i,j) = real(keyrate(mu, nu, L));
                if (z(i,j) < 0)
                    z(i,j) = 0;
                end
        end
    end
    
    [maxZ(L), maxIndex] = max(z(:));
    [maxRowIndex(L), maxColIndex(L)] = ind2sub(size(z), maxIndex);
    
    x = 1:100;
    y = 1:100;
    
    % figure;
    axis tight manual;
    ax = gca;
    
    imagesc(x,y,z);
    colorbar;
    
    xlabel('nu');
    ylabel('mu');
    title(['Distance(km) ', num2str(L)]);
    
    hold on;
    plot(maxColIndex(L), maxRowIndex(L), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); 
    hold off;
    
    drawnow;
    pause(0.1);
    % clf(ax);
end

L = 1:251;

figure
hold on
plot(L,maxRowIndex,'.')
plot(L,maxColIndex,'.')
hold off

figure
semilogy(L,maxZ,'.')