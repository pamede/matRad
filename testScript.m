gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 


[plotX, plotY] = meshgrid(-20:20,-20:20);

gridSize = [0.5, 0.5];

[gridX, gridY] = meshgrid(-20:gridSize(1):20, -20:gridSize(1):20);
gridX = reshape(gridX, numel(gridX), 1);
gridY = reshape(gridY, numel(gridY), 1);

mu = [3, 3];
sigmaTot = 5;
sigmaSub = 1;
thresh = 0;

[ixWeight, finalWeight] = matRad_calcGriddedWeights(mu, gridSize, gridX, gridY, sigmaTot, thresh, sigmaSub);

sum = zeros(size(plotX));   
for i = 1:numel(gridX)
    tmp = finalWeight(i) * gauss(sigmaSub, plotX, plotY, -gridX(i), -gridY(i));
    sum = sum + tmp;
    surf(tmp)
%     waitforbuttonpress;
end

surf(plotX, plotY, sum);
figure
surf(plotX, plotY, gauss(sigmaTot, plotX, plotY, mu(1),mu(2)))
figure
surf(plotX, plotY, gauss(sigmaTot, plotX, plotY, mu(1),mu(2)) - sum)


load scAlderson01

tmp = stf(1).ray(:);
for i = 1:size(tmp, 1)
    x(i) = tmp(i).rayPos_bev(1);
    y(i) = tmp(i).rayPos_bev(3);
end

xLim = [min(x), max(x)];
yLim = [min(y), max(y)];
    
    
