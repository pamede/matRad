clear
close all

gridSize = [2,2];
stf.ray.rayPos_bev = [0, 0, 0];
[gridX, gridY] = matRad_createFineSamplingGrid(stf, gridSize);
% [gridX, gridY] = meshgrid(-5:gridSize(1):5, -5:gridSize(1):5);

mu = [0, 0];
sigmaTot = 5.1832;
thresh = 0;
sigmaSub = 1;

[ixWeight, finalWeight] = matRad_calcGriddedWeights(mu, gridSize, gridX, gridY, sigmaTot, thresh, sigmaSub);


[plotGridX, plotGridY] = meshgrid(-20:1:20, -20:1:20);

gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

erg = gauss(sigmaTot, plotGridX, plotGridY, 0, 0);

surf(erg)

sum = zeros(size(plotGridX));   
for i = 1:numel(gridX)
    tmp = finalWeight(i) * gauss(sigmaSub, plotGridX, plotGridY, gridX(i), gridY(i));
    sum = sum + tmp;
end
figure
surf(sum)

% surf(erg - sum)