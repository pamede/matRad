function [ixWeight, finalWeight] = matRad_calcGriddedWeights(mu, gridSize, gridX, gridY, sigmaTot, sigmaSub)
                    
sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);

gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

ixWeight = [];
finalWeight = [];
for i = 1:numel(gridX)
        weight = integral2(@(x,y) gauss(sigmaHead, x, y, mu(1), mu(2)), ...
                                            gridX(i) - (0.5 * gridSize(1)), gridX(i) + (0.5 * gridSize(1)), ...
                                            gridY(i) - (0.5 * gridSize(2)), gridY(i) + (0.5 * gridSize(2)));
%         if weight > thresh
            ixWeight = [ixWeight, i];
            finalWeight = [finalWeight, weight];
%         end

% %%%%
% intfun = @(x,y) gauss(sigmaHead, x, y, mu(1), mu(2));
% 
% weight = arryfun(@(x,y) integral2(@intfun, x - (0.5 * gridSize(1)), x + (0.5 * gridSize(1)), y - (0.5 * gridSize(2)), y + (0.5 * gridSize(2)),...
%                         gridX, gridY);

end