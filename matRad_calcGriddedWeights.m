function finalWeight = matRad_calcGriddedWeights(mu, gridSize, gridX, gridY, sigmaTot, sigmaSub)
                    
sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);

gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

% finalWeight = arrayfun(@(x,y) integral2(@(x,y) gauss(sigmaHead, x, y, mu(1), mu(2)),... % Platz 2
%                     x - (0.5 * gridSize(1)), x + (0.5 * gridSize(1)),...
%                     y - (0.5 * gridSize(2)), y + (0.5 * gridSize(2))), gridX, gridY);

finalWeight = zeros(size(gridX));         % Platz 1
parfor i = 1:numel(gridX)
        finalWeight(i) = integral2(@(x,y) gauss(sigmaHead, x, y, mu(1), mu(2)), ...
                                            gridX(i) - (0.5 * gridSize(1)), gridX(i) + (0.5 * gridSize(1)), ...
                                            gridY(i) - (0.5 * gridSize(2)), gridY(i) + (0.5 * gridSize(2)));
end

% finalWeight = zeros(size(gridX));         % Platz 3
% for i = 1:numel(gridX)
% finalWeight(i) = quad2d(@(x,y) gauss(sigmaHead, x, y, mu(1), mu(2)), ...
%                                             gridX(i) - (0.5 * gridSize(1)), gridX(i) + (0.5 * gridSize(1)), ...
%                                             gridY(i) - (0.5 * gridSize(2)), gridY(i) + (0.5 * gridSize(2)));
% end                                      
