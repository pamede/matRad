function [weights, xPos] = matRad_calcWeights2D(sigma, N, sigmaSub)
% Russo splitting

sigmaDach = sqrt(sigma^2 - sigmaSub^2);
deltaR = 7 * sigmaDach / N;

weight = @(mu, sigma, i, deltaR) integral(@(x) normpdf(x, mu, sigma), mu + (i - 0.5) * deltaR, mu + (i + 0.5) * deltaR);

weights = zeros(N,1);
xPos = zeros(N,1);
count = 1;
for i = -(N - 1) / 2 : (N - 1) / 2
    weights(count)  = weight(0, sigmaDach, i, deltaR);
    xPos(count) = i * deltaR;
    count = count + 1;
end

