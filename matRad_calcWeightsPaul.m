function [weights, xPos, yPos] = matRad_calcWeightsPaul(sigma, N, sigmaSub)
% Russo splitting

sigmaDach = sqrt(sigma^2 - sigmaSub^2);
deltaR = 7 * sigmaDach / N;

weight = @(mu, sigma, i, deltaR) integral(@(x) normpdf(x, mu, sigma), mu + (i - 0.5) * deltaR, mu + (i + 0.5) * deltaR);

gri = zeros(N,1);
coord = zeros(N,1);
count = 1;
for i = -(N - 1) / 2 : (N - 1) / 2
    gri(count)  = weight(0, sigmaDach, i, deltaR);
    coord(count) = i * deltaR;
    count = count + 1;
end

weightGrid = gri * gri';

weights = zeros(N^2,1);
xPos    = zeros(N^2,1);
yPos    = zeros(N^2,1);

count = 1;
for i = 1:N
    for j = 1:N
        weights(count) = weightGrid(i,j);
        xPos(count) = coord(i);
        yPos(count) = coord(j);
        count = count + 1;
    end
end