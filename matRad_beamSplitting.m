clear all, close all
weight = @(mu, sigma, i, deltaR) integral(@(x) normpdf(x, mu, sigma), mu + (i - 0.5) * deltaR, mu + (i + 0.5) * deltaR);

% plot original gauss
sigma = 10;
mu = 0;

% plot splitted gaussians
N = 21;
sigmaSub = 2;
sigmaDach = sqrt(sigma^2 - sigmaSub^2);
deltaR = 7 * sigmaDach / N;

gri = zeros(N,1);
xVal = zeros(N,1);
count = 1;
for i = -(N - 1) / 2 : (N - 1) / 2
    gri(count)  = weight(mu, sigmaDach, i, deltaR);
    xVal(count) = mu + i * deltaR;
    count = count + 1;
end

yVal = xVal;
weightGrid = gri * gri';

M = 20;
length = 80;
for i = 1:M
    for j = 1:M
       coord(i,j,:) = [i*length/M - length/2, j*length/M - length/2];
    end
end

gauss2dCoord = @(coord, mu, sigma) 1 / (2 * pi * sigma^2) .* exp(-(coord(:,:,1) - mu(1)).^2 / (2 * sigma^2)) .* exp(-(coord(:,:,2) - mu(2)).^2 / (2 * sigma^2));
        
sum = zeros(M,M);
hold on
for i = 1:N
    for j = 1:N
        erg = gauss2dCoord(coord, [xVal(i), yVal(j)], sigmaSub) * weightGrid(i,j);
        surf(erg);
        sum = sum + erg;
    end
end
hold off

figure 
surf(sum - gauss2dCoord(coord, [0, 0], sigma));
figure
surf(sum);