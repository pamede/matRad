clear

%% one dimension
weight = @(mu, sigma, i, deltaR) integral(@(x) normpdf(x, mu, sigma), mu + (i - 0.5) * deltaR, mu + (i + 0.5) * deltaR);

% plot original gauss
sigma = 10;
mu = 0;
x = -40:0.1:40;
plot(x, normpdf(x,mu,sigma));
hold on

% plot splitted gaussians
N = 21;
sigmaSub = 1;
sigmaDach = sqrt(sigma^2 - sigmaSub^2);
deltaR = 7 * sigmaDach / N;
sum = zeros(size(x,2),1);
for i = -(N - 1) / 2 : (N - 1) / 2
    beamletGauss = @(x) weight(mu, sigmaDach, i, deltaR) * normpdf(x, mu + i * deltaR, sigmaSub);
    plot(x, beamletGauss(x));
    sum = sum + beamletGauss(x)';
end
hold off

figure
plot(x, sum)
hold on
plot(x, normpdf(x,mu,sigma));
plot(x, abs(sum - normpdf(x,mu,sigma)'));
hold off
%% two dimensions

gauss2d = @(x, y, mu, sigma) 1 / (2 * pi * sigma^2) * exp(-(sqrt(x^2 + y^2) - mu)^2 / (2 * sigma^2));

for m = 1:M
   tmpX = linspace(-40,40,1000);
   tmpY = 
        

surf(X, Y, gauss2d(X, Y, 0, 10));




















% radialGauss = @(r, mu, sigma) 1 / sqrt(2 * pi * sigma^2) .* exp(-(r - mu).^2 /(2 * sigma^2));
% 
% nTheta = 8;
% nR = 3;
% fn = 1;
% sigma = 5;
% rMax = 3 * sigma;
% dR  = rMax / nR;
% dTheta = 360 / nTheta;
% 
% 
% theta = 0;
% i=1;
% r = 0;
% sigmaSub = 2.1;
% F = fn * dTheta / (2 * pi) *(exp(-(0)^2/(2*sigmaSub^2)) - exp(-(dR/2)^2/(2*sigmaSub^2)));
% erg(i,:) = [r, theta, F];
% i = 2;
% for iTheta = 1:nTheta
%     r = dR;
%     for iR = 1:nR
%         F = fn * dTheta / (2 * pi) *(exp(-(r - dR/2)^2/(2*sigmaSub^2)) - exp(-(r + dR/2)^2/(2*sigmaSub^2)));
%         erg(i,:) = [r, theta, F];
%         r = r + dR;
%         i = i + 1;
%     end
%     theta = theta + dTheta;
% end
% format short; 
% erg