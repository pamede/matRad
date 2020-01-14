 %% example fine smapling 2 dimensions

 clear all, close all
 load slabPhantom1.mat
 load protons_generic_MCsquare.mat
 
 baseData = machine.data(30);
 baseData.LatCutOff.CompFac = 0.999;
 
 SSD = 2000;
 initSigma = matRad_interp1(baseData.initFocus.dist(1,:)', baseData.initFocus.sigma(1,:)',SSD);
 
 %% ray casting
 
sliceCt = ct.cube{1}(:,:,80);
fsDose = zeros(160,160);
 % split gaussian
sigmaSub = 2;
N = 11;
[weights, xPos] = matRad_calcWeights2D(initSigma, N, sigmaSub);





for w = 1:numel(weights)
    weightIx = round(xPos(w) / ct.resolution.x);
    singleRadDepthCast = cumsum(sliceCt(:,80 + weightIx));

    radDepths = [];
    for i = 1:160
        radDepths = [radDepths, singleRadDepthCast];
    end
    radDepths = radDepths * ct.resolution.y;

    latDistCast = abs(((1:160) - 80)* ct.resolution.x + xPos(w));
    latDist = [];
    for i = 1:160
        latDist = [latDist; latDistCast];
    end
    radDepths = reshape(radDepths,160^2,1);
    latDist = reshape(latDist,160^2,1);

    dose = matRad_calcParticleDoseBixel(radDepths, latDist.^2, sigmaSub, baseData);
    fsDose = fsDose + reshape(dose,160,160) * weights(w);
end

plot(-40:0.01:40, normpdf(-40:0.01:40, 0, initSigma));
hold on
sum = zeros(size(-40:0.01:40));
for i = 1:N
    plot(-40:0.01:40, normpdf(-40:0.01:40, xPos(i), sigmaSub) * weights(i));
    sum = sum + normpdf(-40:0.01:40, xPos(i), sigmaSub) * weights(i);
end
plot(-40:0.01:40, sum);
hold off

figure
imagesc(flip(fsDose(40:120,40:120),2));


% imagesc(radDepth);
% figure
% imagesc(latDist);
 
 %% dose calculation