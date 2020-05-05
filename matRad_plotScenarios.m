
matRad_rc
    
% load scSlab01Result
% load scAlderson01Result
% load scLiver01Result
load scLung01Result

anaSum = sum((anaDose - mcDose).^2,'all');
fsSum = sum((anaFsDose - mcDose).^2,'all');
scSum = sum((anaScDose - mcDose).^2,'all');

figure
X = categorical({'analytical','fine sampling','stdCorr'});
X = reordercats(X,{'analytical','fine sampling','stdCorr'});
Y = [anaSum fsSum scSum];
bar(X,Y)
title('total squared dose deviation');



%% plot dose 

maxDose = max([max(anaDose, [], 'all'),max(mcDose, [], 'all'),max(anaFsDose, [], 'all'),max(anaScDose, [], 'all')]);
propPlot = [1,4];
subplot(propPlot(1),propPlot(2),1);
limits = [35, 115, 29, 109];

figure
imagesc(anaDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)));
caxis([0 maxDose]);
title('Standard dose');
xlabel('x [mm]') 
ylabel('y [mm]') 
hold on 
contour(anaDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)),linspace(0,maxDose,16),'color','black');
hold off
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

subplot(propPlot(1),propPlot(2),2);
imagesc(anaFsDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)));
caxis([0 maxDose]);
title('Fine sampling dose');
xlabel('x [mm]') 
ylabel('y [mm]') 
hold on 
contour(anaFsDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)),linspace(0,maxDose,16),'color','black');
hold off
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

subplot(propPlot(1),propPlot(2),3);
imagesc(anaScDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)));
caxis([0 maxDose]);
title('Std corrected dose');
xlabel('x [mm]') 
ylabel('y [mm]') 
hold on 
contour(anaScDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)),linspace(0,maxDose,16),'color','black');
hold off
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

subplot(propPlot(1),propPlot(2),4);
imagesc(mcDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)));
caxis([0 maxDose]);
title('Monte Carlo dose');
xlabel('x [mm]') 
ylabel('y [mm]') 
hold on 
contour(mcDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z)),linspace(0,maxDose,16),'color','black');
hold off
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

%% execute gamma tests
figure
propPlot = [1,3];
gammaTest = [2, 2];
interpolation = 1;

subplot(propPlot(1),propPlot(2),1);
[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z) ,interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'stadard dose'});
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

subplot(propPlot(1),propPlot(2),2);
[gammaCube2,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corrected dose'});
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])

subplot(propPlot(1),propPlot(2),3);
[gammaCube4,gammaPassRateCell] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'fine sampling dose'});
axis(limits)
pbaspect([ct.cubeDim(2), ct.cubeDim(1), ct.cubeDim(3)])