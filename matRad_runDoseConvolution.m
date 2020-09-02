matRad_rc

load SCENARIO_slab_entrance_10mm.mat

%% dose calculation
 % analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(weights,dij);
    anaDose     = resultGUI.physicalDose;
    
 % analytical dose with stdCorr
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false,[1 2.1]);
    resultGUI_SC = matRad_calcCubes(weights,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose     = resultGUI.physicalDoseSC;

  % analytical dose with std convolution approach
    load protons_generic_MCsquare.mat
    weight = arrayfun(@(d) (d/machine.data(38).peakPos).^2, meanRadDepths);
    weight(weight < 0) = 0;
    cStdCtGrid = weight .* cStdCtGrid;
    anaCvDose = matRad_sliceConvnFilter(ct, anaDose,cStdCtGrid);
    resultGUI.physicalDoseConv = anaCvDose; 
   
    
slice = 25;
cMax = max([max(anaDose(:)) max(anaScDose(:))  max(mcDose(:))]);
levels = ((1:8) - 0.5) / 8 * cMax;

figure
subaxis(1,1,1,'MarginTop',0,'MarginBottom',0,'MarginLeft',0,'MarginRight',0);
imagesc(anaDose(:,:,slice))
hold on
contour(anaDose(:,:,slice),levels,'LineWidth',1,'Color', 'k')
colorbar
caxis([0 cMax]);
contour(ct.cubeHU{1}(:,:,25),1,'LineWidth',2,'Color', 'w')
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])
set(gcf,'Position',[100 200 215 480])

figure
subaxis(1,1,1,'MarginTop',0,'MarginBottom',0,'MarginLeft',0,'MarginRight',0);
imagesc(anaCvDose(:,:,slice))
hold on
contour(anaCvDose(:,:,slice),levels,'LineWidth',1,'Color', 'k')
colorbar
caxis([0 cMax]);
contour(ct.cubeHU{1}(:,:,25),1,'LineWidth',2,'Color', 'w')
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])
set(gcf,'Position',[100 200 215 480])

figure
subaxis(1,1,1,'MarginTop',0,'MarginBottom',0,'MarginLeft',0,'MarginRight',0);
imagesc(anaScDose(:,:,slice))
hold on
contour(anaScDose(:,:,slice),levels,'LineWidth',1,'Color', 'k')
colorbar
caxis([0 cMax]);
contour(ct.cubeHU{1}(:,:,25),1,'LineWidth',2,'Color', 'w')
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])
set(gcf,'Position',[100 200 215 480])

figure
subaxis(1,1,1,'MarginTop',0,'MarginBottom',0,'MarginLeft',0,'MarginRight',0);
imagesc(mcDose(:,:,slice))
hold on
contour(mcDose(:,:,slice),levels,'LineWidth',1,'Color', 'k')
colorbar
caxis([0 cMax]);
contour(ct.cubeHU{1}(:,:,25),1,'LineWidth',2,'Color', 'w')
pbaspect([ct.cubeDim(2) ct.cubeDim(1) ct.cubeDim(3)])
set(gcf,'Position',[100 200 215 480])
