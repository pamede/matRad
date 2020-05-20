% matRad example script for running scenarios
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_rc

% load patient data, i.e. ct, voi, cst
load scSlab01
% load scAlderson01
% load scLiver01
% load scLung01
% load scLung02
% load scProstate01.mat

%% dose calculation
 % analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(weights,dij);
    anaDose     = resultGUI.physicalDose;
    
%  % analytical dose with fine sampling
%     pln.propDoseCalc.anaMode = 'fineSampling';
%     pln.propDoseCalc.fineSampling.N = 11;
%     dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
%     resultGUI_FS = matRad_calcCubes(weights,dijFS);
%     resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
%     anaFsDose     = resultGUI.physicalDoseFS;
    
 % analytical dose with stdCorr
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(weights,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose     = resultGUI.physicalDoseSC;

  % analytical dose with std convolution approach
    load protons_generic_MCsquare.mat
    weight = arrayfun(@(d) (d/machine.data(38).peakPos).^5, meanRadDepths);
    weight(weight < 0) = 0;
    cStdCtGrid = weight .* cStdCtGrid;
    anaCvDose = matRad_sliceConvnFilter(ct, anaDose,cStdCtGrid);
    resultGUI.physicalDoseConv = anaCvDose; 

 % Monte Carlo dose
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights, 1e6);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;
    
    
maxDose = max(max(anaCvDose, [], 'all'));

figure
imagesc(anaCvDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))');
caxis([0 maxDose]);
hold on 
contour(anaCvDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))',linspace(0,maxDose,16),'color','black');
contour(ct.cube{1}(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))',linspace(0,1,1),'color','white','linewidth',3);

hold off
pbaspect([ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3)])
xlabel('x [mm]') 
ylabel('y [mm]') 


% weight = arrayfun(@(d) interp1(machine.data(36).depths, (machine.data(36).Z.profileORG-machine.data(36).Z.profileORG(1))/max(machine.data(36).Z.profileORG),d), meanRadDepths);
% weight = arrayfun(@(d) (d/machine.data(36).peakPos).^5, meanRadDepths);
% weight(weight < 0) = 0;
% newCstd1 = weight .* cStdCtGrid;
% newCstd = imgaussfilt3(newCstd1,1);
% tmpCstd = sqrt(newCstd1);
% tmpCstd = imgaussfilt3(tmpCstd,4);
% tmpCstd = newCstd1;






 %% plot doses

% matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaScDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaFsDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaCvDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');

