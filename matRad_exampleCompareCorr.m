
matRad_rc

% load patient data, i.e. ct, voi, cst
% load lung01MC.mat
load slab01MC.mat

pln.machine = 'matRadBDL_APM';

mcDose = mcDose * 1.1252;

%% dose calculation
% analytical dose
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(resultGUI.w,dij);
    anaDose   = resultGUI.physicalDose;
        
% analytical dose with fine sampling
    pln.propDoseCalc.anaMode = 'fineSampling';
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 21;
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(resultGUI.w,dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI.physicalDoseFS;
    
% analytical dose with std correction
%     pln.machine = 'matRadBDL_APM';
    pln.propDoseCalc.anaMode = 'stdCorr';    
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(resultGUI.w,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose   = resultGUI.physicalDoseSC;
    
%% execute gamma tests
gammaTest = [2, 2];
interpolation = 1;

[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'stadard dose'});
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),1,'color','white');
hold off

[gammaCube2,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corrected dose'});
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),1,'color','white');
hold off

[gammaCube4,gammaPassRateCell] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(ct.cubeDim(3)/2),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'fine sampling dose'});
hold on
contour(ct.cube{1}(:,:,round(ct.cubeDim(3)/2)),1,'color','white');
hold off