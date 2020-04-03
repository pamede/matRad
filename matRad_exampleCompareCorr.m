
matRad_rc

% load patient data, i.e. ct, voi, cst
% load lung01MC.mat
load slab01MC.mat


%% dose calculation
% analytical dose with fine sampling
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 21;
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(resultGUI.w,dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI.physicalDoseFS;

% analytical dose with std correction
pln.propDoseCalc.stdCorr = true;
