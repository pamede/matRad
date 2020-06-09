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
% load Lung01
% load scLung02
% load scProstate01.mat


pln.machine = 'generic_MCsquare';





stf = stf(1);





%% dose calculation
 % analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    resultGUI.physicalDose_FS = matRad_calcParticleDose(ct,stf,pln,cst);
    anaDose     = resultGUI.physicalDose;
    
% Monte Carlo dose
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights, 1e6);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;
    
    matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');

 % analytical dose with fine sampling
    pln.propDoseCalc.anaMode = 'fineSampling';
    pln.propDoseCalc.fineSampling.N = 21;
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(weights,dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose     = resultGUI.physicalDoseFS;
    
 % analytical dose with stdCorr
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(weights,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose     = resultGUI.physicalDoseSC;

 

 %% plot doses

matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaScDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaFsDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');

