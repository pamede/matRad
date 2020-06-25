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

clear
matRad_rc

% load patient data, i.e. ct, voi, cst
% load scLung01small.mat
load alderson01.mat
% load PHANTOM_control.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     
pln.machine         = 'generic_MCsquare';

stf = stf(1);

%% dose calculation
 % analytical dose without fine sampling
    tic
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    toc
    resultGUI = matRad_calcCubes(ones(sum([stf(:).totalNumOfBixels]),1),dij);
    anaDose     = resultGUI.physicalDose;

 % Monte Carlo dose
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), 1e7);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;

 %% plot doses

matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRadGUI