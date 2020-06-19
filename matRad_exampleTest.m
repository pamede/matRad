% matRad example script for doing fine sampling in dose calculation and
% comparing it to Monte Carlo simulation
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
% load scLung01small.mat
% clear pln stf

% load patient data, i.e. ct, voi, cst
% load Slab01.mat
% load PHANTOM_slab_entrance_10mm.mat
% load PHANTOM_control.mat
load lung01.mat
clear pln stf resultGUI

% meta information for treatment plan
pln.radiationMode   = 'protons';     
pln.machine         = 'generic_MCsquare';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 200; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 200;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
% load protons_generic_MCsquare
% stf.ray.energy = machine.data(39).energy;

%% dose calculation
 % analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    tic
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    toc
    resultGUI = matRad_calcCubes(ones(sum([stf(:).totalNumOfBixels]),1),dij);
    anaDose     = resultGUI.physicalDose;

 % Monte Carlo dose
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), 1e6);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;

 %% plot doses

matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% ma    tRadGUI
