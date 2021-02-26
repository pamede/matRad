% matRad example script 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team. 
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
load TG119.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon

% define radiation machine through used matRad base data file
pln.machine          = 'generic_MCsquare';

pln.numOfFractions  = 1;

% beam geometry settings
pln.propStf.bixelWidth              = 10; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 10;
pln.propStf.gantryAngles            = [0, 90]; % [?] 
pln.propStf.couchAngles             = [0, 90]; % [?]
pln.propStf.numOfBeams              = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter               = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 5; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 5; % [mm]
%pln.propDoseCalc.doseGrid.resolution = ct.resolution;

%Activate/Deactive the minimal correction for air WEPL given the SSD 
%differences between base data fit and treatment plan
%pln.propDoseCalc.airOffsetCorrection = false;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

% pln.propMC.proton_engine    = 'MCsquare';
pln.propMC.proton_engine    = 'TOPAS'; %Requires separate topas installation
pln.propMC.uncQuantParam    = [0.01, 3.5];
% uncQuantParam(1): relative energy error/spread (relative to mean energy)
% uncQuantParam(2): absolut spot size error (added in quadrature) in [mm]
% uncQuantParam(3): beam divergence  (optional, default = 0) in [radian]
% uncQuantParam(4): beam correlation (optional, default = 0)


%Enable LET calculation - Does not work yet for MCsquare data
pln.propDoseCalc.calcLET    = false;

%Enable/Disable use of range shifter (has effect only when we need to fill 
%up the low-range region)
pln.propStf.useRangeShifter = true;  

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);

%% dose calculation

dij = matRad_calcParticleDose(ct, stf, pln, cst); %Calculate particle dose influence matrix (dij) with analytical algorithm
%dij = matRad_calcParticleDoseMC(ct,stf,pln,cst,1e4); %Calculate particle dose influence matrix (dij) with MC algorithm (slow!!)


resultGUI = matRad_fluenceOptimization(dij,cst,pln); %Optimize
%resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij); %Use uniform weights


%% MC calculation
%resultGUI_recalc = matRad_calcDoseDirect(ct,stf,pln,cst,resultGUI.w);       %Recalculate particle dose analytically
resultGUI_recalc = matRad_calcDoseDirectMC(ct,stf,pln,cst,resultGUI.w,1000);  %Recalculate particle dose with MC algorithm
resultGUI.physicalDose_MC = resultGUI_recalc.physicalDose;



