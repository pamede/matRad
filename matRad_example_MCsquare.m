% matRad script
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

%load HEAD_AND_NECK
% load TG119.mat
%load PROSTATE.mat
%load LIVER.mat
%load BOXPHANTOM
load BOXPHANTOM.mat
% load BOXPHANTOM_NARROW_NEW.mat
% load phantomTest.mat


% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
pln.machine         = 'matRadBDL';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 3; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 3;
pln.propStf.gantryAngles    = 0; % [?] 
pln.propStf.couchAngles     = 0; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% pln.propStf.isoCenter       = [0,0,0];

                            
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);


%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    %dij = matRad_calcPhotonDoseVmc(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,10000);
%     resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1),1000);
   
end

resultGUI = matRad_calcCubes(ones(dij.totalNumOfBixels,1),dij);
%resultGUI_MC = matRad_calcCubes(resultGUI.w,dijMC);
resultGUI.physicalDose_MC = resultGUI_MC.physicalDose;
resultGUI.physicalDose_diff = (resultGUI.physicalDose - resultGUI.physicalDose_MC);

mcDose = reshape(resultGUI.physicalDose_MC, ct.cubeDim);
anaDose = reshape(resultGUI.physicalDose, ct.cubeDim);

[gammaCube,gammaPassRateCell] = matRad_gammaIndex(mcDose,anaDose,ct.cubeDim,[2,2],round(ct.cubeDim(3)/2),0,'global',cst);


% mcIDD = sum(sum(resultGUI.physicalDose_MC,2),3);
% 
% plot(mcIDD);
% hold on
% 
% anaIDD = sum(sum(resultGUI.physicalDose,2),3);
% anaIDD = anaIDD * 0.32^2;
% 
% plot(anaIDD);
% hold off
