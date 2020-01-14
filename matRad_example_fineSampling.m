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
load slabPhantom1.mat
% load PROSTATE.mat
% load LIVER.mat
% load BOXPHANTOM
% load BOXPHANTOMv3.mat
% load BOXPHANTOM_NARROW_NEW.mat
% load phantomTest.mat


% meta information for treatment plan
pln.radiationMode   = 'protons';     % either photons / protons / carbon
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
pln.propDoseCalc.doseGrid.resolution.x = ct.resolution.x; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = ct.resolution.y; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = ct.resolution.z; % [mm]
pln.propDoseCalc.airOffsetCorrection = true;

% optimization settings
pln.propOpt.optimizer       = 'IPOPT';
pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                      % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
stf = matRad_computeSSD(stf,ct);
load protons_generic_MCsquare;
stf.ray.energy = machine.data(30).energy;

% stf = matRad_splitBeam(stf, machine);


% 
%  % calculate dose
%  figure
% for i = 1:stf.numOfRays
%     [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(stf.isoCenter,ct.resolution,...
%                                     stf.ray(i).rayPos,stf.ray(i).targetPoint,ct.cube);
%      plot(cumsum(rho{1}));
%      display(stf.ray(i).targetPoint(2));
%      hold on   
% end
% hold off
% dose = matRad_calcParticleDoseBixel(radDepths, radialDist_sq, sigmaIni_sq, baseData);                                

%% dose calculation
   dij = matRad_calcParticleDose(ct,stf,pln,cst);
%     dijMC = matRad_calcParticleDoseMC(ct,stf,pln,cst,1000000);
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum(stf(:).totalNumOfBixels),1),1000000);

% w  = [stf.ray(:).weight]';
w = ones(size(stf.numOfRays));
% w(1) = 1;
resultGUI = matRad_calcCubes(w,dij);

% 
% stf = matRad_splitBeam(stf, machine);
% dij = matRad_calcParticleDose(ct,stf,pln,cst);
% w  = [stf.ray(:).weight]';
% resultGUI_FS = matRad_calcCubes(w,dij);



imagesc(resultGUI.physicalDose(40:120,40:120,80));
figure
imagesc(resultGUI_MC.physicalDose(40:120,40:120,80));

