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
% load Alderson01
% load scLiver01
% load scLung01
% load scLung02
% load scProstate01.mat

%% dose calculation
 % analytical dose without fine sampling
%     pln.propDoseCalc.anaMode = 'standard';
%     dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
%     resultGUI = matRad_calcCubes(weights,dij);
%     anaDose     = resultGUI.physicalDose;
%     
%  % analytical dose with fine sampling
%     pln.propDoseCalc.anaMode = 'fineSampling';
%     pln.propDoseCalc.fineSampling.N = 23;
%     dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
%     resultGUI_FS = matRad_calcCubes(weights,dijFS);
%     resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
%     anaFsDose     = resultGUI.physicalDoseFS;
%     
 % analytical dose with stdCorr
    pln.propDoseCalc.anaMode = 'stdCorr';
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(weights,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose     = resultGUI.physicalDoseSC;

%  % Monte Carlo dose
%     resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,weights, 1e6);
%     resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
%     mcDose      = resultGUI.physicalDoseMC;

 %% plot doses
 
maxDose = max(max(anaScDose, [], 'all'));

figure
imagesc(anaScDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))');
caxis([0 maxDose]);
hold on 
contour(anaScDose(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))',linspace(0,maxDose,16),'color','black');
contour(ct.cube{1}(:,:,round(stf(1).isoCenter(3)/ct.resolution.z))',linspace(0,1,1),'color','white','linewidth',3);

hold off
pbaspect([ct.cubeDim(1), ct.cubeDim(2), ct.cubeDim(3)])
xlabel('x [mm]') 
ylabel('y [mm]') 

%%
Gauss    = @(x,mu,SqSigma) 1./(sqrt(2*pi.*SqSigma)).*exp(-((x - mu).^2./(2.*SqSigma)));

sumGauss = @(x,mu,SqSigma,w) (1./sqrt(2*pi*ones(numel(x),1) .* SqSigma') .* ...
    exp(-bsxfun(@minus,x,mu').^2 ./ (2* ones(numel(x),1) .* SqSigma' ))) * w;


load protons_generic_MCsquare.mat
baseData = machine.data(70);
% plot(machine.data(end).depths,machine.data(end).Z.profileORG,'linewidth', 1.5);
plot(machine.data(end).depths, sumGauss(machine.data(end).depths,baseData.Z.mean,baseData.Z.width.^2,baseData.Z.weight),'linewidth', 1.2);
hold on
for i = 1:2:20
    plot(machine.data(end).depths, sumGauss(machine.data(end).depths,baseData.Z.mean,baseData.Z.width.^2+i^2,baseData.Z.weight),'linewidth', 1.2);
end
xlabel('Depth [mm]') 
ylabel('Dose [a.u.]')
xlim([0 300]);
ylim([0 22]);


% matRad_compareDose(anaDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaScDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');
% matRad_compareDose(anaFsDose, mcDose, ct, cst, [1, 1, 0] , 'on', pln, [2,2], 1, 'global');

