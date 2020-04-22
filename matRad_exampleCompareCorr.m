
matRad_rc

% load patient data, i.e. ct, voi, cst
% load alderson01MC.mat
% load liver01MC.mat
load lung01MC.mat 
% load lung02MC.mat
% load prostate01MC.mat

clear anaDose mcDose resultGUI

num = 5;
beamNum = 1;
stf = stf(beamNum);
nRay = [];
for i = 1:num
    nRay = [nRay, round(rand * size(stf.ray, 2))];
end
nRay = sort(nRay);
stf.numOfRays = num;
stf.ray = stf.ray(nRay);

for i = 1:num
    ix = floor(rand * size(stf.ray(i).energy,2)) + 1;
    stf.ray(i).energy = stf.ray(i).energy(ix);
    stf.ray(i).focusIx = stf.ray(i).focusIx(ix);
    stf.ray(i).rangeShifter = stf.ray(i).rangeShifter(ix);
end

stf.numOfBixelsPerRay = ones(1,num);
stf.totalNumOfBixels = num;

pln.propStf.couchAngles = pln.propStf.couchAngles(beamNum);
pln.propStf.gantryAngles =  pln.propStf.gantryAngles(beamNum);
w = ones(num,1);


%% dose calculation    
% analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(w,dij);
    anaDose     = resultGUI.physicalDose;

%  % analytical dose with fine sampling
    pln.propDoseCalc.anaMode = 'fineSampling';    
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 21;
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(w,dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI_FS.physicalDose;
    
    % analytical dose with std correction
    pln.propDoseCalc.anaMode = 'stdCorr';    
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(w,dijSC);
    resultGUI.physicalDoseSC = resultGUI_SC.physicalDose;
    anaScDose   = resultGUI_SC.physicalDose;
    
    % analytical dose with std correction
    pln.propDoseCalc.anaMode = 'both';    
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 21;
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    dijBO = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_BO = matRad_calcCubes(w,dijBO);
    resultGUI.physicalDoseBO = resultGUI_BO.physicalDose;
    anaBoDose   = resultGUI_BO.physicalDose;

 % Monte Carlo dose
    pln.propDoseCalc.anaMode = 'standard';
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,w,1000000);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;

       
%% execute gamma tests
ana = sum(abs(anaDose - mcDose), 'all') / max(mcDose(:));
anaFs = sum(abs(anaFsDose - mcDose), 'all') / max(mcDose(:));
anaSc = sum(abs(anaScDose - mcDose), 'all') / max(mcDose(:));
anaBo = sum(abs(anaBoDose - mcDose), 'all') / max(mcDose(:));

gammaTest = [2,2];
interpolation = 1;

figure
[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z) ,interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'stadard dose'});


figure
[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaBoDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z) ,interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corr fine sampling dose'});


figure
[gammaCube2,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corrected dose'});

figure
[gammaCube4,gammaPassRateCell] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'fine sampling dose'});
