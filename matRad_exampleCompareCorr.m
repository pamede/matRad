
matRad_rc

% load patient data, i.e. ct, voi, cst
% load lung01MC.mat 
load alderson01MC.mat
% load lung02.mat
num = 159;
stf = stf(1);
stf.ray = stf.ray(num);
stf.numOfRays = 1;
stf.numOfBixelsPerRay = stf.numOfBixelsPerRay(num);
stf.totalNumOfBixels = sum([stf.numOfBixelsPerRay(:)]);
pln.propStf.numOfBeams = pln.propStf.numOfBeams(1); 
pln.propStf.couchAngles = pln.propStf.couchAngles(1); 
pln.propStf.gantryAngles = pln.propStf.gantryAngles(1);

% load slab01MC.mat
 
%% dose calculation    
% analytical dose without fine sampling
    pln.propDoseCalc.anaMode = 'standard';
    dij = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI = matRad_calcCubes(1,dij);
    anaDose     = resultGUI.physicalDose;

%  % analytical dose with fine sampling
    pln.propDoseCalc.anaMode = 'fineSampling';    
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 11;
    pln.propDoseCalc.fineSampling.sigmaSub = 2;
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(1,dijFS);
    anaFsDose   = resultGUI_FS.physicalDose;
    
    % analytical dose with std correction
    pln.propDoseCalc.anaMode = 'stdCorr';    
    dijSC = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_SC = matRad_calcCubes(1,dijSC);
    anaScDose   = resultGUI_SC.physicalDose;

 % Monte Carlo dose
    pln.propDoseCalc.anaMode = 'standard';
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,1, 1000000);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    
%% execute gamma tests
gammaTest = [2, 2];
interpolation = 2;

figure
[gammaCube1,gammaPassRateCell] = matRad_gammaIndex(anaDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z) ,interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'stadard dose'});

figure
[gammaCube2,gammaPassRateCell] = matRad_gammaIndex(anaScDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'std corrected dose'});

figure
[gammaCube4,gammaPassRateCell] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,round(stf(1).isoCenter(3)/ct.resolution.z),interpolation,'global',cst);
title({[num2str(gammaPassRateCell{1,2}) '% of points > ' num2str(gammaTest(1)) '% pass gamma criterion (' num2str(gammaTest(1)) '%/ ' num2str(gammaTest(2)) 'mm)'];'fine sampling dose'});
