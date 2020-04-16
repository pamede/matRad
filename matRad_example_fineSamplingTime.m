matRad_rc
load slab01MC.mat
gammaTest = [2, 2];
interpolation = 0;
erg = [];
counter = 1;

for sigma = 0.5:0.5:3
    
    erg(1,counter) = sigma
    % analytical dose with fine sampling
    pln.propDoseCalc.anaMode = 'fineSampling';
    pln.propDoseCalc.fineSampling.method = 'russo';
    pln.propDoseCalc.fineSampling.N = 21;
    pln.propDoseCalc.fineSampling.sigmaSub = sigma;

    tic
    dijFS = matRad_calcParticleDose(ct,stf,pln,cst,false);
    resultGUI_FS = matRad_calcCubes(resultGUI.w,dijFS);
    resultGUI.physicalDoseFS = resultGUI_FS.physicalDose;
    anaFsDose   = resultGUI.physicalDoseFS;
    erg(2,counter) = toc;
    

    [~, gamma] = matRad_gammaIndex(anaFsDose,mcDose,[ct.resolution.x,ct.resolution.y,ct.resolution.z],gammaTest,1,interpolation,'local',cst);
    erg(3,counter) = gamma{1,2};
    counter = counter + 1;
end