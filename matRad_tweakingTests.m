clear
matRad_rc

% load patient data, i.e. ct, voi, cst
load WATERPHANTOM_res3.mat

% meta information for treatment plan
pln.radiationMode   = 'protons';     
pln.machine         = 'Generic';

pln.numOfFractions  = 30;

% beam geometry settings
pln.propStf.bixelWidth      = 500; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.longitudinalSpotSpacing = 500;
pln.propStf.gantryAngles    = [0]; % [?] 
pln.propStf.couchAngles     = [0]; % [?]
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propStf.isoCenter(2)    = 0;

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

ixEnergy = 45;

%% generate steering file
stf = matRad_generateStf(ct,cst,pln);
load protons_Generic.mat
stf.ray.energy = machine.data(ixEnergy).energy;

%% dose calculation
 % initial analytical dose calculation
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
    resultGUI = matRad_calcCubes(ones(sum([stf(:).totalNumOfBixels]),1),dij);
    anaDose     = resultGUI.physicalDose;
    anaIDD = sum(anaDose,[2 3]);

 % initial Monte Carlo dose calculation
    resultGUI_MC = matRad_calcDoseDirectMC(ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), 1e4);
    resultGUI.physicalDoseMC = resultGUI_MC.physicalDose;
    mcDose      = resultGUI.physicalDoseMC;
    mcIDD = sum(mcDose,[2 3]);
    
 % energy dependent Monte Carlo calculation 
    load protons_temporaryMachine.mat
    mean        = machine.data(ixEnergy).monteCarloData.MeanEnergy;
    spread      = machine.data(ixEnergy).monteCarloData.EnergySpread;
    spotsize    = machine.data(ixEnergy).monteCarloData.SpotSize1x;
    divergence   = machine.data(ixEnergy).monteCarloData.Divergence1x;
    correlation = machine.data(ixEnergy).monteCarloData.Correlation1x;
    N      = 1e5;

    tmpFig = figure;
	set(gcf,'Position',[100 400 1500 500])
    
%     testFcn = @(a,b,c) display('Hallo');
%     'OutputFcn', testFcn
    
    objectiveFunctionEnergy = @(mean, spread) matRad_calcMCsquareObjectiveEnergy(ct, stf, pln, ...
                        cst, N, anaDose, mean, spread, tmpFig);
    optionsEnergy = optimset('Display','iter', 'TolFun', 0.01, 'TolX', 0.005, 'PlotFcns',@optimplotfval);
    x = fminsearch(@(x) objectiveFunctionEnergy(x(1), x(2)), [mean, spread], optionsEnergy); 
    
    foundMean = x(1);
    foundSpread = x(2);
    
%  optics dependent Monte Carlo calculation

    N      = 5e5;
    objectiveFunctionOptic = @(spot, div, corr) matRad_calcMCsquareObjectiveOptics(ct, stf, pln, cst, N, anaDose, ...
                                            foundMean, foundSpread, spot, div, corr, tmpFig);
    optionsOptic = optimset('Display','iter', 'TolFun', 0.003, 'TolX', 0.01, 'PlotFcns',@optimplotfval);
    x = fminsearch(@(x) objectiveFunctionOptic(x(1), x(2), x(3)), [spotsize, divergence, correlation], optionsOptic);
    
    foundSpotsize = 	x(1);
    foundDivergence = 	x(2);
    foundCorrelation = 	x(3);

    
    

 % all parameter dependent Monte Carlo calculation
    objectiveFunctionAll = @(mean, spread, spot, div, corr) matRad_calcMCsquareObjectiveAll(ct, stf, pln, cst, ...
                                        5e6, anaDose, mean, spread, spot, div, corr, tmpFig);
    optionsAll = optimset('Display','iter', 'TolFun', 0.003, 'TolX', 0.01, 'PlotFcns',@optimplotfval);
    x = fminsearch(@(x) objectiveFunctionAll(x(1), x(2), x(3), x(4), x(5)), [foundMean, foundSpread, foundSpotsize, foundDivergence, foundCorrelation], optionsAll);
 
    finalMean   = x(1);
    finalSpread = x(2);
    finalSpot   = x(3);
    finalDiv    = x(4);
    finalCorr   = x(5);
    
    
 %%
 
    resultFig = figure;
    [F, mcIDD, mcFWHM] = matRad_calcMCsquareObjectiveAll(ct, stf, pln, cst, 1e6, anaDose, ...
                    finalMean, finalSpread, finalSpot, finalDiv, finalCorr, resultFig);

    
    