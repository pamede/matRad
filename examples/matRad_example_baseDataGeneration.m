%% Example: Fit analytical baseData to mcSquare simulation 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Create CT and CST via MatRad_BaseDataGenerator
clear 
cubeDim =       [100, 100, 100];    % dimesions of CT cube (x,y,z)
cubeRes =       [  3,   3,   3];  	% resolution CT cube in each dimension in [mm]    
targetPos =     [150, 150, 150];    % center of target in [mm]
targetSize =    [ 10,  10,  30];    % target size in [mm]

% call BaseDataGenerator constructor, creating the ct and cst struct
BDG = MatRad_BaseDataGenerator(cubeDim, cubeRes, targetPos, targetSize);

% extract exemplary emittance and energy spectrum parameterization from
% matRad base data file
load protons_generic_emmitance.mat
energySpectrum = [machine.data(:).energySpectrum];
tmp  = [machine.data(:).initFocus];
emittanceData  = [tmp(:).Emittance];

% define machine specific parameters, bear in mind. nozzleToIso is
% equivalent to nozzle to skin distance, as isocenter in simulation is set
% at phantom entrance
SAD         = 2028.5;
nozzleToIso = 420.0;
focalPoints = [2218, 1839];

% define energy steps for fitted base data set
fitEnergies = [75, 125, 175];

% number of particles in Monte Carlo simulations
N = 1e5;

% call fitting funtion using previous emittance parameterization
fitMachine = BDG.fitAnalyticalDataViaEmittance(nozzleToIso, SAD, focalPoints, ...
                            energySpectrum, emittanceData, fitEnergies, N)
                        
% call fitting funtion using existing MCsquare base data file
fitMachine1 = BDG.fitAnalyticalDataViaFile(nozzleToIso, SAD, ...
                            'BDL_matRad.txt', fitEnergies, N);

