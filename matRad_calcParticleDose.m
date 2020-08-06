function [dij, counter] = matRad_calcParticleDose(ct,stf,pln,cst)
% matRad particle dose calculation wrapper
% 
% call
%   dij = matRad_calcParticleDose(ct,stf,pln,cst)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   pln:            matRad plan meta information struct
%   cst:            matRad cst struct
%   calcDoseDirect: boolian switch to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] http://iopscience.iop.org/0031-9155/41/8/005
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

global matRad_cfg;
matRad_cfg =  MatRad_Config.instance();

% init dose calc
matRad_calcDoseInit;

% helper function for energy selection
round2 = @(a,b)round(a*10^b)/10^b;

if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
   
        alphaDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
        betaDoseTmpContainer  = cell(numOfBixelsContainer,dij.numOfScenarios);
        for i = 1:dij.numOfScenarios
            dij.mAlphaDose{i}    = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
            dij.mSqrtBetaDose{i} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
        end
        
elseif isequal(pln.propOpt.bioOptimization,'const_RBExD') && strcmp(pln.radiationMode,'protons')
            dij.RBE = 1.1;
            matRad_cfg.dispInfo('matRad: Using a constant RBE of %g\n',dij.RBE);   
end

if isfield(pln,'propDoseCalc') && ...
   isfield(pln.propDoseCalc,'calcLET') && ...
   pln.propDoseCalc.calcLET
  if isfield(machine.data,'LET')
    letDoseTmpContainer = cell(numOfBixelsContainer,dij.numOfScenarios);
    % Allocate space for dij.dosexLET sparse matrix
    for i = 1:dij.numOfScenarios
        dij.mLETDose{i} = spalloc(dij.doseGrid.numOfVoxels,numOfColumnsDij,1);
    end
  else
    warndlg('LET not available in the machine data. LET will not be calculated.');
  end
end

if isfield(pln,'propDoseCalc') && ~isfield(pln.propDoseCalc, 'airOffsetCorrection') 
    pln.propDoseCalc.airOffsetCorrection = true;
end

% generates tissue class matrix for biological optimization
if (isequal(pln.propOpt.bioOptimization,'LEMIV_effect') || isequal(pln.propOpt.bioOptimization,'LEMIV_RBExD')) ... 
        && strcmp(pln.radiationMode,'carbon')
    
    if   isfield(machine.data,'alphaX') && isfield(machine.data,'betaX')
            
        matRad_cfg.dispInfo('matRad: loading biological base data... ');
        vTissueIndex = zeros(size(VdoseGrid,1),1);
        dij.ax       = zeros(dij.doseGrid.numOfVoxels,1);
        dij.bx       = zeros(dij.doseGrid.numOfVoxels,1);

        cst = matRad_setOverlapPriorities(cst);
    
        % resizing cst to dose cube resolution 
        cst = matRad_resizeCstToGrid(cst,dij.ctGrid.x,dij.ctGrid.y,dij.ctGrid.z,...
                                         dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z);
        % retrieve photon LQM parameter for the current dose grid voxels
        [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,dij.doseGrid.numOfVoxels,1,VdoseGrid);

        for i = 1:size(cst,1)

            % check if cst is compatiable 
            if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 

                % check if base data contains alphaX and betaX
                IdxTissue = find(ismember(machine.data(1).alphaX,cst{i,5}.alphaX) & ...
                                 ismember(machine.data(1).betaX,cst{i,5}.betaX));

                % check consitency of biological baseData and cst settings
                if ~isempty(IdxTissue)
                    isInVdoseGrid = ismember(VdoseGrid,cst{i,4}{1});
                    vTissueIndex(isInVdoseGrid) = IdxTissue;
                else
                    matRad_cfg.dispError('biological base data and cst inconsistent\n');
                end
                    
            else
                    vTissueIndex(row) = 1;
                    matRad_cfg.dispInfo(['matRad: tissue type of ' cst{i,2} ' was set to 1\n']);          
            end
        end
        matRad_cfg.dispInfo('done.\n');

    else
        
        matRad_cfg.dispError('base data is incomplement - alphaX and/or betaX is missing');
        
    end
    
% issue warning if biological optimization not possible
elseif sum(strcmp(pln.propOpt.bioOptimization,{'LEMIV_effect','LEMIV_RBExD'}))>0 && ~strcmp(pln.radiationMode,'carbon') ||...
       ~strcmp(pln.radiationMode,'protons') && strcmp(pln.propOpt.bioOptimization,'const_RBExD')
    warndlg([pln.propOpt.bioOptimization ' optimization not possible with ' pln.radiationMode '- physical optimization is carried out instead.']);
    pln.propOpt.bioOptimization = 'none';      
end

% lateral cutoff for raytracing and geo calculations
effectiveLateralCutoff = matRad_cfg.propDoseCalc.defaultGeometricCutOff;

matRad_cfg.dispInfo('matRad: Particle dose calculation...\n');
counter = 0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beam = struct;

% book keeping
dij.bixelNum = [];
dij.rayNum = [];
dij.beamNum = [];
sigmaSub = 1;

% set physical dose to zero 
dij.physicalDose{1} = [];

for i = 1:length(stf) % loop over all beams
    
    % init beam
    matRad_calcDoseInitBeam;
    
    matRad_cfg.dispInfo('matRad: calculate fine sampling weights... ');
    f = waitbar(0,['Calculating weights, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    
    % create fine sampling grid
    gridsize = [1.25, 1.25];
    [gridX, gridY] = matRad_createFineSamplingGrid(stf(i), gridsize);
    
    % initialise weight container
    weightContainer.weights = [];
    weightContainer.eneNum = [];
    weightContainer.bixelNum = [];
    weightContainer.rayNum = [];
    weightContainer.beamNum = [];
    
    for j = 1:stf(i).numOfRays % loop over all rays
        
        waitbar(j/stf(i).numOfRays,f,['Calculating weights, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
        
        % calculate initial sigma and energy for all bixels on ray j
        sigmaInitial = matRad_calcSigmaIni(machine.data,stf(1).ray(j),stf(i).ray(j).SSD);
        rayEnergies = stf(i).ray(j).energy;
        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
            
            % calculate weights for every bixel for every grid ray
            % previously created
            weights = matRad_calcGriddedWeights([stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)], ...
                                            gridsize, -gridX, -gridY, sigmaInitial(k), sigmaSub);
            [~, e] = intersect([machine.data(:).energy], rayEnergies(k));
            weightContainer.weights    = [weightContainer.weights, weights];
            
            % book keeping for bixel identification
            weightContainer.eneNum   = [weightContainer.eneNum   ,e];
            weightContainer.bixelNum = [weightContainer.bixelNum ,k];
            weightContainer.rayNum   = [weightContainer.rayNum   ,j];
            weightContainer.beamNum  = [weightContainer.beamNum  ,i];   
        end
    end
    close(f);
    matRad_cfg.dispInfo('done.\n');
    
    % throw out all weights that summed up make less than $weightCutOff$ of
    % all weights 
    sortingWeights = [];
    weightCutOff = 0.003; % cutoff level for weight selection
    for a = 1:size(weightContainer.weights, 2)
        
        [tmpW, tmpIx] = sort(weightContainer.weights(:,a));
        wSum  = cumsum(tmpW);
        useIx = cumsum(tmpW) > wSum(end) * weightCutOff;
        sortingWeights = [sortingWeights; tmpIx(useIx), ones(size(tmpIx(useIx))) * weightContainer.eneNum(a)];
    end
    calculateIx = unique(sortingWeights,'rows');

    % rearange stucture of selected weights to allow for better indexing in
    % for loops 
    [ab,ac] = unique(calculateIx(:,1));
    neededIx = cell(size(ab,1),2);
    for a = 1:size(ab,1)
        neededIx{a,1} = ab(a);
        if a == size(ab,1)
            neededIx{a,2} = calculateIx(ac(a):end,2);
        else
            neededIx{a,2} = calculateIx(ac(a):ac(a+1)-1,2);
        end
    end
    
    % Determine lateral cutoff
    matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');
    densityThreshold = matRad_cfg.propDoseCalc.ssdDensityThreshold;
    cutOffLevel = matRad_cfg.propDoseCalc.defaultLateralCutOff;
    cuttOffLevel = 1;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    energies = unique(weightContainer.eneNum);
    maxE = max(weightContainer.eneNum(:));
    maxLateralCutoffDoseCalc = max(machine.data(maxE).LatCutOff.CutOff);
    matRad_cfg.dispInfo('done.\n');

    % calculate SSD needed for air offset correction, for every grid ray
    SSDs = matRad_computeGridSSD(pln.propStf.gantryAngles(i), pln.propStf.couchAngles(i), stf(i),ct, gridX, gridY);

    % loop over all grid rays that need to be calculated according to
    % weight selection before
    rayCounter = 1;    
    f = waitbar(0,['Calculating unweighted grid dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    gridDoseContainer = cell(size(gridX,1), size(energies,1), 2);
    for ixGridRay = [neededIx{:,1}]
        waitbar(ixGridRay/size(weightContainer.weights,1),f,['Calculating unweighted grid dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);

        % calculate radial distances according to lateral cutoff
        [ix,currRadialDist_sq,~,~,~,~] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                                         stf(i).sourcePoint_bev, ...
                                                         [2*gridX(ixGridRay), -stf(i).sourcePoint_bev(2), 2*gridY(ixGridRay)], ...
                                                         machine.meta.SAD, ...
                                                         find(~isnan(radDepthVdoseGrid{1})), ...
                                                         maxLateralCutoffDoseCalc);
        
        % calculate projected coordinates, representing a ray cast on the
        % now selected new grid ray, using the first ray cast done in
        % matRad_calcDoseInitBeam
        projCoords = matRad_projectOnComponents(VdoseGrid(ix), size(radDepthsMat{1}), stf(i).sourcePoint_bev,...
                                    -stf(i).sourcePoint_bev, stf(i).isoCenter,...
                                    [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                    gridX(ixGridRay), gridY(ixGridRay), rotMat_system_T);
         
        % throw out Coordinates that would result in NaN values in radDepths                        
        projCoords(projCoords(:,1) <= dij.doseGrid.resolution.x,1) = dij.doseGrid.resolution.x + 0.0001;
        projCoords(projCoords(:,2) <= dij.doseGrid.resolution.y,2) = dij.doseGrid.resolution.y + 0.0001;
        projCoords(projCoords(:,3) <= dij.doseGrid.resolution.z,3) = dij.doseGrid.resolution.z + 0.0001;
                       
        % calculate air offset, that accounts for matRad transporting
        % through vacuum from nozzle to skin and not through air
        if ~isfield(machine.meta, 'fitAirOffset') 
            fitAirOffset = 0;
        else
            fitAirOffset = machine.meta.fitAirOffset;
        end
        if ~isfield(machine.meta, 'BAMStoIsoDist') 
            BAMStoIsoDist = 400;
        else
            BAMStoIsoDist = machine.meta.BAMStoIsoDist;
        end

        nozzleToSkin = ((SSDs(ixGridRay).SSD + BAMStoIsoDist) - machine.meta.SAD);
        dR = 0.0011 * (nozzleToSkin - fitAirOffset);
             
        % interpolate actual radDepths for grid ray using projected
        % coordinates
        radDepths = interp3(radDepthsMat{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                        projCoords(:,2,:)./dij.doseGrid.resolution.y,...
                        projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest') + dR;
         
        % sanity check after applying air offset correction
        radDepths(radDepths < 0) = 0;            
        
        % loop over all needed energies on this particular grid ray
        for ixEne = [neededIx{rayCounter,2}]'
            
           % find depth depended lateral cut off
            if cutOffLevel >= 1
                currIx = radDepths <= machine.data(ixEne).depths(end);
            elseif cutOffLevel < 1 && cutOffLevel > 0
                % perform rough 2D clipping
                currIx = radDepths <= machine.data(ixEne).depths(end) & ...
                     currRadialDist_sq <= max(machine.data(ixEne).LatCutOff.CutOff.^2);

                % peform fine 2D clipping  
                if length(machine.data(ixEne).LatCutOff.CutOff) > 1
                    currIx(currIx) = matRad_interp1((machine.data(ixEne).LatCutOff.depths)',...
                        (machine.data(ixEne).LatCutOff.CutOff.^2)', radDepths(currIx)) >= currRadialDist_sq(currIx);
                end
            else
                matRad_cfg.dispError('cutoff must be a value between 0 and 1')
            end
            
            currRadDepths = radDepths(currIx);
            
            bixelDose = matRad_calcParticleDoseBixel(...
                    currRadDepths, ...
                    currRadialDist_sq(currIx), ...
                    sigmaSub^2, ...
                    machine.data(ixEne));
            counter = counter + 1;
            
            % fill grid dose container
            [~, xe] = intersect(energies, ixEne);
            gridDoseContainer{ixGridRay,xe,1} = bixelDose;
            gridDoseContainer{ixGridRay,xe,2} = VdoseGrid(ix(currIx));
        end
        
    rayCounter = rayCounter + 1;
    end
    close(f);    
    
    bixelContainer = cell(size(weightContainer.weights,1),1);
    
    % calculate dose for every bixel, using grid dose container and making
    % everything sparse
    f = waitbar(0,['Calculating weighted dij, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    for t = 1:size(weightContainer.weights,2)
       waitbar(t/size(weightContainer.weights,2),f,['Calculating weighted dij, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);

       weights = weightContainer.weights(:,t); 
       [~, xe] = intersect(energies, weightContainer.eneNum(t));
       rayDose = gridDoseContainer(:,xe,1);
       rayIx   = gridDoseContainer(:,xe,2);
       for g = 1:size(rayDose,1)
           if isempty(rayIx{g,1}) 
                bixelContainer{g} = zeros(0,3);   
           else
                bixelContainer{g} = [rayIx{g,1}, ones(size(rayIx{g,1},1),1), rayDose{g,1} * weights(g)];
           end
       end
       
       IJV = cell2mat(bixelContainer);
       dij.physicalDose{1} = [dij.physicalDose{1}, sparse(IJV(:,1),IJV(:,2),IJV(:,3),prod(dij.doseGrid.dimensions),1)];
       dij.bixelNum = [dij.bixelNum; weightContainer.bixelNum(t)];
       dij.rayNum   = [dij.rayNum;   weightContainer.rayNum(t)];
       dij.beamNum  = [dij.beamNum;   weightContainer.beamNum(t)];
    end
    close(f);
end
