function dij = matRad_calcParticleDose(ct,stf,pln,cst)
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

matRad_cfg.dispInfo('matRad: calculate fine sampling weights... ');
physicalDose = [];
for i = 1:length(stf) % loop over all beams
    weightedGrid.energy = [];
    gridsize = [1,1];
    [gridX, gridY] = matRad_createFineSamplingGrid(stf(i), gridsize);

    counter = 1;
    for j = 1:stf(i).numOfRays % loop over all rays
        j/stf(i).numOfRays * 100
        sigmaInitial = matRad_calcSigmaIni(machine.data,stf(1).ray(j),stf(i).ray(j).SSD);
        rayEnergies = stf(i).ray(j).energy;
        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
            energy = rayEnergies(k);
            if isempty(intersect([weightedGrid(:).energy], energy))
                ixEnergy = size([weightedGrid.energy], 2) + 1;
                weightedGrid(ixEnergy).energy = energy;
%                 weightedGrid(ixEnergy).weights =  zeros(size(gridX));
                weightedGrid(ixEnergy).weights =  zeros(size(gridX,1), stf(i).totalNumOfBixels);
            else
                [~, ixEnergy]= intersect([weightedGrid(:).energy], energy);
            end
            weights = matRad_calcGriddedWeights([stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)], ...
                                            gridsize, -gridX, -gridY, sigmaInitial(k), sigmaSub);
%             weights(weights < max(weights)/100) =0;

%             weightedGrid(ixEnergy).weights = [weightedGrid(ixEnergy).weights(:)] + weights';  
            weightedGrid(ixEnergy).weights(:,counter) = weights';

            dij.bixelNum = [dij.bixelNum; k];
            dij.rayNum   = [dij.rayNum;   j];
            dij.beamNum  = [dij.beamNum;  i];

            counter = counter + 1;
        end
    end
    matRad_cfg.dispInfo('done.\n');

    % convert voxel indices to real coordinates using iso center of beam i
    xCoordsV       = xCoordsV_vox(:)*ct.resolution.x-stf(i).isoCenter(1);
    yCoordsV       = yCoordsV_vox(:)*ct.resolution.y-stf(i).isoCenter(2);
    zCoordsV       = zCoordsV_vox(:)*ct.resolution.z-stf(i).isoCenter(3);
    coordsV        = [xCoordsV yCoordsV zCoordsV];

    xCoordsVdoseGrid = xCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.x-stf(i).isoCenter(1);
    yCoordsVdoseGrid = yCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.y-stf(i).isoCenter(2);
    zCoordsVdoseGrid = zCoordsV_voxDoseGrid(:)*dij.doseGrid.resolution.z-stf(i).isoCenter(3);
    coordsVdoseGrid  = [xCoordsVdoseGrid yCoordsVdoseGrid zCoordsVdoseGrid];

    % Get Rotation Matrix
    % Do not transpose matrix since we usage of row vectors &
    % transformation of the coordinate system need double transpose

    rotMat_system_T = matRad_getRotationMatrix(stf(i).gantryAngle,stf(i).couchAngle);

    % Rotate coordinates (1st couch around Y axis, 2nd gantry movement)
    rot_coordsV         = coordsV*rotMat_system_T;
    rot_coordsVdoseGrid = coordsVdoseGrid*rotMat_system_T;

    rot_coordsV(:,1) = rot_coordsV(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsV(:,2) = rot_coordsV(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsV(:,3) = rot_coordsV(:,3)-stf(i).sourcePoint_bev(3);

    rot_coordsVdoseGrid(:,1) = rot_coordsVdoseGrid(:,1)-stf(i).sourcePoint_bev(1);
    rot_coordsVdoseGrid(:,2) = rot_coordsVdoseGrid(:,2)-stf(i).sourcePoint_bev(2);
    rot_coordsVdoseGrid(:,3) = rot_coordsVdoseGrid(:,3)-stf(i).sourcePoint_bev(3);

    % calculate geometric distances
    geoDistVdoseGrid{1}= sqrt(sum(rot_coordsVdoseGrid.^2,2));

    % Calculate radiological depth cube
    matRad_cfg.dispInfo('matRad: calculate radiological depth cube... ');
    [radDepthVctGrid, radDepthsMat] = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,effectiveLateralCutoff);
    
    % interpolate radiological depth cube to dose grid resolution
    radDepthVdoseGrid = matRad_interpRadDepth...
    (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);

    % interpolate radiological depth cube used for fine sampling to dose grid resolution
    radDepthsMat{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, radDepthsMat{1}, ...
                                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
      
    %sub set of voxels for which radiological depth calculations are available
    availableIx = find(~isnan(radDepthVdoseGrid{1}));
          
    cutOffLevel = matRad_cfg.propDoseCalc.defaultLateralCutOff;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    
    doseContainer = sparse(size(radDepthVdoseGrid{1},1),stf(i).totalNumOfBixels);
    
    projCoords = matRad_projectOnComponents(VdoseGrid(availableIx), size(radDepthsMat{1}), stf(i).sourcePoint_bev,...
                                        [0, -stf(i).sourcePoint_bev(2), 0], stf(i).isoCenter,...
                                        [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                        gridX(:), gridY(:), rotMat_system_T);

    % interpolate radiological depths at projected
    % coordinates
    radDepths = interp3(radDepthsMat{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                        projCoords(:,2,:)./dij.doseGrid.resolution.y,...
                        projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest');
    matRad_cfg.dispInfo('done.\n');      
                    
    for ixEne = 1:size(weightedGrid, 2)
        ixEne/size(weightedGrid, 2)*100
        energy = weightedGrid(ixEne).energy;
    
        [~, energyIx] = intersect([machine.data(:).energy], energy);  
        
        
    
        for ixGrid = 1:size(weightedGrid(ixEne).weights, 1)
            
            sourcePoint_bev = [0 -stf(i).SAD 0]';
            targetPoint_bev = [2 * gridX(ixGrid) stf(i).SAD 2 * gridY(ixGrid)]';
            
            [~, rad_distancesSq] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                              sourcePoint_bev, ...
                                              targetPoint_bev, ...
                                              stf(i).SAD, ...
                                              radDepths, ...
                                              1000);
            
            currIx = radDepths(:,:,ixGrid) < machine.data(energyIx).depths(end);  
            
            % calculate particle dose for bixel k on ray j of beam i
            bixelDose = matRad_calcParticleDoseBixel(...
                radDepths(currIx,:,ixGrid), ...
                rad_distancesSq(currIx), ...
                sigmaSub^2, ...
                machine.data(energyIx));   
            doseContainer(availableIx(currIx),:) = doseContainer(availableIx(currIx),:) + sparse(weightedGrid(ixEne).weights(ixGrid,:) .* bixelDose);
        end    
    end 
    physicalDose = [physicalDose, doseContainer];
    clear doseContainer;
end
dij.physicalDose{1} = physicalDose;


