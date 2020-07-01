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

CalcCounter = 0;
dij.bixelNum = [];
dij.rayNum = [];
dij.beamNum = [];
sigmaSub = 1;
matRad_cfg.dispInfo('matRad: calculate fine sampling weights... ');
physicalDose = [];
beamWiseBixelCounter = 0;
% tic
for i = 1:length(stf) % loop over all beams
    
    % init beam
    matRad_calcDoseInitBeam;
    
    f = waitbar(0,['Calculating weights, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    
    weightedGrid.energy = [];
    gridsize = [2, 2];
    [gridX, gridY] = matRad_createFineSamplingGrid(stf(i), gridsize);
    weightContainer.weights = [];
    weightContainer.eneNum = [];
    weightContainer.bixelNum = [];
    weightContainer.rayNum = [];
    weightContainer.beamNum = [];
    
    counter = 1;
    for j = 1:stf(i).numOfRays % loop over all rays
        waitbar(j/stf(i).numOfRays,f,['Calculating weights, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
        sigmaInitial = matRad_calcSigmaIni(machine.data,stf(1).ray(j),stf(i).ray(j).SSD);
        rayEnergies = stf(i).ray(j).energy;
        for k = 1:stf(i).numOfBixelsPerRay(j) % loop over all bixels per ray
            weights = matRad_calcGriddedWeights([stf(i).ray(j).rayPos_bev(1) stf(i).ray(j).rayPos_bev(3)], ...
                                            gridsize, -gridX, -gridY, sigmaInitial(k), sigmaSub);
            [~, e] = intersect([machine.data(:).energy], rayEnergies(k));
            weightContainer.weights    = [weightContainer.weights, weights];
            weightContainer.eneNum   = [weightContainer.eneNum   ,e];
            weightContainer.bixelNum = [weightContainer.bixelNum ,k];
            weightContainer.rayNum   = [weightContainer.rayNum   ,j];
            weightContainer.beamNum  = [weightContainer.beamNum  ,i];   
        end
    end
    close(f);
    matRad_cfg.dispInfo('done.\n');
    
    % Determine lateral cutoff
    matRad_cfg.dispInfo('matRad: calculate lateral cutoff...');
    cutOffLevel = matRad_cfg.propDoseCalc.defaultLateralCutOff;
    cuttOffLevel = 1;
    visBoolLateralCutOff = 0;
    machine = matRad_calcLateralParticleCutOff(machine,cutOffLevel,stf(i),visBoolLateralCutOff);
    matRad_cfg.dispInfo('done.\n');
        
    % missing air offset correction
    dR = 0;
    
    energies = unique(weightContainer.eneNum);
    maxE = max(weightContainer.eneNum(:));
    maxLateralCutoffDoseCalc = max(machine.data(maxE).LatCutOff.CutOff);
    
    counter = 1;
%     tmp = cell(max(weightContainer.rayNum), max(weightContainer.bixelNum));
    
    f = waitbar(0,['Calculating unweighted grid dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    test = [];
    for ixGridRay = 1:size(weightContainer.weights,1)
        waitbar(ixGridRay/size(weightContainer.weights,1),f,['Calculating unweighted grid dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
%         ixGridRay
        [ix,currRadialDist_sq,~,~,~,~] = matRad_calcGeoDists(rot_coordsVdoseGrid, ...
                                                         stf(i).sourcePoint_bev, ...
                                                         [2*gridX(ixGridRay), -stf(i).sourcePoint_bev(2), 2*gridY(ixGridRay)], ...
                                                         machine.meta.SAD, ...
                                                         find(~isnan(radDepthVdoseGrid{1})), ...
                                                         maxLateralCutoffDoseCalc);
                                                     
        projCoords = matRad_projectOnComponents(VdoseGrid(ix), size(radDepthsMat{1}), stf(i).sourcePoint_bev,...
                                    -stf(i).sourcePoint_bev, stf(i).isoCenter,...
                                    [dij.doseGrid.resolution.x dij.doseGrid.resolution.y dij.doseGrid.resolution.z],...
                                    gridX(ixGridRay), gridY(ixGridRay), rotMat_system_T);
         
%         [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf(i).isoCenter, ...
%                                  ct.resolution, ...
%                                  stf(i).sourcePoint, ...
%                                  stf(i).ray(j).targetPoint, ...
%                                  {ct.cube{1}});
%         ixSSD = find(rho{1} > densityThreshold,1,'first');
% 
%         % calculate SSD
%         SSD{j} = double(2 * stf(i).SAD * alpha(ixSSD));
            
            
%         if ~isfield(machine.meta, 'fitAirOffset') 
%             fitAirOffset = 0;
%         else
%             fitAirOffset = machine.meta.fitAirOffset;
%         end
% 
%         if ~isfield(machine.meta, 'BAMStoIsoDist') 
%             BAMStoIsoDist = 400;
%         else
%             BAMStoIsoDist = machine.meta.BAMStoIsoDist;
%         end
% 
%         nozzleToSkin = ((stf(i).ray(j).SSD + BAMStoIsoDist) - machine.meta.SAD);
%         dR = 0.0011 * (nozzleToSkin - fitAirOffset);
   
                                
        radDepths = interp3(radDepthsMat{1},projCoords(:,1,:)./dij.doseGrid.resolution.x,...
                        projCoords(:,2,:)./dij.doseGrid.resolution.y,...
                        projCoords(:,3,:)./dij.doseGrid.resolution.z,'nearest') + dR;
                  
        for ixEne = energies
%             weights = weightContainer.weights(ixGridRay,:);
%             [~ ,index] = intersect(weightContainer.eneNum, ixEne);
%             beamNum  = weightContainer.beamNum(index);
%             bixelNum = weightContainer.bixelNum(index);
%             rayNum   = weightContainer.rayNum(index);
            
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
                
            if isempty(bixelDose)
                test = [test; gridX(ixGridRay),gridY(ixGridRay)];
                
                1+1;
            end
            
%             tmp{rayNum, bixelNum} = [tmp{rayNum, bixelNum}, sparse(VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1)];
%             tmp{counter,1} = sparse(VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);
            
            [~, xe] = intersect(energies, ixEne);
            tmp{ixGridRay,xe,1} = bixelDose;
            tmp{ixGridRay,xe,2} = VdoseGrid(ix(currIx));
            
            % Save dose for every bixel in cell array
%             doseTmpContainer{counter,1} = sparse(VdoseGrid(ix(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);
%             tmp = [tmp; counter, ixEne];
        counter = counter + 1;
        end
        

    end
    close(f);
%     scatter(test(:,1), test(:,2))
    
    
    dij.physicalDose{1} = [];
    bixelContainer = cell(size(weightContainer.weights,1),1);
    
    f = waitbar(0,['Calculating weighted dij, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
    for t = 1:size(weightContainer.weights,2)
       waitbar(t/size(weightContainer.weights,2),f,['Calculating weighted dij, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);

       weights = weightContainer.weights(:,t); 
       [~, xe] = intersect(energies, weightContainer.eneNum(t));
       rayDose = tmp(:,xe,1);
       rayIx   = tmp(:,xe,2);
       for g = 1:size(rayDose,1)
           bixelContainer{g} = [rayIx{g,1}, ones(size(rayIx{g,1},1),1), rayDose{g,1} * weights(g)];
       end
       IJV = cell2mat(bixelContainer);
       dij.physicalDose{1} = [dij.physicalDose{1}, sparse(IJV(:,1),IJV(:,2),IJV(:,3),prod(dij.doseGrid.dimensions),1)];
       dij.bixelNum = [dij.bixelNum; weightContainer.bixelNum(t)];
       dij.rayNum   = [dij.rayNum;   weightContainer.rayNum(t)];
       dij.beamNum  = [dij.beamNum;   weightContainer.beamNum(t)];
%        tete = [bixelContainer, ones(size(weightContainer.weights,1),1), rayIx];
%        
%        physicalDose = [physicalDose, sparse(
           
%        bixelContainer{t,:} = rayDose;
%        cell2mat( tmp(:,xe,:));
    end
    close(f);
end
                                                 
                                                 
   
    
    
    
    
    
%     
%     
%     %     dijFS =  sparse(prod(dij.doseGrid.dimensions), size(gridX,1));
%     dijPosition = beamWiseBixelCounter;
%     f = waitbar(0,['Calculating fine sampling dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
%     for ixEne = 1:size(weightedGrid, 2)
%         
%         dijFS  =  sparse(prod(dij.doseGrid.dimensions), size(gridX,1));
%         dijFS1 =  {[]};
%         
%         energy = weightedGrid(ixEne).energy;
%     
%         [~, energyIx] = intersect([machine.data(:).energy], energy);  
%         
%         gridIx = [];
% %         tic
%         for ixBixel = 1:size(weightedGrid(ixEne).weights,2)
%             cutOff = 0.05;
%     %         bixelFactor = 1+cutOff/(1-cutOff);
%             [weights, ixRe] = sort(weightedGrid(ixEne).weights(:,ixBixel));
%             cumWeights = cumsum(weights);
%             cutIx = cumWeights > cutOff * cumWeights(end);
%             gridIx = [gridIx; ixRe(cutIx)];
%         end
%         gridIx = unique(gridIx);
% %         toc
%        
% %         sum(weights)
% %         sum(weights(cutIx)) * bixelFactor
%         
% %         scatter(gridX, gridY);
% %         hold on
% %         scatter(gridX(gridIx), gridY(gridIx));
% %         hold off
% %         
% %         weightSum = sum(weightedGrid(ixEne).weights,'all')/size(gridX,1);
% %         neededIx = gridWeights > weightSum / 10 ;
% %         gridArray = find(neededIx);
%         
%         for ixGrid = gridIx'
%             
% %             tmp = zeros(prod(dij.doseGrid.dimensions),1);
% %             tmp = full(dijFS(:,ix));
% % %             tmp(VdoseGrid(availableIx)) = radDepths(:,1,1);
% %             tmo1 = reshape(tmp,dij.doseGrid.dimensions);
% %             imagesc(tmo1(:,:,12));
% %             tic
%             radialDistSq = matRad_calcRadialDistanceSq(rot_coordsVdoseGrid, [0 -stf(i).SAD 0], [0 stf(i).SAD 0],gridX(ixGrid),gridY(ixGrid));           
% %             toc
% %             rad_distancesSq = radialDistSq(availableIx,ixGrid);
%             rad_distancesSq = radialDistSq(availableIx);
% 
%             currIx = radDepths(:,:,ixGrid) < machine.data(energyIx).depths(end);  
%             
%             % calculate particle dose for bixel k on ray j of beam i
% %             tic
%             bixelDose = matRad_calcParticleDoseBixel(...
%                 radDepths(currIx,1,ixGrid), ...
%                 rad_distancesSq(currIx), ...
%                 sigmaSub^2, ...
%                 machine.data(energyIx));
% %             bixeltim = toc
%             CalcCounter = CalcCounter + 1;
% %             tic
% %             tic
% %             dijFS(:,ixGrid+1)                                           = sparse(VdoseGrid(availableIx(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1); % 10 times slower than bixel dose
% %             t1 = toc
% %             tic
%             dijFS1{ixGrid,1} = sparse(VdoseGrid(availableIx(currIx)),1,bixelDose,dij.doseGrid.numOfVoxels,1);
% %             t2 = toc
%             %             doseTmpContainer{mod(counter-1,numOfBixelsContainer)+1,1} =  sparse(VdoseGrid(ix),1,totalDose,dij.doseGrid.numOfVoxels,1);
% %             other fs code dij filling
% %             dijtim = toc
% %             doseContainer(VdoseGrid(availableIx(currIx)),:) = doseContainer(VdoseGrid(availableIx(currIx)),:) ...
% %                         + sparse(weightedGrid(ixEne).weights(ixGrid,:) .* bixelDose); 
%         end 
% %         tic
%         waitbar(ixEne/size(weightedGrid, 2),f,['Calculating fine sampling dose, beam ' num2str(i) ' of ' num2str(length(stf)) '...']);
% %         tic
%         dijPosition = dijPosition(end)+1:dijPosition(end)+size(weightedGrid(ixEne).weights,2);
% %         dij.physicalDose{1}(:,dijPosition) = dijFS * weightedGrid(ixEne).weights;
% % tic
%         for ixCell = 1:size(dijFS1,1)
%             if ~isempty(dijFS1{ixCell,1})
%                 dij.physicalDose{1}(:,dijPosition) = dij.physicalDose{1}(:,dijPosition) + dijFS1{ixCell,1} * weightedGrid(ixEne).weights(ixCell,:);
%             end
%         end
% % t3 =        toc
%         dij.bixelNum = [dij.bixelNum; weightedGrid(ixEne).bixelNum'];
%         dij.rayNum   = [dij.rayNum;   weightedGrid(ixEne).rayNum'];
%         dij.beamNum  = [dij.beamNum;  weightedGrid(ixEne).beamNum'];
% %         filltime = toc
%     end 
%     close(f);
%     matRad_cfg.dispInfo('done.\n'); 
%     beamWiseBixelCounter = beamWiseBixelCounter + stf(i).totalNumOfBixels;
%     CalcCounter
% end
% 
% 
% 
% 
