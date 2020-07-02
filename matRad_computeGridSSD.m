function arraySSD = matRad_computeGridSSD(gantryAngle,couchAngle, stf,ct, gridX, gridY)
% matRad SSD calculation
% 
% call
%   stf = matRad_computeSSD(stf,ct,mode)
%
% input
%   ct:             ct cube
%   stf:            matRad steering information struct
%   mode:           optional parameter specifying how to handle multiple
%                   cubes to compute one SSD
% output
%   stf:            matRad steering information struct
%
% References
%
%   -
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
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
matRad_cfg = MatRad_Config.instance();

% set density threshold for SSD computation
densityThreshold = matRad_cfg.propDoseCalc.ssdDensityThreshold;

rotMat_vectors_T = transpose(matRad_getRotationMatrix(gantryAngle,couchAngle));

SSD = cell(1,size(gridX,1));
for j = 1:size(gridX,1)
    [alpha,~,rho,~,~] = matRad_siddonRayTracer(stf.isoCenter, ...
                         ct.resolution, ...
                         stf.sourcePoint, ...
                         [2*gridX(j), -stf.sourcePoint_bev(2), 2*gridY(j)] * rotMat_vectors_T, ...
                         {ct.cube{1}});
    ixSSD = find(rho{1} > densityThreshold,1,'first');

    % calculate SSD
    SSD{j} = double(2 * stf.SAD * alpha(ixSSD));
    arraySSD(j).SSD = SSD{j};            
end

% try to fix SSD by using SSD of closest neighbouring ray
SSDnotSet = find(cellfun('isempty',SSD));
if ~isempty(SSDnotSet)
    rayPos_bev = reshape([gridX, zeros(size(gridX)), gridY]',[3 size(gridX,1)])';
    for j = SSDnotSet
         arraySSD(j).SSD =  matRad_closestNeighbourSSD(rayPos_bev, SSD, rayPos_bev(j,:));
    end
end


% default setting only use first cube
function bestSSD = matRad_closestNeighbourSSD(rayPos, SSD, currPos)
    vDistances = sum((rayPos - repmat(currPos,size(rayPos,1),1)).^2,2);
    [~, vIdx]   = sort(vDistances);
    for ix = vIdx'
        bestSSD = SSD{ix};
        % if SSD has been found, bestSSD is not empty
        if ~any(isempty(bestSSD))
            break
        end
    end
    if any(isempty(bestSSD))
%         error('Could not fix SSD calculation.');
    end
end






end
