matRad_cfg.dispInfo('Beam %d of %d:\n',i,dij.numOfBeams);

% remember beam and bixel number
if calcDoseDirect
    dij.beamNum(i)    = i;
    dij.rayNum(i)     = i;
    dij.bixelNum(i)   = i;
end

bixelsPerBeam = 0;

if strcmp(anaMode, 'stdCorr')
    newIx = (1:prod(ct.cubeDim))';
    [yCoordsV_vox, xCoordsV_vox, zCoordsV_vox] = ind2sub(ct.cubeDim,newIx);
end

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

if strcmp(anaMode, 'stdCorr')
    
    [radDepthVctGrid, radDepthsMat] = matRad_rayTracing(stf(i),ct,newIx,rot_coordsV,300);
    
    B = imrotate3(radDepthsMat{1},stf(i).gantryAngle,[0,0,1],'cubic');
    
    rotRadDepths = B;
    rotRadDepths(rotRadDepths == 0) = NaN;

    cStdCtGrid   = zeros(size(rotRadDepths));
    meanRadDepths = zeros(size(rotRadDepths)); 


    sizeRot = size(rotRadDepths);
    save = [];

    kernelSize = min(sizeRot(2), sizeRot(3)) + (mod(min(sizeRot(2), sizeRot(3)),2) -  1);
    res = [ct.resolution.y, ct.resolution.z];
    dist = 0;
    cube = rotRadDepths;
    
    for ii = 1:sizeRot(1)
        if ii == 1
            lastC = 0;
            secondToLastC = 0;
            tmpCstd = 0;
        elseif ii > 1
            secondToLastC = 0;
            lastC = sum(cStdCtGrid(ii-1,:,:),'all')/numel(cStdCtGrid(ii-1,:,:));
        elseif ii > 2
            secondToLastC = sum(cStdCtGrid(ii-2,:,:),'all')/numel(cStdCtGrid(ii-1,:,:));
            lastC = sum(cStdCtGrid(ii-1,:,:),'all')/numel(cStdCtGrid(ii-1,:,:));
        end
   
        slice = reshape(cube(ii,:,:),sizeRot(2),sizeRot(3));
        
        notNaN = find(~isnan(slice));
        slice(isnan(slice)) = 0;
        distSlice = sum(slice(notNaN))/ numel(slice(notNaN));
                
        nPadd = floor(kernelSize/2);
        paddedSlice = zeros(size(slice,1) + 2 * nPadd, size(slice,2) + 2 * nPadd);
        paddedSlice(nPadd + 1:sizeRot(2) + nPadd, nPadd + 1:sizeRot(3) + nPadd) = slice;
        
        paddedSlice(1:nPadd, 1:nPadd) = slice(1,1);
        paddedSlice(sizeRot(2) + nPadd + 1:end, 1:nPadd) = slice(end,1);
        paddedSlice(1:nPadd, sizeRot(3) + nPadd + 1:end) = slice(1,end);
        paddedSlice(sizeRot(2) + nPadd + 1:end, sizeRot(3) + nPadd + 1:end) = slice(end);
        
        paddedSlice(1:nPadd,nPadd + 1:sizeRot(3) + nPadd) = repmat(slice(1,:),nPadd,1);
        paddedSlice(sizeRot(2) + nPadd + 1:end,nPadd + 1:sizeRot(3) + nPadd) = repmat(slice(end,:),nPadd,1);
        paddedSlice(nPadd + 1:sizeRot(2) + nPadd,1:nPadd) = repmat(slice(:,1),1,nPadd);
        paddedSlice(nPadd + 1:sizeRot(2) + nPadd ,sizeRot(3) + nPadd + 1:end) = repmat(slice(:,end),1,nPadd);

        if isnan(dist)
            sigma = 0;
        else
            sigma = 1.2 * (1 + 0.5 * (lastC - secondToLastC));
        end
        
        save = [save, sigma];
        
        kernel = matRad_create2dimGaussKernel(kernelSize, sigma, res);

        tmpCstd     = (convn(paddedSlice.^2, kernel, 'same') - convn(paddedSlice, kernel, 'same').^2);
        tmpCstd(tmpCstd < 0) = 0;
        tmpCstd = sqrt(tmpCstd);
        tmpRadDepth = convn(paddedSlice, kernel, 'same');
        cStdCtGrid(ii,:,:)   = tmpCstd((kernelSize+1)/2:size(tmpCstd,1)-(kernelSize-1)/2,(kernelSize+1)/2:size(tmpCstd,2)-(kernelSize-1)/2);
        meanRadDepths(ii,:,:) = tmpRadDepth((kernelSize+1)/2:size(tmpRadDepth,1)-(kernelSize-1)/2,(kernelSize+1)/2:size(tmpRadDepth,2)-(kernelSize-1)/2);    
    end
    
    B = imrotate3(meanRadDepths,-stf(i).gantryAngle,[0,0,1]);
    s = (size(B) - ct.cubeDim)/2;
    meanRadDepths = B(1+s(1):size(B,1)-s(1),1+s(2):size(B,2)-s(2),1+s(3):size(B,3)-s(3));
    
    B = imrotate3(cStdCtGrid,-stf(i).gantryAngle,[0,0,1]);
    s = (size(B) - ct.cubeDim)/2;
    cStdCtGrid = B(1+s(1):size(B,1)-s(1),1+s(2):size(B,2)-s(2),1+s(3):size(B,3)-s(3));
    

    meanRadDepthsMat{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, meanRadDepths, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');

    cStdCtGridMat{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, cStdCtGrid, ...
                                dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');

    assignin('base', 'meanRadDepths' ,meanRadDepths);
    assignin('base', 'cStdCtGrid' ,cStdCtGrid);
    
    meanRadDepths = {meanRadDepths(VctGrid)};
    cStdCtGrid = {cStdCtGrid(VctGrid)};
    radDepthVctGrid{1} = radDepthVctGrid{1}(VctGrid);
    
else
    [radDepthVctGrid, radDepthsMat] = matRad_rayTracing(stf(i),ct,VctGrid,rot_coordsV,effectiveLateralCutoff);
end

matRad_cfg.dispInfo('done.\n');

% interpolate radiological depth cube to dose grid resolution
radDepthVdoseGrid = matRad_interpRadDepth...
    (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,radDepthVctGrid);
 
if strcmp(anaMode, 'stdCorr')
    meanRadDepthVdoseGrid = matRad_interpRadDepth...
        (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,meanRadDepths);
    cStdCtGridVdoseGrid = matRad_interpRadDepth...
        (ct,1,VctGrid,VdoseGrid,dij.doseGrid.x,dij.doseGrid.y,dij.doseGrid.z,cStdCtGrid);
end

if strcmp(anaMode, 'fineSampling')
    % interpolate radiological depth cube used for fine sampling to dose grid resolution
    radDepthsMat{1} = matRad_interp3(dij.ctGrid.x,  dij.ctGrid.y,   dij.ctGrid.z, radDepthsMat{1}, ...
                                    dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z,'nearest');
end

% limit rotated coordinates to positions where ray tracing is availabe
rot_coordsVdoseGrid = rot_coordsVdoseGrid(~isnan(radDepthVdoseGrid{1}),:);
 