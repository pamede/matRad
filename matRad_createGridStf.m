function newStf = matRad_createGridStf(ct, stf, gridX, gridY, gridsize)
 
newStf.gantryAngle     = stf.gantryAngle;
newStf.couchAngle      = stf.couchAngle;
newStf.bixelWidth      = stf.bixelWidth;
newStf.radiationMode   = stf.radiationMode;
newStf.SAD             = stf.SAD;
newStf.isoCenter       = stf.isoCenter;



newStf.numOfRays = size(gridX, 1);





energies = sort(unique([stf.ray(:).energy]));

for i = 1:size(gridX, 1)
    tmpRay.rayPos_bev = [gridX(i) 0 gridY(i)];
    tmpRay.targetPoint_bev = [2 * gridX(i) stf.SAD 2 * gridY(i)];
    
    rotMat_system_T = matRad_getRotationMatrix(stf.gantryAngle,stf.couchAngle);
    tmpRay.rayPos = tmpRay.rayPos_bev * rotMat_system_T;
    tmpRay.targetPoint = tmpRay.targetPoint_bev * rotMat_system_T;
    tmpRay.energy = energies;
    tmpRay.focusIx = ones(size(tmpRay.energy));
    for j = 1:size(energies,2)
        tmpRay.rangeShifter(j).ID = 0;
        tmpRay.rangeShifter(j).eqThickness = 0;
        tmpRay.rangeShifter(j).sourceRashiDistance = 0;
    end
    newStf.ray(i) = tmpRay;
end

newStf.sourcePoint_bev = stf.sourcePoint_bev;
newStf.sourcePoint     = stf.sourcePoint;
newStf.numOfBixelsPerRay = repmat(size(energies, 2), size(gridX, 1),1)';
newStf.longitudinalSpotSpacing = stf.longitudinalSpotSpacing;
newStf.totalNumOfBixels = size(energies, 2) * size(gridX, 1);
newStf = matRad_computeSSD(newStf, ct);

end