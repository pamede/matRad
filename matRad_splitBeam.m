function newStf = matRad_splitBeam(stf, machine)

newStf = stf;
N = 41;

for i = 1:stf.numOfRays
    tmpRay = stf.ray(i);
    initSigma = matRad_calcInitSigma(tmpRay, stf.ray(i).SSD, machine);
    
    [weights, xPos, yPos] = matRad_calcWeights(initSigma, N, 3);
    
    newRays = [];
    for j = 1:size(xPos,1)
        newRay = tmpRay;
        newRay.rayPos_bev = [ newRay.rayPos_bev(1) + xPos(j) ...
                                newRay.rayPos_bev(2) ...           
                                newRay.rayPos_bev(3) + yPos(j)];
        newRay.targetPoint_bev = [ 2*newRay.rayPos_bev(1) ...
                                    stf.SAD ...
                                    2*newRay.rayPos_bev(3)];
        newRay.rayPos = newRay.rayPos_bev;
        newRay.targetPoint = newRay.targetPoint_bev; %only correct when angles are zero
        newRay.weight = weights(j);
        newRays = [newRays, newRay];
    end

    newStf.ray = newRays;
    newStf.numOfRays = size(newRays,2);
    newStf.numOfBixelsPerRay = ones(1,size(newRays,2));
    newStf.totalNumOfBixels = sum(newStf.numOfBixelsPerRay);
end
