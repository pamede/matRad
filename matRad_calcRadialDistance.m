function radialDist = matRad_calcRadialDistance(rotCoord, Pstart, Pend, gridX, gridY)

for j = 1:size(gridX,1)
    PstartG = [Pstart(1) + gridX(j), Pstart(2), Pstart(3) + gridY(j)];
    PendG = [Pend(1) + gridX(j), Pend(2), Pend(3) + gridY(j)];
    for i = 1:size(rotCoord,1)
        a = PstartG - PendG;
        b = rotCoord(i,:) - PendG;
        radialDist(i,j) = norm(cross(a,b)) / norm(a);
    end
end