function [gridX, gridY] = matRad_createFineSamplingGrid(stf, gridsize)
 
margin = [20, 20];
positions = [stf.ray(:).rayPos_bev];
xLimits = [min(positions(1:3:end))-margin(1), max(positions(1:3:end))+margin(1)];
yLimits = [min(positions(3:3:end))-margin(2), max(positions(3:3:end))+margin(1)];

xValues = [fliplr(0:-gridsize(1):xLimits(1)), gridsize(1):gridsize(1):xLimits(2)];
yValues = [fliplr(0:-gridsize(2):yLimits(1)), gridsize(2):gridsize(2):yLimits(2)];

[gridX, gridY] = meshgrid(xValues, yValues);
gridX = reshape(gridX, numel(gridX), 1);
gridY = reshape(gridY, numel(gridY), 1);

end