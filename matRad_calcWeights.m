function [finalWeight, sigmaBeamlet, posX, posY, numOfSub] = matRad_calcWeights(sigmaTot, sigmaSub, N)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates weights for a fine sampling pencil beam
% algorithm
%
% call
%   [finalWeight, sigmaBeamlet, posX, posY, numOfSub] =
%                               matRad_calcWeights(sigmaTot, sigmaSub, N)
%
% input
%   sigmaTot:       the standard deviation of the lateral spread of the pencil
%                   beam
%   sigmaSub:       is the sigma of the gaussian of the sub-beams
%   N:              number of subsample beams shells. Means we have a
%                   grid of NxN sub beams representing the total beam
%
% output
%   finalWeight:    is the array of the weights of the sub-pencil beams. It
%                   runs over the same index as posx and posy
%   
%   posx & posy:    are the positions of the sub-pencil beams, returned as
%                   meshgrid-matrices if method is 'square' and as vectors 
%                   if method is 'circle'
%   numOfSub:       number of sub-pencil beams
%
% References
%   [1] https://iopscience.iop.org/article/10.1088/0031-9155/61/1/183
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);

gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

R = 3.5 * sigmaHead;
dR = 2 * R / N;

counter = 1;
for iX = -(N - 1) / 2 : (N - 1) / 2
    for iY = -(N - 1) / 2 : (N - 1) / 2
        finalWeight(counter) = integral2(@(x,y) gauss(sigmaHead, x, y, 0, 0), ...
                        (iX - 0.5) * dR, (iX + 0.5) * dR, ...
                    	(iY - 0.5) * dR, (iY + 0.5) * dR);
        posX(counter) = iX * dR;
        posY(counter) = iY * dR;
        sigmaBeamlet(counter) = sigmaSub;

        counter = counter + 1;
    end
end

finalWeight = finalWeight';
posX = posX';
posY = posY';
sigmaBeamlet = sigmaBeamlet';
        
        
numOfSub = N * N;


