function radialDistSq = matRad_calcRadialDistanceSq(rotCoord, Pstart, Pend, gridX, gridY)
% PstartG = [Pstart(1) + gridX, Pstart(2) * ones(size(gridX)), Pstart(3) + gridY];
% PendG =   [Pend(1)   + gridX, Pend(2)   * ones(size(gridX)), Pend(3)   + gridY];

% a = Pstart - Pend;
% b = rotCoord - Pend;
% 
% C = vecnorm(cross(repmat(a,size(b,1),1),b,2),2,2) / norm(a);

% radialDist(i,j) = norm(cross(,b)) / norm(a);

a = Pstart - Pend;
normA = norm(a);


PstartG = [Pstart(1) + gridX, Pstart(2) * ones(size(gridX)), Pstart(3) + gridY];
PendG =   [Pend(1)   + gridX, Pend(2)   * ones(size(gridX)), Pend(3)   + gridY];


b1 = rotCoord(:,1) - PendG(:,1)';
b2 = rotCoord(:,2) - PendG(:,2)';
b3 = rotCoord(:,3) - PendG(:,3)';


% tic
% a2b2 = abs(normA) * (b1^2 + b2^2 + b3^2);
% ab = a(1) * b1 + a(2) * b2 + a(3) * b3;
% cross1 = 
% toc


cross1 = (a(2) * b3 - a(3)*b2).^2;
cross2 = (a(3) * b1 - a(1)*b3).^2;
cross3 = (a(1) * b2 - a(2)*b1).^2;



radialDistSq = (cross1+cross2+cross3) / (normA^2);

% tic
% b1 = bsxfun(@minus, rotCoord(:,1), PendG(:,1)');
% b2 = bsxfun(@minus, rotCoord(:,2), PendG(:,2)');
% b3 = bsxfun(@minus, rotCoord(:,3), PendG(:,3)');
% toc

% % z = zeros(123066,441,3);
% % z(:,:,1) = cross1;
% % z(:,:,2) = cross2;
% % z(:,:,3) = cross3;
% % 
% % abs = vecnorm(z, 2, 3);
% % 
% % 
% % 
% % t = repmat(reshape(a,1,1,3), size(b,1),size(b,2),1);
% % 
% % 
% % % norm cross / norm(a)
% % 
% % tmp = vecnorm(cross(repmat(a,size(b,1),1),b,2),2,2) / norm(a);
% 
% 
% 
% for j = 1:size(gridX,1)
%     PstartG = [Pstart(1) + gridX(j), Pstart(2), Pstart(3) + gridY(j)];
%     PendG = [Pend(1) + gridX(j), Pend(2), Pend(3) + gridY(j)];
%     radialDist(:,j) = radDistTest(rotCoord, PstartG, PendG);
% end