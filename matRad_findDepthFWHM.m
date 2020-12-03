function FWHM = matRad_findDepthFWHM(ct, dose, IDD)

dose = sum(dose, 3);
laterDist = -ct.resolution.y/2 + (1:ct.cubeDim(2)) * ct.resolution.x;
zero = find(IDD < max(IDD)/25);
for i = 1:size(dose, 1)
   if any(zero == i)
       FWHM(i) = 0;
   else
       tmpProfile = dose(i,:);
       maxV = max(tmpProfile);
       tmpProfile = tmpProfile - 0.5 * maxV;  
       [~, ix] = sort(abs(tmpProfile));
       ix = sort(ix(1:4));
       w1 = interp1(tmpProfile(ix(1:2)), laterDist(ix(1:2)), 0, 'linear', 'extrap' );
       w2 = interp1(tmpProfile(ix(3:4)), laterDist(ix(3:4)), 0, 'linear', 'extrap' );
       FWHM(i) = w2 - w1;
   end
end

% FWHM = FWHM(FWHM > 0);