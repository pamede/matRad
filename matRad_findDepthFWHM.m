function FWHM = matRad_findDepthFWHM(ct, dose, IDD)

dose = sum(dose, 3);
zero = find(IDD < max(IDD)/25);
for i = 1:size(dose, 1)
   if any(zero == i)
       FWHM(i) = 0;
%        FWHM1(i) = 0;
   else
       tmpProfile = dose(i,:);
       laterDist = -ct.resolution.y/2 + (1:ct.cubeDim(2)) * ct.resolution.x;

%        plot(tmpProfile)
       
       [maxV, maxIx]= max(tmpProfile);

        profileScale = tmpProfile./maxV;
        nonZero = find(profileScale);
        profileScale = profileScale(nonZero);
        laterDist = laterDist(nonZero);
       
        maxi = find(profileScale == 1);
        
        left    = profileScale(1:maxi);
        leftx   = laterDist(1:maxi);
        right   = profileScale(maxi+1:end);
        rightx  = laterDist(maxi+1:end);

        sp1 = interp1(left,   leftx,  0.5);
        sp2 = interp1(right,  rightx, 0.5);
        FWHM(i) = sp2-sp1;
%         
%         tmpProfile = tmpProfile - 0.5 * maxV; 
% 
%         [~, ix] = sort(abs(tmpProfile));
%         ix = sort(ix(1:4))
%         w1 = interp1(tmpProfile(ix(1:2)), laterDist(ix(1:2)), 0, 'linear', 'extrap' );
%         w2 = interp1(tmpProfile(ix(3:4)), laterDist(ix(3:4)), 0, 'linear', 'extrap' );
%         FWHM1(i) = w2 - w1;
   end
end
%     [~, ixValidFWHM] = max(FWHM);
%     FWHM = FWHM(1:ixValidFWHM);
    
%     plot(FWHM);
%     hold on
%     plot(FWHM1);

%     plot(smoothdata(FWHM, 'lowess'))
%     hold off

% FWHM = FWHM(FWHM > 0);