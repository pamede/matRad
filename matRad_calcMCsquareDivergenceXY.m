function [F, mcIDD, mcFWHM] = matRad_calcMCsquareDivergenceXY(ct, stf, pln, cst, N, anaDose, ...
                                            mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY, figureHandle)
    
                                        
                                            
    anaIDD = sum(anaDose,[2 3]);
    anaFWHM = matRad_findDepthFWHM(ct, anaDose, anaIDD);
    [~, ixValidFWHM] = max(anaFWHM);
    anaFWHM = anaFWHM(1:ixValidFWHM);

    optMode = 'all';
    resultGUI = matRad_calcDoseDirectMCvarE(optMode, [mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY], ...
                                    ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), N);
    mcDose = resultGUI.physicalDose;
    mcIDD = sum(resultGUI.physicalDose,[2 3]); 
    mcFWHM = matRad_findDepthFWHM(ct, mcDose, mcIDD);
    mcFWHM = mcFWHM(1:ixValidFWHM);  
    
    anaDiff = diff(anaFWHM);    
    mcDiff  = diff(mcFWHM);
    
    
    
    
    [~, ixValidDiff] = max(anaDiff);
    anaDiff = anaDiff(1:ixValidDiff);
    mcDiff = mcDiff(1:ixValidDiff);

    anaDiffS = smoothdata(anaDiff, 'movmean');
    mcDiffS  = smoothdata(smoothdata(mcDiff, 'movmean', 10) ,'movmean', 10);
    
    F = 1e5 * sum((anaDiffS - mcDiffS).^2, 'all') / numel(anaFWHM);

%     figure(figureHandle);
    plot(anaDiff, 'Color', 'k','LineStyle', '--')
    hold on
    plot(anaDiffS, 'Color', 'k', 'LineStyle', '-')
    plot(mcDiff, 'Color', 'b','LineStyle', '--')
    plot(mcDiffS, 'Color', 'b', 'LineStyle', '-')
    hold off
    title(F)



    
    