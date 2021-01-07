function [F, mcIDD, mcFWHM] = matRad_calcMCsquareObjectiveAllXY(ct, stf, pln, cst, N, anaDose, ...
                                            mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY, foundFocalX, foundFocalY, figureHandle)
    
 
    anaIDD = sum(anaDose,[2 3]);
    anaFWHM1 = matRad_findDepthFWHM(ct, anaDose, anaIDD, 1);
    anaFWHM2 = matRad_findDepthFWHM(ct, anaDose, anaIDD, 2);

    [~, ixValidFWHM] = max(anaFWHM1);
    anaFWHM1 = anaFWHM1(1:ixValidFWHM);
    anaFWHM2 = anaFWHM2(1:ixValidFWHM);
    
    optMode = 'all';
    resultGUI = matRad_calcDoseDirectMCvarE(optMode, [mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY, foundFocalX, foundFocalY], ...
                                    ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), N);
    mcDose = resultGUI.physicalDose;
    mcIDD = sum(resultGUI.physicalDose,[2 3]); 
    mcFWHM1 = matRad_findDepthFWHM(ct, mcDose, mcIDD, 1);
    mcFWHM2 = matRad_findDepthFWHM(ct, mcDose, mcIDD, 2);

    mcFWHM1 = mcFWHM1(1:ixValidFWHM);
    mcFWHM2 = mcFWHM2(1:ixValidFWHM);
    
    
    F11 = 100 * sum((anaFWHM1 - mcFWHM1).^2, 'all') / numel(anaFWHM1);
    F12 = 100 * sum((anaFWHM2 - mcFWHM2).^2, 'all') / numel(anaFWHM2);

    F1 = F11 + F12;
    
    F2 = 5 * sum((anaIDD - mcIDD).^2, 'all') / numel(anaIDD) * 1e5;
    
    figure(figureHandle);
    subplot(1,2,1)
    plot(anaFWHM1, 'Color', 'k','LineStyle', '-')
    hold on
    plot(anaFWHM2, 'Color', 'b','LineStyle', '-')

    plot(mcFWHM1, 'Color', 'k','LineStyle', '--')
    plot(mcFWHM2, 'Color', 'b','LineStyle', '--')
    hold off
    title(F1)
    
    subplot(1,2,2)
    plot(anaIDD);
    hold on
    plot(mcIDD);
    hold off
    title(F2)
   
    F = 100 * (F1 + F2);

    
    