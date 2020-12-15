function [F, mcIDD, mcFWHM] = matRad_calcMCsquareObjectiveOpticsXY(ct, stf, pln, cst, N, anaDose, ...
                                            mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY, figureHandle)
    
                                        
                                            
    anaIDD = sum(anaDose,[2 3]);
    anaFWHM1 = matRad_findDepthFWHM(ct, anaDose, anaIDD, 1);
    anaFWHM2 = matRad_findDepthFWHM(ct, anaDose, anaIDD, 2);

    [~, ixValidFWHM] = max(anaFWHM1);
    anaFWHM1 = anaFWHM1(1:ixValidFWHM);
    anaFWHM2 = anaFWHM2(1:ixValidFWHM);
    
    optMode = 'all';
    resultGUI = matRad_calcDoseDirectMCvarE(optMode, [mean, spread, spotsizeX, divergenceX, correlationX, spotsizeY, divergenceY, correlationY], ...
                                    ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), N);
    mcDose = resultGUI.physicalDose;
    mcIDD = sum(resultGUI.physicalDose,[2 3]); 
    mcFWHM1 = matRad_findDepthFWHM(ct, mcDose, mcIDD, 1);
    mcFWHM2 = matRad_findDepthFWHM(ct, mcDose, mcIDD, 2);

    mcFWHM1 = mcFWHM1(1:ixValidFWHM);
    mcFWHM2 = mcFWHM2(1:ixValidFWHM);

%     smoothMcFWHM = smoothdata(mcFWHM, 'lowess');
    
%     F = 100 * sum((anaFWHM - smoothMcFWHM).^2, 'all') / numel(anaFWHM);
    F1 = 100 * sum((anaFWHM1 - mcFWHM1).^4, 'all') / numel(anaFWHM1);
    F2 = 100 * sum((anaFWHM2 - mcFWHM2).^4, 'all') / numel(anaFWHM2);

    F = F1 + F2;
    
    figure(figureHandle);
    plot(anaFWHM1, 'Color', 'k','LineStyle', '-')
    hold on
    plot(anaFWHM2, 'Color', 'b','LineStyle', '-')

    plot(mcFWHM1, 'Color', 'k','LineStyle', '--')
    plot(mcFWHM2, 'Color', 'b','LineStyle', '--')
    hold off
    title(F)


    
    