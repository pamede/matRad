function [F, mcIDD, mcFWHM] = matRad_calcMCsquareObjectiveAll(ct, stf, pln, cst, N, anaDose, ...
                                            mean, spread, spotsize, divergence, correlation, figureHandle)
    
    anaIDD = sum(anaDose,[2 3]);
    anaFWHM = matRad_findDepthFWHM(ct, anaDose, anaIDD);
    [~, ixValidFWHM] = max(anaFWHM);
    anaFWHM = anaFWHM(1:ixValidFWHM);

    optMode = 'all';
    resultGUI = matRad_calcDoseDirectMCvarE(optMode, [mean, spread, spotsize, divergence, correlation], ...
                                    ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), N);
    mcDose = resultGUI.physicalDose;
    mcIDD = sum(resultGUI.physicalDose,[2 3]); 
    mcFWHM = matRad_findDepthFWHM(ct, mcDose, mcIDD);
    mcFWHM = mcFWHM(1:ixValidFWHM);
    smoothMcFWHM = smoothdata(mcFWHM, 'lowess', 'SmoothingFactor', 0.35);
    
    F1 = sum((anaFWHM - mcFWHM).^2, 'all') / numel(anaFWHM);
    F2 = 5 * sum((anaIDD - mcIDD).^2, 'all') / numel(anaIDD) * 1e5;
    
    figure(figureHandle);
    subplot(1,2,1)
    plot(anaFWHM);
    hold on
    plot(mcFWHM);
%     plot(smoothMcFWHM)
    hold off
    title(F1)
    
    subplot(1,2,2)
    plot(anaIDD);
    hold on
    plot(mcIDD);
    hold off
    title(F2)
   
    F = 100 * (F1 + F2);

    
    