function [F, mcIDD, mcFWHM] = matRad_calcMCsquareObjectiveOptics(ct, stf, pln, cst, N, anaDose, ...
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
    smoothMcFWHM = smoothdata(mcFWHM, 'lowess');
    
    F = 100 * sum((anaFWHM - smoothMcFWHM).^2, 'all') / numel(anaFWHM);
    
    figure(figureHandle);
    plot(anaFWHM);
    hold on
    plot(mcFWHM);
    plot(smoothMcFWHM)
    hold off
    title(F)


    
    