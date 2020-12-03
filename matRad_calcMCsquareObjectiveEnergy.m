function [F, mcIDD] = matRad_calcMCsquareObjectiveEnergy(ct, stf, pln, cst, N, anaDose, mean, spread, figureHandle)

    anaIDD = sum(anaDose,[2 3]);
    anaIDD = anaIDD(anaIDD > 0);
    
    optMode = 'energy';
    resultGUI = matRad_calcDoseDirectMCvarE(optMode, [mean, spread], ct,stf,pln,cst,ones(sum([stf(:).totalNumOfBixels]),1), N);
    mcIDD = sum(resultGUI.physicalDose,[2 3]); 
    mcIDD = mcIDD(1:numel(anaIDD));
    
    
    F = 500 * sum((anaIDD - mcIDD).^2, 'all') / numel(anaIDD) * 1e5;

    figure(figureHandle);
    plot(anaIDD);
    hold on
    plot(mcIDD);
    hold off
    title(F)