function initSigma = matRad_calcInitSigma(ray, SSD, machine)

[~ ,energyIx, ~] = intersect([machine.data(:).energy],ray.energy);
initSigma = [];
for i = 1:size(energyIx,1)
initSigma = [initSigma, matRad_interp1(machine.data(energyIx(i)).initFocus.dist(ray.focusIx(i),:)', machine.data(energyIx(i)).initFocus.sigma(ray.focusIx(i),:)',SSD)];
end

end