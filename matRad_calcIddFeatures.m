function [peakPos, peakHeight, entryDose, r80] = matRad_calcIddFeatures(ct, IDD)

[peakHeight, peakInd] = max(IDD);
peakPos = peakInd * ct.resolution.y - ct.resolution.y/2;
[~, r80ind] = min(abs(IDD(peakInd:end) - 0.8 * peakHeight));
diff = IDD(peakInd + r80ind - 1) - IDD(peakInd + r80ind);
r80ind = r80ind + (IDD(peakInd + r80ind - 1) - 0.8 * peakHeight)  / diff;
r80 = r80ind * ct.resolution.y - ct.resolution.y/2;
entryDose = IDD(1);