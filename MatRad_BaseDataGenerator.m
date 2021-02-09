classdef MatRad_BaseDataGenerator
    
    properties
        ct
        cst
        bdl_path
    end
    
    properties (SetAccess = private)
        spacedEnergySpectrum
        spacedEmittance
    end
    
    methods
        function obj = MatRad_BaseDataGenerator(cubeDim, cubeRes, center, targetSize)
            
            % create CT
            ixOAR = 1;
            ixPTV = 2;

            % define general VOI properties
            cst{ixOAR,1} = 0;
            cst{ixOAR,2} = 'contour';
            cst{ixOAR,3} = 'OAR';

            cst{ixPTV,1} = 1;
            cst{ixPTV,2} = 'target';
            cst{ixPTV,3} = 'TARGET';

            % define optimization parameter for both VOIs
            cst{ixOAR,5}.TissueClass = 1;
            cst{ixOAR,5}.alphaX      = 0.1000;
            cst{ixOAR,5}.betaX       = 0.0500;
            cst{ixOAR,5}.Priority    = 2;
            cst{ixOAR,5}.Visible     = 1;
            cst{ixOAR,5}.visibleColor     = [1 0 0];
            cst{ixOAR,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
            cst{ixOAR,6}{1,1}.parameters{1}  = 5;
            cst{ixOAR,6}{1,1}.penalty     = 100;


            cst{ixPTV,5}.TissueClass = 1;
            cst{ixPTV,5}.alphaX      = 0.1000;
            cst{ixPTV,5}.betaX       = 0.0500;
            cst{ixPTV,5}.Priority    = 1;
            cst{ixPTV,5}.Visible     = 1;
            cst{ixPTV,5}.visibleColor     = [0 1 0];
            cst{ixPTV,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
            cst{ixPTV,6}{1,1}.parameters{1}  = 60;
            cst{ixPTV,6}{1,1}.penalty     = 800;

           

            ct.cube{1} = ones(cubeDim) * 1;
            ct.cube{1}(1,1,1) = 0; 

            ct.resolution.x = cubeRes(1);
            ct.resolution.y = cubeRes(2);
            ct.resolution.z = cubeRes(3);

            ct.cubeDim = cubeDim;

            ct.numOfCtScen  = 1;


            %% Create a cubic phantom

            % create an ct image series with zeros
            ct.cubeHU{1} = ones(ct.cubeDim) * 0;
            ct.cubeHU{1}(1,1,1) = -1000; 

            ct.hlut = [1,0;0,-1024];

            % create body of full phantom size
            mask = ones(ct.cubeDim);
            cst{1,4}{1} = find(mask == 1);

            %create target
            centerP_corr = round(center ./ cubeRes);
            height_corr = round(targetSize(1) / cubeRes(1));
            width_corr  = round(targetSize(2) / cubeRes(2));
            depth_corr  = round(targetSize(3) / cubeRes(3));


            mask = zeros(ct.cubeDim);

            for i=-ceil(height_corr/2)+1:floor(height_corr/2)
                for j=-ceil(width_corr/2):floor(width_corr/2)
                    for k=-ceil(depth_corr/2):floor(depth_corr/2)
                        coordX  = centerP_corr(1) + i;
                        coordY  = centerP_corr(2) + j;
                        coordZ  = centerP_corr(3) + k;
                        coordX(coordX<1) = 1;  
                        coordY(coordY<1) = 1;  
                        coordZ(coordZ<1) = 1; 
                        mask(coordX, coordY, coordZ) = 1;
                    end
                end
            end
            cst{2,4}{1} = find(mask == 1);
            
            obj.ct  = ct;
            obj.cst = cst;
        end   
        
        function fitMachine = fitAnalyticalDataViaFile(obj, nozzleToIso, SAD, ...
                            givenFile, fitEnergies, N)
            
            filepath = ['MCsquare\bin\BDL\' givenFile];

            pln.radiationMode   = 'protons';     % either photons / protons / carbon
            pln.machine         = 'generic_MCsquare';

            pln.propStf.bixelWidth      = 500; % [mm] / also corresponds to lateral spot spacing for particles
            pln.propStf.longitudinalSpotSpacing = 500;
            pln.propStf.gantryAngles    = 0; % [?] 
            pln.propStf.couchAngles     = 0; % [?]
            pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
            pln.propStf.isoCenter       = [obj.ct.cubeDim(2) / 2 * obj.ct.resolution.y,0, obj.ct.cubeDim(3) / 2 * obj.ct.resolution.z];
            pln.propDoseCalc.doseGrid.resolution.x = obj.ct.resolution.x; % [mm]
            pln.propDoseCalc.doseGrid.resolution.y = obj.ct.resolution.y; % [mm]
            pln.propDoseCalc.doseGrid.resolution.z = obj.ct.resolution.z; % [mm]

            % optimization settings
            pln.propOpt.optimizer       = 'IPOPT';
            pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                                  % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
            pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
            pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

            %% generate steering file
            stf = matRad_generateStf(obj.ct,obj.cst,pln);
            
            dataMC = readtable(filepath);             
            dataMC = dataMC{11:end,:};
            tmp = [];
            for i = 1:size(dataMC,1)    
                tmp = [tmp; strsplit(dataMC{i})];
            end
            energyMC    = str2double(tmp(:,1));
            spotMC      = (str2double(tmp(:,6)) + str2double(tmp(:,9)))  / 2;
            divMC       = (str2double(tmp(:,7)) + str2double(tmp(:,10))) / 2;
            corMC       = (str2double(tmp(:,8)) + str2double(tmp(:,11))) / 2;
            spread      = str2double(tmp(:,3));
            
            count = 1;
            for currentEnergy = fitEnergies
                
                mcData.divNozzle = interp1(energyMC,divMC,currentEnergy); 
                mcData.corNozzle = interp1(energyMC,corMC,currentEnergy); 
                mcData.spotNozzle = interp1(energyMC,spotMC,currentEnergy);
                mcData.z = nozzleToIso;
    
                stf.ray.energy = currentEnergy;
                resultGUI = matRad_calcDoseDirectMC(obj.ct,stf,pln,obj.cst,ones(sum(stf(:).totalNumOfBixels),1), ...
                                                    N,'fitFile.txt');
                
                tmpMachine = matRad_fitBaseData(resultGUI.physicalDose, obj.ct.resolution, currentEnergy, mcData);
                
                tmpMachine.energySpectrum.mean = currentEnergy;
                tmpMachine.energySpectrum.spread = interp1(energyMC,spread,currentEnergy);                
                tmpMachine.initFocus.Emittance.spotsize.x1 = mcData.spotNozzle;
                tmpMachine.initFocus.Emittance.spotsize.y1 = mcData.spotNozzle;
                tmpMachine.initFocus.Emittance.spotsize.x2 = mcData.spotNozzle;
                tmpMachine.initFocus.Emittance.spotsize.y2 = mcData.spotNozzle;
                tmpMachine.initFocus.Emittance.divergence.x1 = mcData.divNozzle;
                tmpMachine.initFocus.Emittance.divergence.y1 = mcData.divNozzle;
                tmpMachine.initFocus.Emittance.divergence.x2 = mcData.divNozzle;
                tmpMachine.initFocus.Emittance.divergence.y2 = mcData.divNozzle;
                tmpMachine.initFocus.Emittance.correlation.x1 = mcData.corNozzle;
                tmpMachine.initFocus.Emittance.correlation.y1 = mcData.corNozzle;
                tmpMachine.initFocus.Emittance.correlation.x2 = mcData.corNozzle;
                tmpMachine.initFocus.Emittance.correlation.y2 = mcData.corNozzle;
                tmpMachine.initFocus.Emittance.weight.first   = 1;
                tmpMachine.initFocus.Emittance.weight.second = 2;
                
                fitMachine.data(count) = tmpMachine;

                count = count + 1;
            end
            
            fitMachine.meta.radiationMode = 'protons';
            fitMachine.meta.dataType = 'doubleGauss';
            fitMachine.meta.created_on = date;
            fitMachine.meta.created_by = 'Paul Anton Meder';
            fitMachine.meta.SAD = SAD;
            fitMachine.meta.BAMStoIsoDist = nozzleToIso;
            fitMachine.meta.machine = 'Generic';
            fitMachine.meta.LUT_bxWidthminFWHM = [1, Inf; 8 ,8];
            fitMachine.meta.fitAirOffset = nozzleToIso;
        end
        
        function fitMachine = fitAnalyticalDataViaEmittance(obj, nozzleToIso, SAD, focalPoints, ...
                            energySpectrum, emittanceData, fitEnergies, N)
                        
            obj = interpolateToSpacing(obj, energySpectrum, emittanceData, fitEnergies);
            
            filepath = 'MCsquare\bin\BDL\fitFile.txt';
            writeMCsquareData(obj, nozzleToIso, focalPoints, energySpectrum, emittanceData, filepath);

            pln.radiationMode   = 'protons';     % either photons / protons / carbon
            pln.machine         = 'generic_MCsquare';

            pln.propStf.bixelWidth      = 500; % [mm] / also corresponds to lateral spot spacing for particles
            pln.propStf.longitudinalSpotSpacing = 500;
            pln.propStf.gantryAngles    = 0; % [?] 
            pln.propStf.couchAngles     = 0; % [?]
            pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
            pln.propStf.isoCenter       = [obj.ct.cubeDim(2) / 2 * obj.ct.resolution.y,0, obj.ct.cubeDim(3) / 2 * obj.ct.resolution.z];
            pln.propDoseCalc.doseGrid.resolution.x = obj.ct.resolution.x; % [mm]
            pln.propDoseCalc.doseGrid.resolution.y = obj.ct.resolution.y; % [mm]
            pln.propDoseCalc.doseGrid.resolution.z = obj.ct.resolution.z; % [mm]

            % optimization settings
            pln.propOpt.optimizer       = 'IPOPT';
            pln.propOpt.bioOptimization = 'none'; % none: physical optimization;             const_RBExD; constant RBE of 1.1;
                                                  % LEMIV_effect: effect-based optimization; LEMIV_RBExD: optimization of RBE-weighted dose
            pln.propOpt.runDAO          = false;  % 1/true: run DAO, 0/false: don't / will be ignored for particles
            pln.propOpt.runSequencing   = false;  % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below

            %% generate steering file
            stf = matRad_generateStf(obj.ct,obj.cst,pln);
            
            divMC  = [];
            corMC  = [];
            spotMC = [];
            for iX = 1:size(energySpectrum,2)
                spotMC      = [spotMC, emittanceData(iX).spotsize.x1];
                divMC       = [divMC,  emittanceData(iX).divergence.x1];
                corMC       = [corMC,  emittanceData(iX).correlation.x1];
            end
            
            count = 1;
            for currentEnergy = fitEnergies
                
                mcData.divNozzle = interp1([energySpectrum.mean],divMC,currentEnergy); 
                mcData.corNozzle = interp1([energySpectrum.mean],corMC,currentEnergy); 
                mcData.spotNozzle = interp1([energySpectrum.mean],spotMC,currentEnergy);
                mcData.z = nozzleToIso;
    
                stf.ray.energy = currentEnergy;
                resultGUI = matRad_calcDoseDirectMC(obj.ct,stf,pln,obj.cst,ones(sum(stf(:).totalNumOfBixels),1), ...
                                                    N,'fitFile.txt');
                
                tmpMachine = matRad_fitBaseData(resultGUI.physicalDose, obj.ct.resolution, currentEnergy, mcData);
                
                tmpMachine.energySpectrum = obj.spacedEnergySpectrum(count);
                tmpMachine.initFocus.Emmitance = obj.spacedEmittance(count).Emittance;
                fitMachine.data(count) = tmpMachine;
                
                count = count + 1;
            end
            
            fitMachine.meta.radiationMode = 'protons';
            fitMachine.meta.dataType = 'doubleGauss';
            fitMachine.meta.created_on = date;
            fitMachine.meta.created_by = 'Paul Anton Meder';
            fitMachine.meta.SAD = SAD;
            fitMachine.meta.BAMStoIsoDist = nozzleToIso;
            fitMachine.meta.machine = 'Generic';
            fitMachine.meta.LUT_bxWidthminFWHM = [1, Inf; 8 ,8];
            fitMachine.meta.fitAirOffset = nozzleToIso;
                
        end
        
    end
    
        
    methods (Access = protected)
       
        function obj = writeMCsquareData(obj, nozzleToIso, focalPoints, energySpectrum, emittanceData, filepath)
            %function that writes a data file containing Monte Carlo base
            %data for a simulation with MCsquare
                        
            %write MCsqaure data base file
            try
                
                fileID = fopen(filepath,'w');
                
                %Header
                %fprintf(fileID,'--matRad: Beam Model for machine %s (%s)--\n',machine.meta.machine,machine.meta.dataType);
                fprintf(fileID,'--UPenn beam model (double gaussian)--\n');
                
                fprintf(fileID,'Nozzle exit to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n', nozzleToIso);
                
                fprintf(fileID,'SMX to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n', focalPoints(1));
                
                fprintf(fileID,'SMY to Isocenter distance\n');
                fprintf(fileID,'%.1f\n\n', focalPoints(2));
                                   
                numberOfEnergies = size(energySpectrum,2);
                fprintf(fileID,'Beam parameters\n%d energies\n\n', numberOfEnergies);
                
                fprintf(fileID,'NominalEnergy	MeanEnergy	EnergySpread	ProtonsMU	Weight1	SpotSize1x	Divergence1x	Correlation1x	SpotSize1y	Divergence1y	Correlation1y	Weight2	SpotSize2x	Divergence2x	Correlation2x	SpotSize2y	Divergence2y	Correlation2y\n');	

                for k = 1:numberOfEnergies
                    fprintf(fileID, '%g ', energySpectrum(k).mean);
                    fprintf(fileID, '%g ', energySpectrum(k).mean);
                    fprintf(fileID, '%g ', energySpectrum(k).spread);
                    fprintf(fileID, '%g ', 1e+6);
                    fprintf(fileID, '%g ', emittanceData(k).weight.first);
                    fprintf(fileID, '%g ', emittanceData(k).spotsize.x1);
                    fprintf(fileID, '%g ', emittanceData(k).divergence.x1);
                    fprintf(fileID, '%g ', emittanceData(k).correlation.x1);
                    fprintf(fileID, '%g ', emittanceData(k).spotsize.y1);
                    fprintf(fileID, '%g ', emittanceData(k).divergence.y1);
                    fprintf(fileID, '%g ', emittanceData(k).correlation.y1);
                    fprintf(fileID, '%g ', emittanceData(k).weight.second);
                    fprintf(fileID, '%g ', emittanceData(k).spotsize.x2);
                    fprintf(fileID, '%g ', emittanceData(k).divergence.x2);
                    fprintf(fileID, '%g ', emittanceData(k).correlation.x2);
                    fprintf(fileID, '%g ', emittanceData(k).spotsize.y2);
                    fprintf(fileID, '%g ', emittanceData(k).divergence.y2);
                    fprintf(fileID, '%g ', emittanceData(k).correlation.y2);
                    fprintf(fileID, '\t');

                    fprintf(fileID, '\n');
                end
                
                fclose(fileID);
                
                obj.bdl_path = filepath;
                
            catch MException
                error(MException.message);
            end
        end 
           
        function obj = interpolateToSpacing(obj, energySpectrum, emmittance, fitEnergies)
            for i = 1:size(energySpectrum, 2)
                means(i)    = energySpectrum(i).mean;
                spreads(i)  = energySpectrum(i).spread;
                weight1(i)  = emmittance(i).weight.first;
                weight2(i)  = emmittance(i).weight.second;   
                spotX1(i)   = emmittance(i).spotsize.x1;
                spotY1(i)   = emmittance(i).spotsize.y1;
                divX1(i)    = emmittance(i).divergence.x1;
                divY1(i)    = emmittance(i).divergence.y1;
                corrX1(i)   = emmittance(i).correlation.x1;
                corrY1(i)   = emmittance(i).correlation.y1;
                spotX2(i)   = emmittance(i).spotsize.x2;
                spotY2(i)   = emmittance(i).spotsize.y2;
                divX2(i)    = emmittance(i).divergence.x2;
                divY2(i)    = emmittance(i).divergence.y2;
                corrX2(i)   = emmittance(i).correlation.x2;
                corrY2(i)   = emmittance(i).correlation.y2;
            end 
                
            count = 1;
            for currentEnergy = fitEnergies
                obj.spacedEnergySpectrum(count).mean   = currentEnergy;
                obj.spacedEnergySpectrum(count).spread = interp1(means, spreads, currentEnergy);
                obj.spacedEmittance(count).Emittance.weight.first = interp1(means, weight1, currentEnergy);
                obj.spacedEmittance(count).Emittance.weight.second = interp1(means, weight2, currentEnergy);
                obj.spacedEmittance(count).Emittance.spotsize.x1 = interp1(means, spotX1, currentEnergy);
                obj.spacedEmittance(count).Emittance.spotsize.x2 = interp1(means, spotX2, currentEnergy);
                obj.spacedEmittance(count).Emittance.spotsize.y1 = interp1(means, spotY1, currentEnergy);
                obj.spacedEmittance(count).Emittance.spotsize.y2 = interp1(means, spotY2, currentEnergy);
                obj.spacedEmittance(count).Emittance.divergence.x1 = interp1(means, divX1, currentEnergy);
                obj.spacedEmittance(count).Emittance.divergence.x2 = interp1(means, divX2, currentEnergy);
                obj.spacedEmittance(count).Emittance.divergence.y1 = interp1(means, divY1, currentEnergy);
                obj.spacedEmittance(count).Emittance.divergence.y2 = interp1(means, divY2, currentEnergy);
                obj.spacedEmittance(count).Emittance.correlation.x1 = interp1(means, corrX1, currentEnergy);
                obj.spacedEmittance(count).Emittance.correlation.x2 = interp1(means, corrX2, currentEnergy);
                obj.spacedEmittance(count).Emittance.correlation.y1 = interp1(means, corrY1, currentEnergy);
                obj.spacedEmittance(count).Emittance.correlation.y2 = interp1(means, corrY2, currentEnergy);
                count = count + 1;
            end
        end
    end           
end
      


