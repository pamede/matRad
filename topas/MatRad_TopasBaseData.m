classdef MatRad_TopasBaseData < MatRad_MCemittanceBaseData
    % MatRad_TopasBaseData class for calculating TOPAS base data and
    % writing it into a file, formatted for TOPAS to use
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2019 the matRad development team.
    %
    % This file is part of the matRad project. It is subject to the license
    % terms in the LICENSE file found in the top-level directory of this
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
    % of the matRad project, including this file, may be copied, modified,
    % propagated, or distributed except according to the terms contained in the
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        
    end
    
    methods
        function obj = MatRad_TopasBaseData(varargin)
            % set default values for the simulation, if not specified
            if  nargin < 2 % input MatRad_TopasBaseData(machine)
                TopasConfig = struct;  
                constArguments = 1;
            elseif nargin < 3 && isfield(varargin{2}, 'gantryAngle') % input MatRad_TopasBaseData(machine,stf)
                TopasConfig = struct;
                constArguments = [1,2];
            else % input MatRad_TopasBaseData(machine,T) || (machine,stf,T)
                TopasConfig = varargin{end};
                if nargin == 2
                    constArguments = 1;
                elseif nargin == 3
                    constArguments = [1,2];
                end
            end
            
            % Call MatRad_MCemmitanceBaseData constructor
            obj = obj@MatRad_MCemittanceBaseData(varargin{constArguments});           
        end
        
        function obj = uncQuant(obj, uncQuantParam)
            % alter exisitning TOPAS MC parameterization for uncertainty
            % Quantification, uncQuantParam(1) = errorEnergy is relative 
            % (0 to 1), uncQuantParam(2) = errorSigmaXY (is added in 
            % quadrature), uncQuantParam(3) = beam divergence (optional),
            % uncQuantParam(4) = beam correlation (optional)
            
            % assign default values to beam divergence and correlation
            if size(uncQuantParam,2) < 2 || size(uncQuantParam,2) > 4
                error('Wrong number of parameters for uncertainty Quantification!'); 
            elseif size(uncQuantParam,2) < 3
                errorEnergy = uncQuantParam(1);
                errorSigmaXY = uncQuantParam(2);
                div = 0;
                corr = 0;
            elseif size(uncQuantParam,2) < 4
                errorEnergy = uncQuantParam(1);
                errorSigmaXY = uncQuantParam(2);
                div = uncQuantParam(3);
                corr = 0;
            else
                errorEnergy = uncQuantParam(1);
                errorSigmaXY = uncQuantParam(2);
                div = uncQuantParam(3);
                corr = uncQuantParam(4);
            end
            
            
            %  alter entries for energy spread as well as emittance
            %  parameterizeation according to input parameters
            for i = 1:size(obj.monteCarloData,2)
                
                if errorEnergy < 0 || errorEnergy > 1
                % check whether energy error is reasonable

                	error('Selected invalid energy error in uncertainty Quantification!'); 
                else
                    obj.monteCarloData(i).EnergySpread = obj.monteCarloData(i).MeanEnergy * errorEnergy;
                end
                
                obj.monteCarloData(i).SpotSize1x    = sqrt(obj.monteCarloData(i).SpotSize1x^2 + errorSigmaXY^2);
                obj.monteCarloData(i).SpotSize1y    = sqrt(obj.monteCarloData(i).SpotSize1y^2 + errorSigmaXY^2);
                obj.monteCarloData(i).Divergence1x  = div;
                obj.monteCarloData(i).Divergence1y  = div;
                obj.monteCarloData(i).Correlation1x = corr;
                obj.monteCarloData(i).Correlation1y = corr;
            end
        end
        
        function obj = writeTopasData(obj,ct,stf,pln,w)
            %function that writes a data file containing stf specific data
            %for a Monte Carlo simulation with TOPAS
                     
            %look up focus indices
            focusIndex = obj.selectedFocus(obj.energyIndex);
            
            % NozzleAxialDistance
            nozzleAxialDistance_mm = 1500;
            if isfield( obj.machine.meta,'nozzleAxialDistance')
                disp('Using NAD from basedata')
                nozzleAxialDistance_mm =  obj.machine.meta.nozzleAxialDistance;
            elseif isfield(obj.machine.meta,'BAMStoIsoDist')
                disp('Using BAMStoIsoDist from basedata')
                nozzleAxialDistance_mm = obj.machine.meta.BAMStoIsoDist;
            else
                disp('Using default nozzleAxialDistance')
            end
            
            
        end
        
    end
end

