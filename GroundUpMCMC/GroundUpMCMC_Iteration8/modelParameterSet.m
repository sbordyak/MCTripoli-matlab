classdef modelParameterSet
    %MODELPARAMETERS Model parameters - unknowns
    %   Model parameters
    
    properties
        logratios               (:,1) modelParameter
        logIntensityKnots       (:,1) modelParameter
        refIntensities          (:,1) modelParameter
        darkNoise               (:,1) modelParameter
        collectorRelativeEfficiencies (:,1) modelParameter
        ionCounterDeadTime      (:,1) modelParameter
        betaFaraday             (1,1) modelParameter
        betaDaly                (1,1) modelParameter
        upMassTailFaraday       (1,1) modelParameter
        downMassTailFaraday     (1,1) modelParameter
        upMassTailIC            (1,1) modelParameter
        downMassTailIC          (1,1) modelParameter
        interference            (:,1) modelParameter
    end
    
    methods
        function obj = modelParameterSet(analyte, setup)
            %MODELPARAMETERS Construct an instance of this class
            %   create a new set of model parameters
            
            if nargin > 0 % if no arguments, use blank defaults
            % For now, hard-code unknown logratios from the setup file
            % these are the subset of ratios in the analyte object that
            % have both species measured in the method file. 
            
            obj.logratios(1).value = log( ...
                analyte.relativeAbundances(4)/analyte.relativeAbundances(2) ...
                );
            obj.logIntensityKnots(:).value = setup.logmspl;
            obj.refIntensities(:).value = setup.refIntensities;
            end
        
        end % constructor
        
        function modelVector = modelVector(obj)
            %MAKEMODELVECTOR Assemble model vector from unknown parameters
            %   This order hard-coded. Labels for each entry in modelLabels
            modelVector = [
                obj.logratios.value; 
                obj.logIntensityKnots.value;
                obj.refIntensities.value;
                obj.collectorRelativeEfficiencies.value;
                obj.betaFaraday.value;
                obj.betaDaly.value;
                obj.upMassTailFaraday.value;
                obj.downMassTailFaraday.value;
                obj.darkNoise.value;
                obj.upMassTailIC.value;
                obj.downMassTailIC.value
                ];
        end % function modelVector

        function modelLabels = makeModelLabels(obj)
            %MAKEMODELVECTOR Assemble model vector from unknown parameters
            %   This order hard-coded. Labels for each entry in modelLabels
            modelLabels = [
                repelem("logratios", length(obj.logratios))'; 
                repelem("logIntensityKnots", length(obj.logIntensityKnots))';
                repelem("refIntensities", length(obj.refIntensities))';
                repelem("collectorRelativeEfficiencies", length(obj.collectorRelativeEfficiencies))';
                repelem("betaFaraday", length(obj.betaFaraday))';
                repelem("betaDaly", length(obj.betaDaly))';
                repelem("upMassTailFaraday", length(obj.upMassTailFaraday))';
                repelem("downMassTailFaraday", length(obj.downMassTailFaraday))';
                repelem("darkNoise", length(obj.darkNoise))';
                repelem("upMassTailIC", length(obj.upMassTailIC))';
                repelem("downMassTailIC", length(obj.downMassTailIC))'; 
                ];
        end % function modelLabels

    end
end

