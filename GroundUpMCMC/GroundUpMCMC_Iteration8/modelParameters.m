classdef modelParameters
    %MODELPARAMETERS Model parameters - unknowns
    %   Model parameters
    
    properties
        logratios                     double = []
        logIntensityKnots             double = []
        refIntensities                double = []
        darkNoise                     double = []
        collectorRelativeEfficiencies double = []
        ionCounterDeadTime            double = []
        betaFaraday                   double = []
        betaDaly                      double = []
        upMassTailFaraday             double = []
        downMassTailFaraday           double = []
        upMassTailIC                  double = []
        downMassTailIC                double = []
        interference                         = []
    end
    
    methods
        function obj = modelParameters(analyte, setup)
            %MODELPARAMETERS Construct an instance of this class
            %   create a new set of model parameters
            
            % For now, hard-code unknown logratios from the setup file
            % these are the subset of ratios in the analyte object that
            % have both species measured in the method file. 
            obj.logratios = log( ...
                analyte.relativeAbundances(4)/analyte.relativeAbundances(2) ...
                );
            obj.logIntensityKnots = setup.logmspl;
            obj.refIntensities = setup.refIntensities;
        
        end % constructor
        
        function modelVector = makeModelVector(obj)
            %MAKEMODELVECTOR Assemble model vector from unknown parameters
            %   This order hard-coded. Labels for each entry in modelLabels
            modelVector = [
                obj.logratios; 
                obj.logIntensityKnots;
                obj.refIntensities;
                obj.collectorRelativeEfficiencies;
                obj.betaFaraday;
                obj.betaDaly;
                obj.upMassTailFaraday;
                obj.downMassTailFaraday;
                obj.darkNoise;
                obj.upMassTailIC;
                obj.downMassTailIC
                ];
        end % function modelVector

        function modelLabels = makeModelLables(obj)
            %MAKEMODELVECTOR Assemble model vector from unknown parameters
            %   This order hard-coded. Labels for each entry in modelLabels
            modelLabels = [
                "obj.logratios; 
                obj.logIntensityKnots;
                obj.refIntensities;
                obj.collectorRelativeEfficiencies;
                obj.betaFaraday;
                obj.betaDaly;
                obj.upMassTailFaraday;
                obj.downMassTailFaraday;
                obj.darkNoise;
                obj.upMassTailIC;
                obj.downMassTailIC;
                ];
        end 

    end
end

