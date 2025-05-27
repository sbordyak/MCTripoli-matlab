classdef modelParameters
    %MODELPARAMETERS Model parameters - unknowns
    %   Model parameters
    
    properties
        logratios                     double = []
        logIntensity                  double = []
        refVolts                      double = []
        darkNoise                     double = []
        collectorRelativeEfficiencies double = []
        ionCounterDeadTime            double = []
        betaFaraday                   double = []
        betaDaly                      double = []
        upMassTailFaraday       (1,1) double = []
        downMassTailFaraday     (1,1) double = []
        upMassTailIC            (1,1) double = []
        downMassTailIC          (1,1) double = []
        interference                         = []
    end
    
    methods
        function obj = modelParameters(analyte, setup)
            %MODELPARAMETERS Construct an instance of this class
            %   create a new set of model parameters
            
            % For now, hard-code unknown logratios from the setup file
            % these are the subset of ratios in the analyte object that
            % have both species measured in the method file. 
            obj.logratios = log(analyte.relativeAbundances);
            obj.
        
        end
        
        function modelVector = makeModelVector(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

