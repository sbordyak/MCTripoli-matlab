classdef modelParameters
    %MODELPARAMETERS Model parameters - unknowns
    %   Model parameters
    
    properties
        logratios    
        logIntensity 
        refVolts     
        darkNoise
        collectorRelativeEfficiencies
        ionCounterDeadTime
        betaFaraday
        betaDaly
        upMassTailFaraday
        downMassTailFaraday
        upMassTailIC
        downMassTailIC
        interference
    end
    
    methods
        function obj = modelParameters(analyte, setup)
            %MODELPARAMETERS Construct an instance of this class
            %   create a new set of model parameters
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function modelVector = makeModelVector(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

