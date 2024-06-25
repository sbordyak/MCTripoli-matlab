classdef ionBeamStats
    %IONBEAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        integrationTime
        resistorOhms
        intensity struct
        PoissonVarianceInCPS
        totalVarianceInVolts
    end

    properties (Constant)
        kB = 1.380649e-23; % exact, Joules per Kelvin
        T = 290; % temperature, Kelvin
        ionsPerCoulomb = 6241509074460762607.776;
    end
    
    methods
        function obj = ionBeam(integrationTime,resistorOhms,intensity)
            %IONBEAM Construct an instance of this class
            %   Detailed explanation goes here
            obj.integrationTime = integrationTime;
            obj.resistorOhms = resistorOhms;
            obj.intensity = intensity;
            
            obj.PoissonVarianceInCPS = integrationTime * intensity.cps;
        end
        
        function intensity = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

