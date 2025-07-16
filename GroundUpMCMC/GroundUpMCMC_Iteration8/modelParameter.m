classdef modelParameter
    %MODELPARAMETER Unknown or assumed model parameter
    %   isFree = true for parameters free to vary in model
    %   isFree = false for parameters with an assumed value
    %   value is the vector of values

    properties
        isFree     (:,1) logical = false
        value      (:,1) double  = []
    end

    methods
        function obj = modelParameter(value, isFree, indexInSet)
            %MODELPARAMETER Construct an instance of this class
            %   Create a model parameter
            
            if nargin > 0 
                obj.isFree = isFree;
                obj.value = value;
            elseif nargin == 0 % if no arguments, create object with defaults
                obj.isFree = false;
                obj.value = [];
            end % if nargin > 0

        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
    end
end