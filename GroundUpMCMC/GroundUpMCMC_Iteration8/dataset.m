classdef dataset
    %DATASET Measured mass spectrometer data
    %   Includes constructor and vectorizing methods
    
    properties
        intensity       (:,1) double  = []
        time            (:,1) double  = []
        integrationTime (:,1) double  = []
        detectorIndex   (:,1) uint8   = []
        sequenceIndex   (:,1) uint8   = []
        mass            (:,1) double  = []
        isotope         (:,1) uint8   = []
        block           (:,1) uint16  = []
        cycle           (:,1) uint32  = []
        isOP            (:,1) logical = []
        isUsed          (:,1) logical = []
    end
    
    methods
        function obj = dataset(varargin)
            
        if nargin > 0 % if constructred with inputs
            data = cell2struct(varargin(1));
            method = cell2struct(varargin(2));
            % Scott: replicate assembleDataVector(data, method) here
            % to create a dataset

            

        end % if constructed from a data input

        end
       
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

