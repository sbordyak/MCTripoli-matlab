classdef synDataSetup
    %SYNDATASETUP Metadata and methods for synthetic data
    %   Contains all required information to define a synthetic dataset and
    %   class methods to create synthetic datasets

    properties
        
        modelParameters modelParameterSet
        nBlocks         uint32
        BLTimes         (:,1) double
        OPTimes         (:,1) double

    end

    properties (SetAccess = immutable)

        massSpec    massSpecModel
        method      struct

    end

    methods

        function obj = synDataSetup(synDataSetupName)
            %SYNDATASETUP Construct an instance of this class
            %   Hard-code synthetic data parameter properties by name

            switch synDataSetupName

                case "Synthetic_TwoIsotopeStatic_1Mcps_Listric"
                % first synthetic data example, from GroundUp7b
                
                % model parameter set
                NBS982 = referenceMaterial("NBS982");
                mTrue = modelParameterSet();
                mTrue.logratios.value = NBS982.logRatioAbundances(4);
                mTrue.logratios.isFree = 1;
                mTrue.logIntensityKnots.value = [ % nseg = 2, bdeg = 3, t = 106:205
                    14.4846439892590;
                    14.5251376191125;
                    14.4721275963024;
                    14.4292231240153;
                    14.4531727595166];
                mTrue.logIntensityKnots.isFree = true(5, 1);
                mTrue.refIntensities.value = [-1e2; 2e2];
                mTrue.refIntensities.isFree = true(2, 1);
                mTrue.betaFaraday.value = -0.2;
                mTrue.betaFaraday.isFree = false;
                mTrue.upMassTailFaraday.value = 0;
                mTrue.upMassTailFaraday.isFree = false;
                mTrue.downMassTailFaraday.value = 0;
                mTrue.downMassTailFaraday.isFree = false;
                obj.modelParameters = mTrue;

                % other parameters needed for synthetic data generation
                obj.nBlocks = 1;
                obj.BLTimes = 1:100;
                obj.OPTimes = 106:205;

                % massSpec
                obj.massSpec = massSpecModel("PhoenixKansas_1e12");

                % make method struct from TIMSAM and massSpec collector names
                methodsFolder = "./massSpecMethods/";
                methodName = "TwoCollectorStaticPb";
                method = parseTIMSAM(methodsFolder + methodName + ".TIMSAM");
                method = processMethod(method, obj.massSpec.getCollectorNames);
                obj.method = method;


            end % switch synDataSetupName
            



        end

        function dataset = makeSynData(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end