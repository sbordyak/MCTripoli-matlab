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

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)

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
            



        end % constructor method

        function dataset = makeSynData(obj)
            %MAKESYNDATA Make synthetic dataset (ie data vectors)
            %   Uses massSpec, method, and true modelParameterSet
            
            arguments (Output)
                dataset (1,1) dataVectors
            end % arguments

            dataset = obj.method;

        end % makeSynData method 


    end % public methods block


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)


    end % protected methods block


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Access = protected)

        function ionBeam = simulateIonBeam(trueCountRates, integrationTimes, detector)
            %SIMULATEIONBEAM Create synthetic measured ion beam
            %   This simple, pared-down version returns a simulated ion beam.
            %   Noise is the sum shot noise and, if Farday output is selected,
            %   Johnson-Nyquist noise.  No dead time correction.
            %
            %   Inputs:
            %   - trueCountRates: vector of true count rates in cps
            %   - integrationTimes: vector of integration times in cps
            %   - detector: struct containing fields
            %       - type = "F" or "IC" for "Faraday" or "Ion Counter"
            %       - resistance = eg 1e11 or 1e12
            %       - gain = eg 1 or 0.9, fraction of trueCounts measured by detector
            %
            %   Output:
            %   - ionBeam: vector of measured intensities, with
            %              units of cps if type = "F", cps if type = "IC"
            %              normal approximation for ion counter shot noise

            %% Constants

            kB = 1.380649e-23; % exact, Joules per Kelvin
            T = 290; % temperature, Kelvin
            R = detector.resistance; % 1e11; % resistance, ohms

            ionsPerCoulomb = 6241509074460762607.776;
            CPSperVolt = ionsPerCoulomb/R;

            % true beams <= this are Poisson distributed
            smallIonBeamCPS = 100;

            %% Johnson noise

            deltaf = 1./integrationTimes; % bandwidth in Hertz = 1/integration time
            JNvarianceInVolts = 4*kB*T*R*deltaf; % volts^2
            JNvarianceInCPS = JNvarianceInVolts * (CPSperVolt)^2;

            %% Shot noise

            % Poisson variance = total ions = (counts per second) / seconds
            % units of cps^2
            PoissonVarianceInCPS = (trueCountRates*detector.gain) ./ integrationTimes;


            %% Create output

            if detector.type == "F"
                % noise is shot noise + Johnson noise
                % output is in cps
                totalVariance = JNvarianceInCPS + PoissonVarianceInCPS;
                ionBeamStdDevInCPS = sqrt(totalVariance);

                ionBeam = random(...
                    'normal', ...
                    trueCountRates * detector.gain, ...
                    ionBeamStdDevInCPS);

            elseif detector.type == "IC"
                % noise is shot noise only
                % output is in cps
                totalVariance = PoissonVarianceInCPS;
                ionBeamStdDevInCPS = sqrt(totalVariance);

                ionBeamSmall = poissrnd(...
                    trueCountRates * detector.gain ./ integrationTimes)...
                    * integrationTimes;

                ionBeamLarge = random(...
                    'normal', ...
                    trueCountRates * detector.gain, ...
                    ionBeamStdDevInCPS);

                ionBeam = ionBeamSmall*(trueCountRates <= smallIonBeamCPS) + ...
                    ionBeamLarge*(trueCountRates > smallIonBeamCPS);

            else %
                error("unrecognized detector type, use F or IC")

            end % if detector.type, for output


        end % function static method simulateIonBeam


    end % protected methods block
end