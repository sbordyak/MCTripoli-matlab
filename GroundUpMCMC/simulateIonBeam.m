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


end % function