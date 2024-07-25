function ionBeam = simulateIonBeam(countRates, integrationTimes, detector)
%SIMULATEIONBEAM Create synthetic measured ion beam 
%   This simple, pared-down version returns a simulated ion beam.
%   Noise is the sum shot noise and, if Farday output is selected, 
%   Johnson-Nyquist noise.
%
%   Inputs:
%   - countRates: vector of count rates in cps
%   - integrationTimes: vector of integration times in cps
%   - detector: struct containing fields
%       - type = "F" or "IC" for "Faraday" or "Ion Counter"
%       - resistance = eg 1e11 or 1e12
%       - gain = eg 1 or 0.9
%
%   Output:
%   - ionBeam: vector of measured intensities, with
%              units of volts if type = "F", cps if type = "IC"

%% Constants

kB = 1.380649e-23; % exact, Joules per Kelvin
T = 290; % temperature, Kelvin
R = detector.resistance; % 1e11; % resistance, ohms

ionsPerCoulomb = 6241509074460762607.776;


%% Johnson noise

deltaf = 1./integrationTimes; % bandwidth in Hertz = 1/integration time
JNvarianceInVolts = 4*kB*T*R*deltaf; % volts^2


%% Shot noise

voltsPerCPS = R/ionsPerCoulomb;
intensityInVolts = countRates * voltsPerCPS;

% Poisson variance = total ions = counts/second * seconds
PoissonVarianceInCPS = countRates .* integrationTimes; % counts^2
PoissonVarianceInVolts = voltsPerCPS^2 * PoissonVarianceInCPS;


%% Create output

if detector.type == "F"
    % noise is shot noise + Johnson noise
    % output is in volts
    totalVariance = JNvarianceInVolts + PoissonVarianceInVolts;
    ionBeamStdDevInVolts = sqrt(totalVariance);

    ionBeam = random(...
        'normal', ...
        intensityInVolts * detector.gain, ...
        ionBeamStdDevInVolts);

elseif detector.type == "IC"
    % noise is shot noise only
    % output is in cps
    
    ionBeamCounts = random(...
        'Poisson', ...
        PoissonVarianceInCPS * detector.gain);
    ionBeam = ionBeamCounts ./ integrationTimes;

else % 
    error("unrecognized detector type, use F or IC")
    
end % if detector.type, for output


end % function