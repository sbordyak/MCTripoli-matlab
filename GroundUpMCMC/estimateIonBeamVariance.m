function ionBeamVariance = ...
    estimateIonBeamVariance(trueCountRates, integrationTimes, detector)
%ESTIMATEIONBEAMVARIANCE Estimate ion beam variance
%   Estimate the ion beam variance based on mass spec properties.
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
%   - ionBeamVariance: vector of estimated variances in cps^2

%% Constants

kB = 1.380649e-23; % exact, Joules per Kelvin
T = 290; % temperature, Kelvin
R = detector.resistance; % 1e11; % resistance, ohms

ionsPerCoulomb = 6241509074460762607.776;
CPSperVolt = ionsPerCoulomb/R;

%% Johnson noise

deltaf = 1./integrationTimes; % bandwidth in Hertz = 1/integration time
JNvarianceInVolts = 4*kB*T*R*deltaf; % volts^2
JNvarianceInCPS = JNvarianceInVolts * (CPSperVolt)^2;

%% Shot noise

% Poisson variance = total ions = counts/second * seconds
PoissonVarianceInCPS =(trueCountRates*detector.gain) ./ integrationTimes; % counts^2


%% Create output

if detector.type == "F"
    % noise is shot noise + Johnson noise
    % output is in volts
    ionBeamVariance = JNvarianceInCPS + PoissonVarianceInCPS;

elseif detector.type == "IC"
    % noise is shot noise only
    % output is in cps
    
    ionBeamVariance = PoissonVarianceInCPS;

else % 
    error("unrecognized detector type, use F or IC")
    
end % if detector.type, for output


end % function