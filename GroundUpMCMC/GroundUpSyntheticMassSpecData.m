%% make some mass spectrometer data, adding complexity along the way

% start with two types of noise: poisson (counts) and normal (Johnson)
% input: ion beam size in counts per second
% output: voltages on Faradays

%% measurement parameters (property of dataset/method)

% integrationTime = 1; % second
% 
% intensityInCPS = 1e6; % true ion beam intensity, cps
% 
% 
% %% Johnson-Nyquist noise
% 
% kB = 1.380649e-23; % exact, Joules per Kelvin
% T = 290; % temperature, Kelvin
% R = 1e11; % resistance, ohms
% 
% deltaf = 1/integrationTime; % bandwidth in Hertz = 1/integration time
% 
% JNvarianceInVolts = 4*kB*T*R*deltaf; % volts^2
% 
% 
% %% Poisson (shot) noise
% 
% ionsPerCoulomb = 6241509074460762607.776;
% voltsPerCPS = R/ionsPerCoulomb;
% intensityInVolts = intensityInCPS * voltsPerCPS;
% countsPerIntegration = intensityInCPS * integrationTime;
% 
% % Poisson variance = total ions = counts/second * seconds
% PoissonVarianceInCPS = intensityInCPS * integrationTime; % counts^2
% PoissonVarianceInVolts = voltsPerCPS^2 * PoissonVarianceInCPS;
% 
% totalVariance = JNvarianceInVolts + PoissonVarianceInVolts;
% ionBeamStdDevInVolts = sqrt(totalVariance);
% 
% nData = 100;
% syntheticDataInVolts = ...
%     random('poisson', countsPerIntegration, [100, 1])*voltsPerCPS + ...
%     random('normal', 0, ionBeamStdDevInVolts, [100, 1]);
% 
% plot(syntheticDataInVolts, '.', 'MarkerSize',25)
% xlabel('Time (seconds)')
% ylabel('Measured Intensity (Volts)')
% set(gca, 'FontSize', 18)


%% refactor a bit

%% SETUP

% constants
const.kB = 1.380649e-23; % exact, Joules per Kelvin
const.T = 290; % temperature, Kelvin
const.ionsPerCoulomb = 6241509074460762607.776;

% user inputs: mass spec and method
massSpec.R = [1e12, 1e12];
massSpec.detectorTypes = ["F", "F"]; % "F" or "IC"

method.integrationTime = 0.1; % seconds of ion beam integration per data point
method.totalTime = 100; % seconds
method.outputFormat = ["Volts", "Volts"]; % "Volts" or "CPS"

% user inputs: ion beam behavior and sample properties (model parameters) 
ionBeam.monitorBeamIntensityFunctionCPS = @(t) 0e6 * ones(size(t)); % cps
ionBeam.relativeAbundances = [2, 1]; % monitor beam has abundance = 1


%% Calculations

% expand method
method.nIntegrations = method.totalTime/method.integrationTime;
method.times = linspace(method.integrationTime, ...
    method.totalTime, method.nIntegrations)';
method.integrationTimeVector = method.integrationTime * ...
    ones(method.nIntegrations, 1);
method.outputFormat = categorical(method.outputFormat, ["Volts", "CPS"]); 

% calculate Johnson noise for each integration, based on 
% integration times and amplifier resistances
massSpec.JNvarianceInVolts = ... % volts^2
    4*const.kB * const.T * (1./method.integrationTimeVector) * massSpec.R; 
massSpec.voltsPerCPS = massSpec.R/const.ionsPerCoulomb;
massSpec.detectorTypes = categorical(massSpec.detectorTypes, ["F", "IC"]);

% calculate ion beam intensity and total counts for each integration
% This is the model.  Currently: 
ionBeam.intensityInCPS = ionBeam.monitorBeamIntensityFunctionCPS(method.times);
ionBeam.countsPerIntegration = ionBeam.intensityInCPS .* ...
    method.integrationTimeVector .* ionBeam.relativeAbundances;

% create synthetic data as volts, cps, and output according to outputFormat
syndata.volts = ...
    random('poisson', ionBeam.countsPerIntegration) .* ...
    massSpec.voltsPerCPS + ...
    (massSpec.detectorTypes == "F") .* ...
    random('normal', 0, sqrt(massSpec.JNvarianceInVolts));

syndata.cps = syndata.volts ./ massSpec.voltsPerCPS;

syndata.output = ...
    (method.outputFormat == "Volts") .* syndata.volts + ...
    (method.outputFormat == "CPS") .* syndata.cps;

