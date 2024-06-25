%% make some mass spectrometer data, adding complexity along the way

% start with two types of noise: poisson (counts) and normal (Johnson)
% input: ion beam size in counts per second
% output: voltages on Faradays

%% measurement parameters (property of dataset/method)

integrationTime = 1; % second

intensityInCPS = 1e6; % true ion beam intensity, cps


%% Johnson-Nyquist noise

kB = 1.380649e-23; % exact, Joules per Kelvin
T = 290; % temperature, Kelvin
R = 1e11; % resistance, ohms

deltaf = 1/integrationTime; % bandwidth in Hertz = 1/integration time

JNvarianceInVolts = 4*kB*T*R*deltaf; % volts^2


%% Poisson (shot) noise

ionsPerCoulomb = 6241509074460762607.776;
voltsPerCPS = R/ionsPerCoulomb;
intensityInVolts = intensityInCPS * voltsPerCPS;

% Poisson variance = total ions = counts/second * seconds
PoissonVarianceInCPS = intensityInCPS * integrationTime; % counts^2
PoissonVarianceInVolts = voltsPerCPS^2 * PoissonVarianceInCPS;

totalVariance = JNvarianceInVolts + PoissonVarianceInVolts;
ionBeamStdDevInVolts = sqrt(totalVariance);


syntheticDataInVolts = random(...
    'normal', ...
    intensityInVolts, ...
    ionBeamStdDevInVolts, ...
    [100, 1]);

plot(syntheticDataInVolts, '.', 'MarkerSize',25)
xlabel('Time (seconds)')
ylabel('Measured Intensity (Volts)')
set(gca, 'FontSize', 18)