%% Metropolis Hastings MCMC from the ground up

% 7. Simplified (model 1) mass spectrometer data

%% synthetic data setup

% parameters for simulated measurement
setup.nBLIntegrations = 100;
setup.nOPIntegrations = 100;
setup.detector.type = "F";
setup.detector.resistance = 1e11;
setup.detector.gain = 1;
setup.BLIntegrationTimes = ones(setup.nBLIntegrations,1);
setup.OPIntegrationTimes = ones(setup.nOPIntegrations,1);

% true parameters for simulated data
truth.lograb = log(2);   % log(a/b)
truth.logCb = log(1e6);  % log(Cb) log(current of isotope b)
truth.ref1 = -1e2; % detector 1, cps
truth.ref2 =  2e2; % detector 2, cps
truth.model = [truth.lograb; truth.logCb; truth.ref1; truth.ref2];

truth.ca = exp(truth.lograb + truth.logCb) + truth.ref1;
truth.cb = exp(truth.logCb) + truth.ref2;

rng(2) % start random number stream in one spot

% generate random BL and OP data, assemble into data vector
data.BL1 = ...
    simulateIonBeam(truth.ref1, ...
        setup.BLIntegrationTimes, ...
        setup.detector);

data.BL2 = ...
    simulateIonBeam(truth.ref2, ...
        setup.BLIntegrationTimes, ...
        setup.detector);

data.OP1 = ...
    simulateIonBeam(truth.ca, ...
        setup.OPIntegrationTimes, ...
        setup.detector);

data.OP2 = ...
    simulateIonBeam(truth.cb, ...
        setup.OPIntegrationTimes, ...
        setup.detector);

% data vector of measured intensities
data.int = [data.BL1; data.BL2; data.OP1; data.OP2];

% data flags for data.int

% is datum an on-peak measurement? 
data.isOP = [false(2*setup.nBLIntegrations,1); 
              true(2*setup.nOPIntegrations,1)];

% detector index for this measurement
data.det = [1*ones(setup.nBLIntegrations,1);
            2*ones(setup.nBLIntegrations,1);
            1*ones(setup.nOPIntegrations,1);
            2*ones(setup.nOPIntegrations,1)];

% isotope index for this measurement, 1 = a, 2 = b
data.iso = [zeros(2*setup.nBLIntegrations,1);
            1*ones(setup.nOPIntegrations,1);
            2*ones(setup.nOPIntegrations,1)];


%% Solve least squares problem to initialize model parameters

isIsotopeA = data.iso == 1;
isIsotopeB = data.iso == 2;
rough.lograb = mean(log(data.int(isIsotopeA)./data.int(isIsotopeB)));
rough.logCb = max(1, mean(real(log(data.int(isIsotopeB)))));

inBL1 = ~data.isOP & data.det == 1;
inBL2 = ~data.isOP & data.det == 2;
rough.ref1 = mean(data.int(inBL1));
rough.ref2 = mean(data.int(inBL2));

m0 = [rough.lograb; rough.logCb; rough.ref1; rough.ref2];
[modelCurrent, llCurrent] = fminunc(@(m) -loglikLeastSquares(m, data, setup), m0);



%% initialize model parameters and likelihoods

% some maximum likelihood calculations to get close
X = [setup.x, ones(setup.ndata,1)];
beta = X \ data; % maximum likelihood solution of model parameters
noiseMean = var(data - [setup.x, ones(setup.ndata,1)]*beta);
betaCov = noiseMean * inv(X'*X);
noiseVariance = 2*noiseMean^2/(setup.ndata - 1);

modelCurrent = [beta; noiseMean];
noiseCurrent = noiseMean;

% initialize data variance vector, dhat vector, and log-likelihood
dvarCurrent = noiseCurrent * ones(setup.ndata, 1);
dhatCurrent = evaluateModel(modelCurrent, setup);
llCurrent = loglik(dhatCurrent, data, dvarCurrent);

setup.proposalCov = blkdiag(betaCov, noiseVariance);
setup.nMC = 1e7; % number of MCMC trials
setup.seive = 20;

setup.nmodel = length(modelCurrent); % number of model parameters


%% MCMC

modelSamples = nan([setup.nmodel, setup.nMC/setup.seive], "double");
for iMC = 1:setup.nMC

    % save off current model
    if ~mod(iMC,setup.seive)
        modelSamples(:, iMC/setup.seive) = modelCurrent;
    end
    
    % propose a new model
    modelProposed = modelCurrent + ...
                    mvnrnd( zeros(setup.nmodel,1), ...
                            setup.proposalCov )';

    % calculate residuals for old and new models
    dhatProposed = evaluateModel(modelProposed, setup);
    
    % create data covariance with current and proposed noise terms
    dvarProposed = updateDataVariance(modelProposed, setup);
    
    % calculate log-likelihoods of current and proposed samples
    llProposed = loglik(dhatProposed, data, dvarProposed);

    % difference in log-likelihood = log(likelihood ratio)
    delta_ll = llProposed - llCurrent;
    
    % probability of keeping the proposed model
    keep = min(1, exp(delta_ll)); 

    % keep the proposed model with probability = keep
    if keep > rand(1) % if we're keeping the proposed model

        % the proposed model becomes the current one
        modelCurrent = modelProposed; 
        dhatCurrent = dhatProposed;
        dvarCurrent = dvarProposed;
        llCurrent  = llProposed;

    end % if keep


end % for iMC = 1:nMC

%% MCMC results

results.mean = mean(modelSamples, 2);
results.meanSlope      = results.mean(1);
results.meanYIntercept = results.mean(2);
% correlation plot
figure('Position', [50, 50, 800, 600])
plotmatrix(modelSamples')

figure('Position', [900 50, 700, 600])
scatter(modelSamples(1,:), modelSamples(2,:), 'Marker', 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2, ...
    'MarkerEdgeColor','none', 'SizeData', 10)
hold on


%% theoretical OLS results

% repeated from initialization routine
X = [setup.x, ones(setup.ndata,1)];
beta = X \ data;
betaCov = var(data - [setup.x, ones(setup.ndata,1)]*beta) * inv(X'*X);

% ellipse from OLS fit
thetas = linspace(0, 2*pi, 100); circpoints = 2*[cos(thetas); sin(thetas)];
cholcov = chol(betaCov, "lower"); 
ellipsepoints = cholcov*circpoints + beta;
plot(ellipsepoints(1,:), ellipsepoints(2,:), 'r', 'LineWidth', 2)
plot(beta(1), beta(2), '.r', 'MarkerSize', 25)

%% plot data and sample fits
figure('Position', [20, 20, 1200, 1000])
xmin = -1;
xmax = 11;

nlines = 1e3;

linesymin = modelSamples(1,1:nlines)*xmin + modelSamples(2,1:nlines);
linesymax = modelSamples(1,1:nlines)*xmax + modelSamples(2,1:nlines);

line([xmin; xmax], [linesymin; linesymax], 'Color', [0.8500 0.3250 0.0980, 0.05])
hold on
plot(setup.x, data, '.', 'MarkerSize', 20, 'MarkerEdgeColor', 'k')
hax = gca;
set(hax, 'FontSize', 24)


%% evaluate this particular model

function dhat = evaluateModel(m, setup)

    % y = slope * x       + y-intercept
    % dhat = m(1) * setup.x + m(2);

    lograb = m(1);
    logCb  = m(2) * ones(setup.nOPIntegrations, 1);
    ref1   = m(3);
    ref2   = m(4);

    BL1 = ref1*ones(setup.nBLIntegrations,1);
    BL2 = ref2*ones(setup.nBLIntegrations,1);
    OP1 = exp(lograb + logCb) + ref1;
    OP2 = exp(logCb) + ref2;

    dhat = [BL1; BL2; OP1; OP2];

end % function evaluateModel


%% update data covariance matrix

function dvar = updateDataVariance(m, setup)

    % estimate ion beam variances from model parameters
    lograb = m(1);
    logCb = m(2);

    BL1var = estimateIonBeamVariance(...
        zeros(setup.nBLIntegrations,1), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    BL2var = estimateIonBeamVariance(...
        zeros(setup.nBLIntegrations,1), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    OP1var = estimateIonBeamVariance(...
        exp(lograb + logCb), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    OP2var = estimateIonBeamVariance(...
        exp(logCb), ...
        setup.BLIntegrationTimes, ...
        setup.detector);

    dvar = [BL1var; BL2var; OP1var; OP2var];

end % updateDataCovariance


%% calculate log-likelihood (including -1/2 term)
% ignore constant terms (for now)

function ll = loglik(dhat, data, dvar)
    
    % enforce bounds on model prior
    % if m(3) <= 0
    %     ll = -Inf;
    %     return
    % end

    residuals = (data.int - dhat);
    chiSqTerms = residuals.^2 ./ dvar;
    ll = -1/2 * sum(chiSqTerms) - 1/2 * sum(log(dvar));

end % function loglik

%% caluclate log-likelihood as a function of model parameters
% minimize ll for least squares 

function ll = loglikLeastSquares(m, data, setup)

    dhat = evaluateModel(m, setup);
    dvar = updateDataVariance(m, setup);
    ll = loglik(dhat, data, dvar);

end % function loglikLeastSquares