%% Metropolis Hastings MCMC from the ground up

% 7. Simplified (model 1) mass spectrometer data

%% synthetic data setup

% parameters for simulated measurement
setup.nBLIntegrations = 1e2;
setup.nOPIntegrations = 1e2;
setup.detector.type = "F";
setup.detector.resistance = 1e11;
setup.detector.gain = 1;
setup.BLIntegrationTimes = ones(setup.nBLIntegrations,1);
setup.OPIntegrationTimes = ones(setup.nOPIntegrations,1);

% true parameters for simulated data
truth.lograb = log(0.03);   % log(a/b)
truth.logCb = log(1e6);  % log(Cb) log(current of isotope b)
truth.ref1 = -1e2; % detector 1, cps
truth.ref2 =  2e2; % detector 2, cps
truth.model = [truth.lograb; truth.logCb; truth.ref1; truth.ref2];

truth.ca = exp(truth.lograb + truth.logCb) + truth.ref1;
truth.cb = exp(truth.logCb) + truth.ref2;

rng(); % start random number stream in one spot

% generate random BL and OP data, assemble into data vector
data.BL_det1 = ...
    simulateIonBeam(truth.ref1, ...
        setup.BLIntegrationTimes, ...
        setup.detector);

data.BL_det2 = ...
    simulateIonBeam(truth.ref2, ...
        setup.BLIntegrationTimes, ...
        setup.detector);

data.OP_det1 = ...
    simulateIonBeam(truth.ca, ...
        setup.OPIntegrationTimes, ...
        setup.detector);

data.OP_det2 = ...
    simulateIonBeam(truth.cb, ...
        setup.OPIntegrationTimes, ...
        setup.detector);

% data vector of measured intensities
data.int = [data.BL_det1; data.BL_det2; data.OP_det1; data.OP_det2];

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


%% Solve maximum likelihood problem to initialize model parameters

isIsotopeA = data.iso == 1;
isIsotopeB = data.iso == 2;
rough.lograb = mean(log(data.int(isIsotopeA)./data.int(isIsotopeB)));
rough.logCb = max(1, mean(real(log(data.int(isIsotopeB)))));

inBL_det1 = ~data.isOP & data.det == 1;
inBL_det2 = ~data.isOP & data.det == 2;
rough.ref1 = mean(data.int(inBL_det1));
rough.ref2 = mean(data.int(inBL_det2));

m0 = [rough.lograb; rough.logCb; rough.ref1; rough.ref2];
[modelCurrent, negLogLik] = fminunc(@(m) -loglikLeastSquares(m, data, setup), m0);
llCurrent = -negLogLik;
dvarCurrent = updateDataVariance(modelCurrent, setup);
dhatCurrent = evaluateModel(modelCurrent, setup);

% build Jacobian matrix
% derivatives of [BL_det1, BL_det2, OP_det1, OP_det2]
% with respect to [log(a/b), log(Cb), ref1, ref2]
G = zeros(length(data.int), 4);

% derivative of ca wrt log(a/b)
G(isIsotopeA, 1) = exp(modelCurrent(1)+modelCurrent(2));

% derivative of ca wrt log(Cb)
G(isIsotopeA, 2) = exp(modelCurrent(1)+modelCurrent(2));

% derivative of cb wrt log(Cb)
G(isIsotopeB, 2) = exp(modelCurrent(2));

G(data.det == 1, 3) = 1; % derivative wrt ref1
G(data.det == 2, 4) = 1; % derivative wrt ref2

% least squares model parameter covariance matrix
CM = inv(G'*diag(1./dvarCurrent)*G);


%% initialize model parameters and likelihoods


setup.proposalCov = CM;
setup.nMC = 1e6; % number of MCMC trials
setup.seive = 20;

setup.nmodel = length(modelCurrent); % number of model parameters


%% Test 1: perturb modelCurrent to ensure convergence
% perturb modelCurrent by setup.perturbation standard deviations

setup.perturbation = 0;
modelCurrent = modelCurrent + setup.perturbation*randn(4,1).*sqrt(diag(CM));
dhatCurrent = evaluateModel(modelCurrent, setup);
dvarCurrent = updateDataVariance(modelCurrent, setup);
llCurrent = loglik(dhatCurrent, data, dvarCurrent);


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

% convert logratio and log-intensity to ratio and intensity
outputSamples = modelSamples;
outputSamples(1:2,:) = exp(outputSamples(1:2,:)); % for ratio and intensity

results.mean = mean(outputSamples, 2);
% results.meanSlope      = results.mean(1);
% results.meanYIntercept = results.mean(2);
% correlation plot
figure('Position', [50, 50, 800, 600])
plotmatrix(outputSamples')

% figure('Position', [900 50, 700, 600])
% scatter(modelSamples(1,:), modelSamples(2,:), 'Marker', 'o', ...
%     'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2, ...
%     'MarkerEdgeColor','none', 'SizeData', 10)
% hold on


%% theoretical OLS results

% FROM ITERATION 6:
% % repeated from initialization routine
% X = [setup.x, ones(setup.ndata,1)];
% beta = X \ data;
% betaCov = var(data - [setup.x, ones(setup.ndata,1)]*beta) * inv(X'*X);
% 
% % ellipse from OLS fit
% thetas = linspace(0, 2*pi, 100); circpoints = 2*[cos(thetas); sin(thetas)];
% cholcov = chol(betaCov, "lower"); 
% ellipsepoints = cholcov*circpoints + beta;
% plot(ellipsepoints(1,:), ellipsepoints(2,:), 'r', 'LineWidth', 2)
% plot(beta(1), beta(2), '.r', 'MarkerSize', 25)

%% plot data and sample fits
% figure('Position', [20, 20, 1200, 1000])
% xmin = -1;
% xmax = 11;
% 
% nlines = 1e3;
% 
% linesymin = modelSamples(1,1:nlines)*xmin + modelSamples(2,1:nlines);
% linesymax = modelSamples(1,1:nlines)*xmax + modelSamples(2,1:nlines);
% 
% line([xmin; xmax], [linesymin; linesymax], 'Color', [0.8500 0.3250 0.0980, 0.05])
% hold on
% plot(setup.x, data, '.', 'MarkerSize', 20, 'MarkerEdgeColor', 'k')
% hax = gca;
% set(hax, 'FontSize', 24)


%% evaluate this particular model

function dhat = evaluateModel(m, setup)

    % y = slope * x       + y-intercept
    % dhat = m(1) * setup.x + m(2);

    lograb = m(1);
    logCb  = m(2) * ones(setup.nOPIntegrations, 1);
    ref1   = m(3);
    ref2   = m(4);

    BL_det1 = ref1*ones(setup.nBLIntegrations,1);
    BL_det2 = ref2*ones(setup.nBLIntegrations,1);
    OP_det1 = exp(lograb + logCb) + ref1;
    OP_det2 = exp(logCb) + ref2;

    dhat = [BL_det1; BL_det2; OP_det1; OP_det2];

end % function evaluateModel


%% update data covariance matrix

function dvar = updateDataVariance(m, setup)

    % estimate ion beam variances from model parameters
    lograb = m(1);
    logCb = m(2);

    BL_det1_var = estimateIonBeamVariance(...
        zeros(setup.nBLIntegrations,1), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    BL_det2_var = estimateIonBeamVariance(...
        zeros(setup.nBLIntegrations,1), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    OP_det1_var = estimateIonBeamVariance(...
        exp(lograb + logCb), ...
        setup.BLIntegrationTimes, ...
        setup.detector);
    OP_det2_var = estimateIonBeamVariance(...
        exp(logCb), ...
        setup.BLIntegrationTimes, ...
        setup.detector);

    dvar = [BL_det1_var; BL_det2_var; OP_det1_var; OP_det2_var];

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