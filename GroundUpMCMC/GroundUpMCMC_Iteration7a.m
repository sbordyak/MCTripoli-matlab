%% Metropolis Hastings MCMC from the ground up

% 7. Simplified (model 1) mass spectrometer data
% 7a. Add multiple chains for QAQC, rafactor script into functions, and
% perform multiple simulations to test if results consistent with synthetic
% data

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
truth.modelParameterNames = ["$\log(a/b)$"; 
                             "$\log(C_b)$"; 
                             "$re\hspace{-1pt}f_1$"; 
                             "$re\hspace{-1pt}f_2$"];
truth.lograb = log(0.3);   % log(a/b), with (a/b) <= 1 preferred
truth.logCb = log(2e6);  % log(Cb) log(current of isotope b)
truth.ref1 = -1e2; % detector 1, cps
truth.ref2 =  2e2; % detector 2, cps
truth.model = [truth.lograb; truth.logCb; truth.ref1; truth.ref2];

truth.ca = exp(truth.lograb + truth.logCb) + truth.ref1;
truth.cb = exp(truth.logCb) + truth.ref2;

%%
truth.dhat.dvar = updateDataVariance(truth.model, setup);
truth.dhat.int = [truth.ref1*ones(setup.nBLIntegrations,1);
                  truth.ref2*ones(setup.nBLIntegrations,1);
                  truth.ca*ones(setup.nOPIntegrations,1);
                  truth.cb*ones(setup.nOPIntegrations,1)];

% is datum an on-peak measurement? 
truth.dhat.isOP = [false(2*setup.nBLIntegrations,1); 
              true(2*setup.nOPIntegrations,1)];

% detector index for this measurement
truth.dhat.det = [1*ones(setup.nBLIntegrations,1);
            2*ones(setup.nBLIntegrations,1);
            1*ones(setup.nOPIntegrations,1);
            2*ones(setup.nOPIntegrations,1)];

% isotope index for this measurement, 1 = a, 2 = b
truth.dhat.iso = [zeros(2*setup.nBLIntegrations,1);
            1*ones(setup.nOPIntegrations,1);
            2*ones(setup.nOPIntegrations,1)];


truth.G = makeG(truth.model, truth.dhat);
truth.CM = inv(truth.G'*diag(1./truth.dhat.dvar)*truth.G);


%% START SIMULATIONS HERE

tic
setup.nSimulations = 1e4;
for iSim = 1:setup.nSimulations

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

mRough = [rough.lograb; rough.logCb; rough.ref1; rough.ref2];
functionToMinimize = @(m) -loglikLeastSquares(m, data, setup);
opts = optimoptions('fminunc', 'Display', 'off');
modelInitial = fminunc(@(m) functionToMinimize(m), mRough, opts);
%llInitial = -negLogLik;
dvarCurrent = updateDataVariance(modelInitial, setup);
%dhatCurrent = evaluateModel(modelInitial, setup);

% build Jacobian matrix
G = makeG(modelInitial, data);

% least squares model parameter covariance matrix
CM = inv(G'*diag(1./dvarCurrent)*G);



%% initialize model parameters and likelihoods

setup.proposalCov = CM;
setup.nMC = 2e4; % number of MCMC trials
setup.seive = 20;

setup.nmodel = length(modelInitial); % number of model parameters


%% Perturb initial model to ensure convergence

setup.perturbation = 10;
setup.nChains = 8;

% perturb modelCurrent by setup.perturbation standard deviations
modelInitialVector = zeros(setup.nmodel, setup.nChains);
llInitialVector = zeros(1, setup.nChains);
for iChain = 1:setup.nChains
    modelInitialVector(:, iChain) = modelInitial + ...
                   setup.perturbation*randn(4,1).*sqrt(diag(CM));
    dhatCurrent = evaluateModel(modelInitialVector(:,iChain), setup);
    dvarCurrent = updateDataVariance(modelInitialVector(:,iChain), setup);
    llInitialVector(iChain) = loglik(dhatCurrent, data, dvarCurrent);
end


%% Run up multiple chains of Metropolis Hastings
nSavedModels = setup.nMC/setup.seive;
modelChains  = nan([setup.nmodel, nSavedModels, setup.nChains], "double");
loglikChains = nan([1, nSavedModels, setup.nChains], "double");

parfor iChain = 1:setup.nChains

[outputModels, outputLogLiks] = ...
    MetropolisHastings(...
    modelInitialVector(:,iChain), ...
    llInitialVector(iChain), ...
    data, setup);

modelChains(:,:,iChain) = outputModels;
loglikChains(:,:,iChain) = outputLogLiks;

end % parfor iChain = 1:nChains


%% Burn in

setup.burnin = 50;
postBurnInChains = modelChains(:,setup.burnin+1:end,:);


%% plot data and sample fits

% aggregate chains
setup.nPostBurnIn = size(postBurnInChains, 2);
mAll = reshape(postBurnInChains, [setup.nmodel, setup.nPostBurnIn*setup.nChains]);
result(iSim).modelMean = mean(mAll,2);
result(iSim).modelCov = cov(mAll');

% calculate chiSqare with true values
result(iSim).r = result(iSim).modelMean - truth.model;
result(iSim).ChiSq = result(iSim).r' * inv(result(iSim).modelCov) * result(iSim).r;

if ~mod(iSim,10)
    disp("Simulation number " + num2str(iSim))
    toc
end % if ~mod(iSim,10)

end % for iSim = 1:nSimulations
toc


%% Inspect chains for agreement

firstChain = modelChains(:,:,1);
firstChainPostBurnIn = firstChain(:,setup.burnin+1:end,:);
firstChainPostBurnIn(1:2,:) = exp(firstChainPostBurnIn(1:2,:));
figure('Position', [50, 50, 800, 600])
plotmatrix(firstChainPostBurnIn')

makeForestPlot(postBurnInChains)
makeECDFs(postBurnInChains)

%% Inspect results for accuracy

inspectSimulationResults(result, truth)

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%
% Current script above, functions below


%% MCMC - MetropolisHastings

function [outputModels, outputLogLiks] = ...
    MetropolisHastings(modelInitial, llInitial, data, setup)

modelCurrent = modelInitial;
llCurrent = llInitial;

outputModels = nan([setup.nmodel, setup.nMC/setup.seive], "double");
outputLogLiks = nan([1, setup.nMC/setup.seive], "double");

for iMC = 1:setup.nMC

    % save off current model
    if ~mod(iMC,setup.seive)
        outputIndex = iMC/setup.seive;
        outputModels(:, outputIndex) = modelCurrent;
        outputLogLiks(outputIndex) = llCurrent;
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
        llCurrent  = llProposed;

    end % if keep


end % for iMC = 1:nMC

end % function MetropolisHastings()


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


%% calculate log-likelihood as a function of model parameters
% minimize ll for least squares 

function ll = loglikLeastSquares(m, data, setup)

    dhat = evaluateModel(m, setup);
    dvar = updateDataVariance(m, setup);
    ll = loglik(dhat, data, dvar);

end % function loglikLeastSquares


%% assemble design matrix G from model and data
% G is the linearized design matrix in d = Gm

function G = makeG(m, data)

isIsotopeA = data.iso == 1;
isIsotopeB = data.iso == 2;

G = zeros(length(data.int), 4);

% derivative of ca wrt log(a/b): dca__dlogra_b
G(isIsotopeA, 1) = exp(m(1)+m(2));

% derivative of ca wrt log(Cb)
G(isIsotopeA, 2) = exp(m(1)+m(2));

% derivative of cb wrt log(Cb)
G(isIsotopeB, 2) = exp(m(2));

G(data.det == 1, 3) = 1; % derivative wrt ref1
G(data.det == 2, 4) = 1; % derivative wrt ref2

end % function G = makeG(m, data)


%% make a forest plot from the sampled model parameters from multiple chains

function makeForestPlot(modelChains)

vbuffer = 0.15; % vertical buffer at top and bottom
confidenceLevels = [0.68 0.95]; % for histogram

nChains = size(modelChains,3);
figure('Position', [1, 1, 1000, 750], 'Units', 'pixels')
colormap(winter(nChains));

nVariables = size(modelChains, 1);
mRows = floor(sqrt(nVariables));
nColumns = ceil(sqrt(nVariables));
tiledlayout(mRows,nColumns) % some function of nVariables in future

intervalArray1 = zeros(2, nChains); % confidenceLevels(1)
intervalArray2 = zeros(2, nChains); % confidenceLevels(2)
medianArray = zeros(1, nChains); % medians
for iVariable = 1:nVariables

    nexttile
    hold on
    
    for iChain = 1:nChains
    
        sortedChain = sort(modelChains(iVariable, :, iChain));
        interval1 = findShortestInterval(sortedChain, confidenceLevels(1));
        interval2 = findShortestInterval(sortedChain, confidenceLevels(2));

        intervalArray1(:,iChain) = interval1;
        intervalArray2(:,iChain) = interval2;

        nSamples = length(sortedChain);
        medianArray(iChain) = ...
            (...
            sortedChain(ceil(nSamples/2)) + ...
            sortedChain(floor(nSamples/2+1)) ...
            ) / 2; % median without sorting again

    end % for iChain

    yCoords = vbuffer + (1-2*vbuffer)/(nChains-1) * (0:(nChains-1));
    
    line(intervalArray1, [yCoords; yCoords], "LineWidth", 5, "Color", "#0072BD")
    line(intervalArray2, [yCoords; yCoords], "LineWidth", 2, "Color", "#0072BD")
    plot(medianArray, yCoords, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', "#A2142F", 'MarkerEdgeColor', 'k')

    ylim([0, 1])
    set(gca,'YTickLabel',[]);

end % for iVariable



end % function makeForestPlot


%%  Compare empirical cdfs

function makeECDFs(modelChains)

nChains = size(modelChains,3);
figure('Position', [1, 1, 1000, 750], 'Units', 'pixels')
cmap = colormap(winter(nChains));

nVariables = size(modelChains, 1);
mRows = floor(sqrt(nVariables));
nColumns = ceil(sqrt(nVariables));
tiledlayout(mRows,nColumns) % some function of nVariables in future

for iVariable = 1:nVariables

    nexttile
    hold on
    
    for iChain = 1:nChains

        [F, x] = ecdf(modelChains(iVariable, :, iChain));
        plot(x, F, 'LineWidth', 2, 'Color', cmap(iChain,:)) 

    end % function makeECDFs

end % for iVariable

end % function makeECDFs


%% Find shortest empirical confidence interval among samples of distribution

function interval = findShortestInterval(sortedSamples, confidenceLevel)
% find the shortest interval containing a confidenceLevel confidence
% interval
% INPUTS: samples: a sorted vector of samples of the distribution
%         confidenceLevel: 0 < confidenceLevel < 1
% OUTPUTS: interval: two-element vector bounding confidence interval

nSamples = length(sortedSamples);

intervalIndices = ceil(confidenceLevel*nSamples);
upperStop = nSamples - intervalIndices;
shortestInterval = Inf;
lowerLimit = 1;
for iLowerLimit = 1:upperStop
    
    intervalWidth = sortedSamples(iLowerLimit+intervalIndices) - ...
                    sortedSamples(iLowerLimit);

    if intervalWidth < shortestInterval
        shortestInterval = intervalWidth;
        lowerLimit = iLowerLimit;
    end % if intervalWidth < shortestInterval

end % for iSample for limit1

interval = [sortedSamples(lowerLimit), ...
            sortedSamples(lowerLimit+intervalIndices)];

end % function findShortestInterval


%% inspect simulation results

function inspectSimulationResults(result, truth)

fh = figure('Position', [100 100 600 500], 'Units', 'pixels', ...
    'Name', 'Simulation Results', 'NumberTitle','off', ...
    'Toolbar', 'none');
tg = uitabgroup(fh, 'Position', [0 0 1 1], 'Units', 'normalized');

% Chi-square plot
tabChiSq = uitab(tg, "Title", "ChiSquare");
axChiSq = axes(tabChiSq, 'OuterPosition', [0 0 1 1], ...
    'Units', 'normalized');

chiSq = [result(:).ChiSq]';
nModelParams = length(result(1).r);

histogram(axChiSq, chiSq, 'normalization', 'pdf')
hold on
xRange = xlim(gca);
xVec = linspace(xRange(1), xRange(2), 500);
plot(axChiSq, xVec, chi2pdf(xVec, nModelParams), 'LineWidth', 2)
set(axChiSq, 'FontSize', 16)
xlabel(axChiSq, '$x$', 'FontSize', 24, 'Interpreter', 'latex')
ylabel(axChiSq, '$f_4(x)$', 'FontSize', 24, 'Interpreter', 'latex')
title(axChiSq, '$\chi^2$ for Simulations', 'FontSize', 26, 'Interpreter', 'latex')

% Matrix Plot: setup
tabMatrixPlot = uitab(tg, "Title", "MatrixPlot");
axMatrixPlot = axes(tabMatrixPlot, 'OuterPosition', [0 0 1 1], ...
    'Units', 'normalized');

simulationMeans = [result(:).modelMean]';
[S,AX,BigAx,H,HAx] = plotmatrix(axMatrixPlot, simulationMeans);
set(AX, {'NextPlot'}, {'add'})

% Matrix Plot: true values on histograms
for rowModelParam = 1:nModelParams

    histYRange = HAx(rowModelParam).YLim;
    line(HAx(rowModelParam), ...
        truth.model(rowModelParam)*ones(1,2), histYRange, ...
        'LineWidth', 2, 'Color', 'r')

    for colModelParam = 1:nModelParams

        if colModelParam == 1
            ylabel(AX(rowModelParam, colModelParam), ...
                truth.modelParameterNames(rowModelParam), ...
                "FontSize", 14, "Interpreter", "latex");
        end

        if rowModelParam == nModelParams
            xlabel(AX(rowModelParam, colModelParam), ...
                truth.modelParameterNames(colModelParam), ...
                "FontSize", 14, "Interpreter", "latex");
        end

        if rowModelParam == colModelParam, continue, end
        
        trueModel = [truth.model(colModelParam); ...
                     truth.model(rowModelParam)];
        plot(AX(rowModelParam, colModelParam), ...
            trueModel(1), trueModel(2),...
            'o', 'MarkerSize', 6, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')

        subCov = truth.CM([colModelParam rowModelParam],...
                          [colModelParam rowModelParam]);
        thetas = linspace(0, 2*pi, 100); 
        circpoints = 2*[cos(thetas); sin(thetas)];
        cholcov = chol(subCov, "lower"); 
        ellipsepoints = cholcov*circpoints + trueModel;
        plot(AX(rowModelParam, colModelParam), ...
            ellipsepoints(1,:), ellipsepoints(2,:), 'r', 'LineWidth', 2)

    end % for colModelParam

end % for iModelParam

end % function inspectSimulationResults

