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
setup.BLTimes = cumsum(setup.BLIntegrationTimes);
setup.OPTimes = max(setup.BLTimes) + 5 + cumsum(setup.OPIntegrationTimes);


% true parameters for simulated data
truth.modelParameterNames = ["$\log(a/b)$"; 
                             "$\log(C_b)$"; 
                             "$re\hspace{-1pt}f_1$"; 
                             "$re\hspace{-1pt}f_2$"];
truth.speciesNames = ["a"; "b"];
truth.lograb = log(0.3);   % log(a/b), with (a/b) <= 1 preferred
truth.logCb = log(2e6);  % log(Cb) log(current of isotope b)
truth.ref1 = -1e2; % detector 1, cps
truth.ref2 =  2e2; % detector 2, cps
truth.model = [truth.lograb; truth.logCb; truth.ref1; truth.ref2];


truth.ca = exp(truth.lograb + truth.logCb) + truth.ref1;
truth.cb = exp(truth.logCb) + truth.ref2;

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
setup.nSimulations = 1e1;
results = []; % hard to initialize a struct.
for iSim = 1:setup.nSimulations

rng(); % start random number stream in one spot

data = syntheticData(truth, setup);

%% Solve maximum likelihood problem to initialize model parameters

maxlik = maxLikelihood(data, setup);

%% initialize model parameters and likelihoods

setup.proposalCov = maxlik.CM;
setup.nMC = 2e4; % number of MCMC trials
setup.seive = 20;
setup.nmodel = length(truth.model); % number of model parameters
setup.nChains = 8;
setup.perturbation = 10;


%% Perturb initial model to ensure convergence

[initModels, initLogLiks] = initializeChains(setup, data, maxlik);


%% Run up multiple chains of Metropolis Hastings

nSavedModels = setup.nMC/setup.seive;
modelChains  = nan([setup.nmodel, nSavedModels, setup.nChains], "double");
loglikChains = nan([1, nSavedModels, setup.nChains], "double");

parfor iChain = 1:setup.nChains

[outputModels, outputLogLiks] = ...
    MetropolisHastings(...
    initModels(:,iChain), ...
    initLogLiks(iChain), ...
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

% update command window
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


%% create synthetic dataset
% generate random BL and OP data, assemble into data vector

function data = syntheticData(truth, setup)

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

data.BLTimes = cumsum(setup.BLIntegrationTimes);
data.OPTimes = max(setup.BLTimes) + 5 + cumsum(setup.OPIntegrationTimes);

end % function syntheticData

%% Maximum likelihood estimate from data

function maxlik = maxLikelihood(data, setup)
% maxlik is a struct with three fields:
%   model is the maximum likelihood model
%   dvar is the expected variance in the data given the model parameters
%   CM is the covariance matrix of the max likelihood model parameters

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
maxlik.model = fminunc(@(m) functionToMinimize(m), mRough, opts);
%llInitial = -negLogLik;
maxlik.dvar = updateDataVariance(maxlik.model, setup);
%dhatCurrent = evaluateModel(modelInitial, setup);

% build Jacobian matrix
G = makeG(maxlik.model, data);

% least squares model parameter covariance matrix
maxlik.CM = inv(G'*diag(1./maxlik.dvar)*G);

end % function maxLikelihood


%% Initialize chains with perturbed max likelihood models
% start MCMC chains at perturbed positions

function [initModels, initLogLiks] = initializeChains(setup, data, maxlik)

% perturb modelCurrent by setup.perturbation standard deviations
initModels = zeros(setup.nmodel, setup.nChains);
initLogLiks = zeros(1, setup.nChains);
for iChain = 1:setup.nChains
    initModels(:, iChain) = maxlik.model + ...
                   setup.perturbation*randn(setup.nmodel,1) .* ...
                   sqrt(diag(maxlik.CM));
    dhatCurrent = evaluateModel(initModels(:,iChain), setup);
    dvarCurrent = updateDataVariance(initModels(:,iChain), setup);
    initLogLiks(iChain) = loglik(dhatCurrent, data, dvarCurrent);
end

end % initializeChains()



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
[~,AX,BigAx,~,HAx] = plotmatrix(axMatrixPlot, simulationMeans);
set(AX, {'NextPlot'}, {'add'})
BigAx.Title.String = "Simulation Results vs. True Values";
BigAx.Title.FontSize = 20;
BigAx.Title.FontWeight = 'normal';

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


%% Inspect models' fit to data

function inspectModelFitToData(data, truth, result, setup)
% INPUTS: data and truth structs from synthetic data generation
% result must be struct with dimension 1, e.g. result(1)


fh = figure('Position', [5 5 1000 700], 'Units', 'pixels', ...
    'Name', 'Data Fits', 'NumberTitle','off');
tg = uitabgroup(fh, 'Position', [0 0 1 1], 'Units', 'normalized');

% Baselines: Data wrangling
nDet = max(data.det);
nBLIntegrations = sum(data.det == 1 & ~data.isOP);
BLmodelIndices = [3 4];

% Baselines: Axes setup
tabBL = uitab(tg, "Title", "Baselines");
axBLAll = axes(tabBL, 'OuterPosition', [0 0 1 1], ...
    'Units', 'normalized');

for iDet = 1:nDet

    axesiBL = subplot(nDet, 1, iDet);
    axesiBL.NextPlot = 'add';
    axesiBL.FontSize = 16;
    modelIndex = BLmodelIndices(iDet);
    BLname = "BL" + num2str(iDet);

    inBL_iDet = ~data.isOP & data.det == iDet;

    % Baselines: Uncertainty in model
    iBL_2sigma = 2*sqrt(result.modelCov(modelIndex,modelIndex)) * ...
        ones(nBLIntegrations,1);
    iBL_Unct_X = [data.BLTimes; data.BLTimes(end:-1:1)];
    iBL_Unct_Y = result.modelMean(modelIndex) + [-iBL_2sigma;  iBL_2sigma];
    patch(axesiBL, 'XData', iBL_Unct_X, 'YData', iBL_Unct_Y, ...
        'FaceColor', "#77AC30", 'EdgeColor', 'none', 'FaceAlpha', 0.3)
    %Baselines: model
    plot(axesiBL, data.BLTimes, ...
        result.modelMean(modelIndex) * ones(nBLIntegrations,1), ...
        '-', 'Color', 'r', 'LineWidth', 2)
    %Baselines: data
    plot(axesiBL, data.BLTimes, data.int(inBL_iDet), '.', ... 
        'MarkerSize', 12, 'Color', lines(1))
    %Baselines: axes labels and title
    xlabel('Time (seconds)', 'FontSize', 20)
    ylabel('Intensity (cps)', 'FontSize', 20)
    text(axesiBL, 0.04, 0.9, BLname, "Units", "normalized", "FontSize", 24);

end % for iBL


% On Peak Intensity Measurements ca and cb: Data wrangling
nIso = max(data.iso);
nSeq = 1; % for now
nOPIntegrations = sum(data.det == 1 & data.isOP);
logRatioModelIndices = 1;
logIntensityModelIndices = 2;

dhat = evaluateModel(result.modelMean, setup);

% On Peak: Axes setup
tabBL = uitab(tg, "Title", "On Peak");
axOPAll = axes(tabBL, 'OuterPosition', [0 0 1 1], ...
    'Units', 'normalized');

for iIso = 1:nIso

    dataIndices = data.iso == iIso;
    axesiIso = subplot(nIso, 1, iIso);
    axesiIso.NextPlot = 'add';
    axesiIso.FontSize = 16;

    isoData = data.int(dataIndices);
    isoDhat = dhat(dataIndices);
    intensityName = "c_{" + truth.speciesNames(iIso) + "}";
    
    %On Peak: dhat uncertainties

    G = makeG(result.modelMean, data);
    CDhat = G*result.modelCov*G';
    iIso_2sigma = 2*sqrt(diag(CDhat(dataIndices,dataIndices)));
    iIso_Unct_X = [data.OPTimes; data.OPTimes(end:-1:1)];
    iIso_Unct_Y = [isoDhat; isoDhat] + [-iIso_2sigma;  iIso_2sigma];
    patch(axesiIso, 'XData', iIso_Unct_X, 'YData', iIso_Unct_Y, ...
        'FaceColor', "#77AC30", 'EdgeColor', 'none', 'FaceAlpha', 0.3)

    %On Peak: dhat
    plot(axesiIso, data.OPTimes, ...
        isoDhat, '-', 'Color', 'r', 'LineWidth', 2)
    %On Peak: data
    plot(axesiIso, data.OPTimes, isoData, '.', 'MarkerSize', 12, ...
        'Color', lines(1))
    %On Peak: axes labels and title
    xlabel('Time (seconds)', 'FontSize', 20)
    ylabel('Intensity (cps)', 'FontSize', 20)
    text(axesiIso, 0.04, 0.9, intensityName, ...
        "Units", "normalized", "FontSize", 24, 'Interpreter', 'tex');

end % for iIso

end % function inspectDataFit


%%

inspectModelFitToData(data, truth, result(2), setup)
