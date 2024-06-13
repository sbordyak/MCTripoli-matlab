%% Metropolis Hastings MCMC from the ground up

% 6. Performance improvements, add seive and data plot

%% synthetic data setup

setup.ndata = 500;

truth.slope = 2;
truth.yIntercept = 4;
truth.yUncertainty = 10; % standard deviation of normally distributed errors
truth.model = [truth.slope; truth.yIntercept]; % slope, y-intercept

rng(1) % start random number stream in one spot

setup.x = random('uniform', 0, 10, [setup.ndata, 1]);
truth.y = truth.slope * setup.x + truth.yIntercept;

data = random('normal', 0, truth.yUncertainty, [setup.ndata, 1]) + truth.y;


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
llCurrent = loglik(dhatCurrent, data, dvarCurrent, modelCurrent);

setup.proposalCov = blkdiag(betaCov, noiseVariance);
setup.nMC = 1e6; % number of MCMC trials
setup.seive = 10;

setup.nmodel = length(modelCurrent); % number of model parameters


%%

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
    llProposed = loglik(dhatProposed, data, dvarProposed, modelProposed);

    % difference in log-likelihood = log(likelihood ratio)
    delta_ll = llProposed - llCurrent;
    
    % probability of keeping the proposed model
    keep = min(1, exp(delta_ll)); % negative one-half from likelihood eqn

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
    dhat = m(1) * setup.x + m(2);

end % function evaluateModel


%% update data covariance matrix

function dvar = updateDataVariance(m, setup)

    dvar = m(3) * ones(setup.ndata, 1);

end % updateDataCovariance


%% calculate log-likelihood (including -1/2 term)
% ignore constant terms (for now)

function ll = loglik(dhat, data, dvar, m)
    
    % enforce bounds on model prior
    if m(3) <= 0
        ll = -Inf;
        return
    end

    residuals = (data - dhat);
    chiSqTerms = residuals.^2 ./ dvar;
    ll = -1/2 * sum(chiSqTerms) - 1/2 * sum(log(dvar));

end % function loglik