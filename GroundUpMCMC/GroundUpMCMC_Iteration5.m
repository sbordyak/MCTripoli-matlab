%% Metropolis Hastings MCMC from the ground up

% 5. Add y uncertainty as model parameter to linear regression
% try to add as an extra model parameter, not separately proposed
% to avoid numerical instability in log(det()), use variances only

%% synthetic data setup

setup.ndata = 20;

truth.slope = 2;
truth.yIntercept = 4;
truth.yUncertainty = 1; % standard deviation of normally distributed errors
truth.model = [truth.slope; truth.yIntercept]; % slope, y-intercept

rng(1) % start random number stream in one spot

setup.x = random('uniform', 0, 10, [setup.ndata, 1]);
truth.y = truth.slope * setup.x + truth.yIntercept;

data = random('normal', 0, truth.yUncertainty, [setup.ndata, 1]) + truth.y;


%% initialize model parameters

modelCurrent = [setup.x, ones(setup.ndata,1)] \ data; 
noiseCurrent = var(data - [setup.x, ones(setup.ndata,1)]*modelCurrent);
modelCurrent = [modelCurrent; noiseCurrent]; % append noise to model
dvarCurrent = noiseCurrent * ones(setup.ndata, 1);

setup.proposalCov = [0.015 0 0; 0 0.3 0; 0 0 0.5]; % proposal distribution
setup.nMC = 1e6; % number of MCMC trials

setup.nmodel = length(modelCurrent); % number of model parameters


%%

modelSamples = nan([setup.nmodel, setup.nMC], "double");
for iMC = 1:setup.nMC

    % save off current model
    modelSamples(:, iMC) = modelCurrent;
    
    % propose a new model
    modelProposed = modelCurrent + ...
                    mvnrnd( zeros(setup.nmodel,1), ...
                            setup.proposalCov )';

    % calculate residuals for old and new models
    dhatCurrent  = evaluateModel(modelCurrent,  setup);
    dhatProposed = evaluateModel(modelProposed, setup);
    
    % create data covariance with current and proposed noise terms
    dvarCurrent  = updateDataVariance(modelCurrent,  setup);
    dvarProposed = updateDataVariance(modelProposed, setup);
    
    % calculate log-likelihoods of current and proposed samples
    llCurrent  = loglik(dhatCurrent,  data, dvarCurrent,  modelCurrent);
    llProposed = loglik(dhatProposed, data, dvarProposed, modelProposed);

    % difference in log-likelihood = log(likelihood ratio)
    delta_ll = llProposed - llCurrent;
    
    % probability of keeping the proposed model
    keep = min(1, exp(delta_ll)); % negative one-half from likelihood eqn

    % keep the proposed model with probability = keep
    if keep > rand(1) % if we're keeping the proposed model

        modelCurrent = modelProposed; % the proposed model becomes the current one

    end % if keep


end % for iMC = 1:nMC

% MCMC results
results.mean = mean(modelSamples, 2);
results.meanSlope      = results.mean(1);
results.meanYIntercept = results.mean(2);
scatter(modelSamples(1,1:100:end), modelSamples(2,1:100:end), 'Marker', 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2, ...
    'MarkerEdgeColor','none', 'SizeData', 10)
hold on

%% theoretical OLS results

X = [setup.x, ones(setup.ndata,1)];
beta = X \ data;
betaCov = inv(X'*X);
% ellipse from OLS fit
thetas = linspace(0, 2*pi, 100); circpoints = 2*[cos(thetas); sin(thetas)];
cholcov = chol(betaCov, "lower"); 
ellipsepoints = cholcov*circpoints + beta;
plot(ellipsepoints(1,:), ellipsepoints(2,:), 'r', 'LineWidth', 2)
plot(beta(1), beta(2), '.r', 'MarkerSize', 25)


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