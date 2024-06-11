%% Metropolis Hastings MCMC from the ground up

% 4. Generalize to an arbitrary model:
%    synthetic data reflects true values in struct 'truth'
%    restrict data to data vector = data, and G(m) = dhat
%    restrict data uncertainties to covariance matrix = dcov,
%    restrict model parameters to model vector = modelCurrent,
%       modelProposed, and modelSamples to aggregate samples
%    the rest of the parameters go in a struct = setup.
%    summary stats and interpretation go in results struct

% Start by reproducing linear regression case

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
dcov = eye(setup.ndata) * truth.yUncertainty^2;

modelCurrent = [setup.x, ones(setup.ndata,1)] \ data; % initial model parameter
setup.proposalCov = [1 0; 0 1]; % proposal distribution std dev
setup.nMC = 1e7; % number of MCMC trials

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
    
    % calculate log-likelihoods of current and proposed samples
    llCurrent  = loglik(dhatCurrent,  data, dcov);
    llProposed = loglik(dhatProposed, data, dcov);

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

%histogram(MCMC.m(1,:), "Normalization", "pdf")
%hold on
%xlimits = xlim(gca);
%xvecForDist = linspace(xlimits(1), xlimits(2), 500);
%plot(xvecForDist, pdf("normal", xvecForDist, MCMC.mean, sqrt(1/n)), '-r', 'LineWidth', 2)
%plot(xvecForDist, ...
%    pdf("normal", xvecForDist, data.mean, sqrt(1/length(data.d))), '-g', 'LineWidth', 2)


%% evaluate this particular model

function dhat = evaluateModel(m, setup)

    % y = slope * x       + y-intercept
    dhat = m(1) * setup.x + m(2);

end % function evaluateModel

%% calculate log-likelihood (including -1/2 term)
% ignore constant terms (for now)

function ll = loglik(dhat, data, dcov)

    residuals = (data - dhat);
    chiSqTerms = residuals.^2 ./ diag(dcov);
    ll = -1/2 * sum(chiSqTerms);

end % function loglik