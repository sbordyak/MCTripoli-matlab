%% Metropolis Hastings MCMC from the ground up

% 6. apply version 6 to superposition model for Erick's Eastern Shelf U-Pb
% ages

%% data setup

%      ABBA-3C  NoPlat PAM-2F WFBU10E CFC-3D
data = [296.0,  291.9  289.0  284.7   283.5]';
data1s = [6.7,    6.0    4.2    6.9     5.6]'/2;
dvar = data1s.^2;

setup.ndata = length(data);

rng(1) % start random number stream in one spot


%% initialize model parameters and likelihoods

% model parameters:
% m(1) = date of ABBA-3C
% m(2) = delta-t from ABBA-3C to NoPlat
% m(3) = detla-t from NoPlat to PAM-2F
% m(4) = delta-t from PAM-2F to WFBU10E
% m(5) = delta-t from WFBU10E to CFC-3D
% so that
% d(1) = m(1)
% d(2) = m(1) - m(2)
% d(3) = m(1) - m(2) - m(3)
% d(4) = m(1) - m(2) - m(3) - m(4)
% d(5) = m(1) - m(2) - m(3) - m(4) - m(5)

modelCurrent = zeros(5,1);
modelCurrent(1) = data(1);
modelCurrent(2:end) = data(1:end-1) - data(2:end);

% initialize data variance vector, dhat vector, and log-likelihood
dhatCurrent = evaluateModel(modelCurrent);
llCurrent = loglik(dhatCurrent, data, dvar, modelCurrent);

setup.proposalCov = 5*eye(5);
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
    dhatProposed = evaluateModel(modelProposed);
    
    % create data covariance with current and proposed noise terms
    % dvarProposed = updateDataVariance(modelProposed, setup);
    
    % calculate log-likelihoods of current and proposed samples
    llProposed = loglik(dhatProposed, data, dvar, modelProposed);

    % difference in log-likelihood = log(likelihood ratio)
    delta_ll = llProposed - llCurrent;
    
    % probability of keeping the proposed model
    keep = min(1, exp(delta_ll)); % negative one-half from likelihood eqn

    % keep the proposed model with probability = keep
    if keep > rand(1) % if we're keeping the proposed model

        % the proposed model becomes the current one
        modelCurrent = modelProposed; 
        dhatCurrent = dhatProposed;
        % dvarCurrent = dvarProposed;
        llCurrent  = llProposed;

    end % if keep


end % for iMC = 1:nMC

plotmatrix(modelSamples')



%% MCMC results

% results.mean = mean(modelSamples, 2);
% results.meanSlope      = results.mean(1);
% results.meanYIntercept = results.mean(2);
% % correlation plot
% figure('Position', [50, 50, 800, 600])
% plotmatrix(modelSamples')
% 
% figure('Position', [900 50, 700, 600])
% scatter(modelSamples(1,:), modelSamples(2,:), 'Marker', 'o', ...
%     'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2, ...
%     'MarkerEdgeColor','none', 'SizeData', 10)
% hold on
% 
% 
% %% theoretical OLS results
% 
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
% 
% %% plot data and sample fits
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


datehat = [modelSamples(1,:);
           modelSamples(1,:) - modelSamples(2,:);
           modelSamples(1,:) - modelSamples(2,:) - modelSamples(3,:);
           modelSamples(1,:) - modelSamples(2,:) - modelSamples(3,:) - modelSamples(4,:);
           modelSamples(1,:) - modelSamples(2,:) - modelSamples(3,:) - modelSamples(4,:) - modelSamples(5,:)];

lb = datehat(2,:); % noPlat
ub = datehat(3,:); % PAM-2F
fracDist = rand(1, size(datehat,2));
WLB = lb - fracDist .* (lb - ub);

%% evaluate this particular model

function dhat = evaluateModel(m)

    % y = slope * x       + y-intercept
    dhat = m(1) - cumsum([0; m(2:end)]);

end % function evaluateModel


%% update data covariance matrix

function dvar = updateDataVariance(m, setup)

    dvar = m(3) * ones(setup.ndata, 1);

end % updateDataCovariance


%% calculate log-likelihood (including -1/2 term)
% ignore constant terms (for now)

function ll = loglik(dhat, data, dvar, m)
    
    % enforce bounds on model prior
    if any(m < 0)
        ll = -Inf;
        return
    end

    residuals = (data - dhat);
    chiSqTerms = residuals.^2 ./ dvar;
    ll = -1/2 * sum(chiSqTerms); % - 1/2 * sum(log(dvar));

end % function loglik