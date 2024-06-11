%% Metropolis Hastings MCMC from the ground up

% 3. Try linear regression, start with y-only errors

ndata = 20;
nmodel = 2; % number of model parameters
model.slope = 2;
model.yintercept = 4;


rng(1) % start random number stream in one spot
data.yUncertainty = 1; % standard deviation of normally distributed errors

data.x = random('uniform', 0, 10, [ndata, 1]);
data.yTrue = model.slope * data.x + model.yintercept;
data.y = random('normal', 0, data.yUncertainty, [ndata, 1]) + data.yTrue;

m_current = [data.x, ones(ndata,1)] \ data.y; % initial model parameter
MCMC.proposalCov = [1 0; 0 1]; % proposal distribution std dev
MCMC.nMC = 1e7; % number of MCMC trials

%%

MCMC.m = nan([nmodel, MCMC.nMC], "double");
for iMC = 1:MCMC.nMC

    % save off current model
    MCMC.m(:, iMC) = m_current;
    
    % propose a new model
    m_propose = m_current + mvnrnd(zeros(nmodel,1), MCMC.proposalCov)';

    % unpack parameters for maximum code clarity
    m_propose_slope = m_propose(1);
    m_propose_yintercept = m_propose(2);
    m_current_slope = m_current(1);
    m_current_yintercept = m_current(2);

    % calculate residuals for old and new models
    dhat_current = m_current_slope * data.x + m_current_yintercept;
    dhat_propose = m_propose_slope * data.x + m_propose_yintercept;
    
    % squared weighted residuals for likelihood function
    r2_current = (data.y - dhat_current).^2 / data.yUncertainty^2; 
    r2_propose = (data.y - dhat_propose).^2 / data.yUncertainty^2;
    
    % sum of the squared residuals \propto -loglikelihood
    ssr_current = sum(r2_current); 
    ssr_propose = sum(r2_propose);

    % difference in ssr (-log of likelihood ratio, neglecting 1/2 factor)
    delta_ssr = ssr_propose - ssr_current;
    
    % probability of keeping the proposed model
    keep = min(1, exp(-delta_ssr/2)); % negative one-half from likelihood eqn

    % keep the proposed model with probability = keep
    if keep > rand(1) % if we're keeping the proposed model

        m_current = m_propose; % the proposed model becomes the current one

    end % if keep


end % for iMC = 1:nMC

% MCMC results
MCMC.mean = mean(MCMC.m, 2);
MCMC.meanSlope      = MCMC.mean(1);
MCMC.meanYIntercept = MCMC.mean(2);
scatter(MCMC.m(1,1:100:end), MCMC.m(2,1:100:end), 'Marker', 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.2, ...
    'MarkerEdgeColor','none', 'SizeData', 10)
hold on

%% theoretical OLS results

X = [data.x, ones(ndata,1)];
beta = [data.x, ones(ndata,1)] \ data.y;
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