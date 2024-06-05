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
    m_propose = m_current + mvnrnd(zeros(nmodel,1), MCMC.proposalCov);

    % unpack parameters for maximum code clarity
    m_propose_slope = m_propose(1);
    m_propose_yintercept = m_propose(2);
    m_current_slope = m_propose(1);
    m_current_yintercept = m_current(2);

    % calculate residuals for old and new models
    dhat_current = m_current * ones(ndata, 1);
    dhat_propose = m_propose * ones(ndata, 1);
    r2_current = (data.d - dhat_current).^2; % squared residuals
    r2_propose = (data.d - dhat_propose).^2;
    ssr_current = sum(r2_current); % sum of the squared residuals = -loglikelihood
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

MCMC.mean = mean(MCMC.m);
MCMC.stdv = std(MCMC.m);

histogram(MCMC.m, "Normalization", "pdf")
hold on
xlimits = xlim(gca);
xvecForDist = linspace(xlimits(1), xlimits(2), 500);
%plot(xvecForDist, pdf("normal", xvecForDist, MCMC.mean, sqrt(1/n)), '-r', 'LineWidth', 2)
plot(xvecForDist, ...
    pdf("normal", xvecForDist, data.mean, sqrt(1/length(data.d))), '-g', 'LineWidth', 2)