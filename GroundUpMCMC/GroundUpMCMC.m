%% Metropolis Hastings MCMC from the ground up

% 1. take the mean of a 1D dataset

n = 2;

rng(1) % start random number stream in one spot
data.d = random('normal', 0, 1, [n, 1]);
data.mean = mean(data.d);
data.stdv = std(data.d);
data.stde = data.stdv/sqrt(n);

m_current = data.mean; % initial model parameter
MCMC.proposalStdv = 1; % proposal distribution std dev
MCMC.nMC = 1e7; % number of MCMC trials

MCMC.m = nan([MCMC.nMC, 1], "double");
for iMC = 1:MCMC.nMC

    % save off current model
    MCMC.m(iMC) = m_current;
    
    % propose a new model
    m_propose = m_current + random("normal", 0, MCMC.proposalStdv);

    % calculate residuals for old and new models
    dhat_current = m_current * ones(n, 1);
    dhat_propose = m_propose * ones(n, 1);
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
plot(xvecForDist, pdf("normal", xvecForDist, MCMC.mean, sqrt(1/n)), '-r', 'LineWidth', 2)