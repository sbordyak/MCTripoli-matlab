% joing pdf for sum of guassian and poisson random variables
% from product of fourier transforms of pdfs (ie characteristic functions)

% lambda is Poisson rate parameter, sigmasq is Gaussian noise parameter
% sigma^2, x is domain of random variable (measured value), k is
% integration parameter (needed for inverse Fourier transform).
integrand = @(lambda, sigmasq, x, k) exp( lambda.*(exp(1i*k)-1) - sigmasq.*k.^2/2 - 1i*k.*x );



lambda = 4e7;
sigmasq = 6.24e6;
nx = 2000;
xvec = linspace(lambda-5*sqrt(lambda)-5*sqrt(sigmasq), lambda+5*sqrt(lambda)+5*sqrt(sigmasq), nx);
tic
px = zeros(nx,1);
for ix = 1:nx
    px(ix) = 1/(2*pi)*real(integral(@(k) integrand(lambda, sigmasq, xvec(ix), k), -Inf, Inf));
end 
toc
plot(xvec, px, 'LineWidth', 3)

% if xvec range is large enough, make sure we get a pdf with area = 1
assert(abs(trapz(xvec, px) - 1) < 1e-5, "did not return a valid pdf over range" )

%% Time the integration evaluation

%p2 = 1/(2*pi)*real(integral(@(k) integrand(10, 1, 10, k), -Inf, Inf));

pGP = @() 1/(2*pi)*real(integral(@(k) integrand(10, 1, 10, k), -Inf, Inf));
timeit(pGP)


% %% try just the Poisson portion:
% integrand2 = @(lambda, x, k) exp( lambda.*(exp(1i*k)-1) - 1i*k.*x );
% 
% xvec = 0:20;
% for ix = 1:21
%     p3(ix) = 1/(2*pi)*real(integral(@(k) integrand2(10, xvec(ix), k), -Inf, Inf));
% end
% 
% % numerically unstable!
