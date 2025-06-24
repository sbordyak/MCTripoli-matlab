function setup = getConfig(setup, method)
% GETCONFIG add additional setup fields for synthetic data
% setup and method structs contain info about analysis, including
% assumed values for parameters

% parameters for simulated measurement
setup.nBLIntegrations = 1e2;
setup.nOPIntegrations = 1e2;

% data and metadata needed specifically for synthetic data generation
setup.BLIntegrationTimes = ones(...
    setup.nBLIntegrations, method.BLIntegrationTimes.integrationTime);
setup.OPIntegrationTimes = ones(...
    setup.nOPIntegrations, method.OPIntegrationTiming.integrationTime);

setup.BLTimes = cumsum(setup.BLIntegrationTimes);
setup.OPTimes = max(setup.BLTimes) + 5 + cumsum(setup.OPIntegrationTimes);