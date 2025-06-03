% GroundUpMCMC. Starting at 7b... 
% - introduce real methods from xml files
% - use the methods to generate synthetic data
% - write the synethetic data to a .txt file
% - read in the data from the .txt file

%% retrieve setup for this run

setupName = "Synthetic_TwoIsotopeStatic_1Mcps_Listric";

setup = getSetup(setupName);
massSpec = massSpecModel("PhoenixKansas_1e12");
mySample = referenceMaterial("NBS982");

%% If synthetic data, create a new dataset


if setup.dataMode == "synthetic"
    
    % unpack method from TIMSAM file
    method = parseTIMSAM(setup.methodsFolder + setup.methodName + ".TIMSAM");
    method = processMethod(method, massSpec.getCollectorNames);

    % populate unknowns object for synthetic data creation
    mTrue = modelParameters(mySample, setup);

end

