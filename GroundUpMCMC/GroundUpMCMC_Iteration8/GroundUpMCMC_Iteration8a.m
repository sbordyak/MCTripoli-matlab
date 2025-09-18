% GroundUpMCMC_8a. Starting at 7b... 
% - introduce real methods from xml files
% - use the methods to generate synthetic data
% - write the synethetic data to a .txt file
% - read in the data from the .txt file

%% implement setup for this run

setupName = "Synthetic_TwoIsotopeStatic_1Mcps_Listric";
dataMode = "synthetic";

% massSpec = massSpecModel("PhoenixKansas_1e12");
% NBS982 = referenceMaterial("NBS982");
% 
% % unpack method from TIMSAM file
% method = parseTIMSAM(setup.methodsFolder + setup.methodName + ".TIMSAM");
% method = processMethod(method, massSpec.getCollectorNames);


%% If synthetic data, create a new dataset


if dataMode == "synthetic"
    
    % create mTrue, then specify method, mass spec, and run parameters
    setup = synDataSetup(setupName);

    % use setup to create synthetic data
    syndata = setup.makeSynData();


else % if real data

    % create dataVectors from input filename by parsing txt file

end

