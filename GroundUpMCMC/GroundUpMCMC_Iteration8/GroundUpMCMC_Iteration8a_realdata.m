% GroundUpMCMC_8a. Starting at 7b... 
% - introduce real methods from xml files
% - use the methods to generate synthetic data
% - write the synethetic data to a .txt file
% - read in the data from the .txt file

%% retrieve setup for this run

%setupName = "Synthetic_TwoIsotopeStatic_1Mcps_Listric";

%setup = getSetup(setupName);
%myReferenceMaterial = referenceMaterial("NBS982");

% unpack method from TIMSAM file
%method = parseTIMSAM(setup.methodsFolder + setup.methodName + ".TIMSAM");
%method = processMethod(method, massSpec.getCollectorNames);

massSpec = massSpecModel("PhoenixKansas_1e12");
dataFolder = "~/Documents/GitHub/TripoliTestData/IsotopxPhoenixTIMS/KU_IGL/IsolinxVersion2/NBS981 230024b.RAW";
methodFile = "Pb NBS981 204-5-6-7-8 Daly 14-1-5-5-2 sec.TIMSAM";

d = dataVectors(dataFolder);
d = loadMethodFile(d, massSpec.getCollectorNames,methodFile);
d = formDataVectors(d);

return
%% If synthetic data, create a new dataset


if setup.dataMode == "synthetic"
    
    % an additional configuration info is required for synthetic data
    setup = getConfig(setup, method);

    % populate unknowns object for synthetic data creation
    mTrue = modelParameterSet(myReferenceMaterial, setup);
    
    

end

