function setup = getSetup(setupName, method)
%GETSETUP Make setup struct of run parameters
%   Setup contains information not found in the data or method files

setup.methodsFolder = "./massSpecMethods/";



switch setupName 

    case "Synthetic_TwoIsotopeStatic_1Mcps_Listric"

    % 100 cycles BL and OP for two-collector static 206-208 on two Faradays 
    % with a concave up, decaying beam described by by a cubic spline 
    % with two segments.

    setup.dataMode = "synthetic"; % "synthetic" or "measured"
    setup.methodName = "TwoCollectorStaticPb";

    setup.logmspl = [     % nseg = 2, bdeg = 3, t = 106:205
    14.4846439892590;
    14.5251376191125;
    14.4721275963024;
    14.4292231240153;
    14.4531727595166;
    ];

    setup.refIntensities = [-1e2; 2e2];
    setup.collRelativeEfficiencies = [1; 1];
    setup.betaFaraday = -0.2;
    setup.upMassTailFaraday = 0;
    setup.downMassTailFaraday = 0;




    


end % switch setupName

end % function
