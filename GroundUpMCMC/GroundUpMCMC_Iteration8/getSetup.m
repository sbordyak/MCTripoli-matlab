function setup = getSetup(setupName)
%GETSETUP Make setup struct of run parameters
%   Setup contains information not found in the data or method files

setup.dataMode = "synthetic"; % "synthetic" or "measured"

if setupName == "Synthetic_TwoIsotopeStatic_1Mcps_Listric"

    logmspl = [     % nseg = 2, bdeg = 3, t = 106:205
    14.4846439892590
    14.5251376191125
    14.4721275963024
    14.4292231240153
    14.4531727595166
    ];

end % if setupName

end % function
