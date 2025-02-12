function method = processMethod(method, collectorNames)
%PROCESSMETHOD Extract useful information from parsed method structure
%
%   Then append the useful info to the method structure

nSeq = size(method.onpeaks, 2);
OPnames = [method.onpeaks.Name]';
nFara = size(collectorNames,2);

% determine collector in axial position (Axial Faraday or Daly/Photomultiplier)
if string(method.settings.AxialColl) == "Axial"
    method.axialPositionDetector.Name = "Axial";
    method.axialPositionDetector.Code = "Ax";
elseif string(method.settings.AxialColl) == "PhotoMultiplier"
    method.axialPositionDetector.Name = "PhotoMultiplier";
    method.axialPositionDetector.Code = "PM";
else
    disp("Axial Position Detector Not Recognized!!")
end

% make blank table for OnPeak Sequence Table "OPTable"
OPT = table('Size', [nSeq, nFara], ...
                'VariableTypes', repelem("string", nFara), ...
                'VariableNames', collectorNames, ...
                'RowNames', OPnames);
OPTable = fillmissing(OPT, 'constant', ""); 
OPMasses = fillmissing(OPT, 'constant', "NaN"); clear OPT

% make blank table for OnPeak integration time table "OPIntegrationTimes"
OPIntegrationTiming = table('Size', [nSeq, 3], ...
    'VariableTypes', ["double", "double", "uint32"], ...
    'VariableNames', ["integrationPeriod", "integrationTime", "integrationsPerCycle"], ...
    'RowNames', OPnames);

% make OP table
for iSeq = 1:nSeq

    seqName = method.onpeaks(iSeq).Name;
    colArrayString = string(method.onpeaks(iSeq).CollectorArray);
    seqAssign = split(colArrayString, ",");
    seqAssign = split(seqAssign, ":");
    if size(seqAssign,2) == 1, seqAssign = seqAssign'; end %if one mass 

    % detectors names are concatenated with sequence name, eg "H3S1"
    activeCollectors = extractBefore(seqAssign(:,2), seqName)';
    
    % assign massID to OP sequence table
    nMasses = size(seqAssign, 1);
    for iMass = 1:nMasses
        OPTable.(activeCollectors(iMass))(seqName) = seqAssign(iMass,1);
    end

    % Parse and assign integration time to OPIntegrationTimes table. 
    % All times in seconds.
    integrationPeriodString = string(method.onpeaks(iSeq).IntegPeriod);
    integrationPeriod = double(extract(integrationPeriodString, digitsPattern))/1000;
    OPIntegrationTiming{OPnames(iSeq),'integrationPeriod'} = integrationPeriod;
    
    integrationTime = str2double(method.onpeaks(iSeq).IntegTime);
    OPIntegrationTiming{OPnames(iSeq),'integrationTime'} = integrationTime;

    OPIntegrationTiming{OPnames(iSeq),'integrationsPerCycle'} = ...
        uint32(integrationTime/integrationPeriod);

    % record approx isotopic masses (N values) in mass table for deltas
    for iMass = 1:nMasses
        massName = extract(seqAssign(iMass,1),  lettersPattern) + ...
        extract(seqAssign(iMass,1),  digitsPattern);
        OPMasses.(activeCollectors(iMass))(seqName) = double(mass.(massName));
    end % for iMass
    
end % for iOP

OPMasses = convertvars(OPMasses, collectorNames, 'double'); % convert to double


% calculate collector deltas (amu differences between collector positions)
AxialPositionDetectorIndex = find(collectorNames == method.axialPositionDetector.Code);

OPMassMatrix = table2array(OPMasses);
nMassesPerSequence = sum(~isnan(OPMassMatrix),2);
nMassPairs = sum(nMassesPerSequence.*(nMassesPerSequence-1)/2); % total pairs
massDiffMatrix = zeros(nMassPairs, nFara-1);
massDiffColumnIdcs = [1:AxialPositionDetectorIndex-1 0 AxialPositionDetectorIndex:nFara-1];
massDiffVector = zeros(nMassPairs,1);

massDiffRow = 0;
for iSeq = 1:nSeq

    OPMassColIdcs = find(~isnan(OPMassMatrix(iSeq,:)));
    for iStartIndex = 1 : (nMassesPerSequence(iSeq)-1)
        for jEndIndex = (iStartIndex+1) : nMassesPerSequence(iSeq)
            
            OPMassStartIndex = OPMassColIdcs(iStartIndex);
            OPMassEndIndex   = OPMassColIdcs(jEndIndex);

            massDiffRow = massDiffRow + 1;
            
            if OPMassStartIndex ~= AxialPositionDetectorIndex
                massDiffMatrix(massDiffRow,massDiffColumnIdcs(OPMassStartIndex)) = -1;
            end % if iStartIndex is not axial position
            if OPMassEndIndex ~= AxialPositionDetectorIndex
                massDiffMatrix(massDiffRow,massDiffColumnIdcs(OPMassEndIndex)) = 1;
            end % if jStartIndex is not axial position

            massDiff = OPMassMatrix(iSeq,OPMassEndIndex) - ...
                       OPMassMatrix(iSeq,OPMassStartIndex);
            massDiffVector(massDiffRow) = massDiff;

        end % for end mass index            

    end % for start mass index

end % for each sequence

warning('off', 'MATLAB:rankDeficientMatrix') % suppress warning
collectorDeltasFromAxPos = (massDiffMatrix\massDiffVector)';
warning('on', 'all'); % restore warning
detectorDeltas = [collectorDeltasFromAxPos(1:AxialPositionDetectorIndex-1) ...
                   0 collectorDeltasFromAxPos(AxialPositionDetectorIndex:end)];


% create a F_ind matrix for further data reduction
% first, find unique MassIDs (stings, label for unique isotopes)
OPTableStack = stack(OPTable, collectorNames);
MassIDs = table2array(unique(OPTableStack(:,2)));
MassIDs(MassIDs == "") = []; % delete blanks
% next, create F_ind with indexes to those MassIDs
F_ind = zeros(nSeq, nFara);
for iFara = 1:nFara
    FaraColumn = OPTable.(collectorNames(iFara));

    for iMassID = 1:length(MassIDs)

        F_ind(:,iFara) = F_ind(:,iFara) + (FaraColumn == MassIDs(iMassID)) * iMassID;

    end % for iMassID

end % for iMassID

% get axial masses for OP, from user-entered axial masses or peak centers
for iSeq = 1:nSeq

    axialMasses.OP(iSeq) = double(string(method.onpeaks(iSeq).AxMass)) + ...
                           double(string(method.onpeaks(iSeq).AxMassOffset));
                           % AxMass + AxMassOffset

end % for iSeq

% handle baselines if present
if isfield(method, 'baselines') % if baselines present

    nBL = size(method.baselines, 2);
    BLnames = string({method.baselines.Name})';
    BLTable = table('Size', [nBL, nFara], ...
        'VariableTypes', repelem("double", nFara), ...
        'VariableNames', collectorNames, ...
        'RowNames', BLnames);

    BLIntegrationTiming = table('Size', [nBL, 1], ...
        'VariableTypes', "double", ...
        'VariableNames', {'integrationTime'}, ...
        'RowNames', BLnames);

    for iBL = 1:nBL

        % save off baseline integration time
        integrationPeriodString = string(method.baselines(iBL).IntegPeriod);
        integrationPeriod = double(extract(integrationPeriodString, digitsPattern))/1000;
        BLIntegrationTiming{BLnames(iBL),'integrationPeriod'} = integrationPeriod;
    
        integrationTime = str2double(method.baselines(iBL).IntegTime);
        BLIntegrationTiming{BLnames(iBL),'integrationTime'} = integrationTime;

        BLIntegrationTiming{BLnames(iBL),'integrationsPerCycle'} = ...
            uint32(integrationTime/integrationPeriod);

        % if baseline mass defined by user-entered axial mass "AxMass"
        if method.baselines(iBL).BLReferences == "MASS"

            AxMass = double(string(method.baselines(iBL).AxMass)) + ...
                double(string(method.baselines(iBL).AxMassOffset));
            BLTable{iBL,collectorNames} = AxMass + detectorDeltas;
            axialMasses.BL(iBL) = AxMass;

        else % if baseline mass defined by offset from a peak-centered mass

            % user-generated reference to a peak ("Species:CollectorSequence")
            BLRefString = string(method.baselines(iBL).BLReferences);
            
            % determine sequence table position (eg H2S3), collector, and sequence name
            refSeqTablePosition = extractAfter(BLRefString,":");
            refCollector = extractBefore(refSeqTablePosition,3);
            refSequence = extractAfter(refSeqTablePosition,2);
            
            % find mass (in amu) referenced in OPMasses table, and its axial mass
            refMassInRefColl = table2array(OPMasses(refSequence,refCollector));
            refMassDistFromAxialAMU = detectorDeltas(collectorNames == refCollector);
            AxMassOffset = double(string(method.baselines(iBL).AxMassOffset));
            
            % AxMass during BL is AxMass during refSequence + AxMassOffset
            AxMass = refMassInRefColl - refMassDistFromAxialAMU + AxMassOffset;
            BLTable{iBL,collectorNames} = AxMass + detectorDeltas;
            axialMasses.BL(iBL) = AxMass;

        end % if baseline mass defined by "AxMass"

    end % for iBL

end

% save 
method.OPTable = OPTable;
method.OPMasses = OPMasses;
method.OPIntegrationTiming = OPIntegrationTiming;
method.F_ind = F_ind;
method.BLTable = BLTable;
method.BLIntegrationTimes = BLIntegrationTiming;
method.MassIDs = MassIDs;
method.detectorDeltas = detectorDeltas;
method.axialMasses = axialMasses;

end % function processMethod

