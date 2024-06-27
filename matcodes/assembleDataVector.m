function [d, data] = assembleDataVector(data,method)
%ASSEMBLEDATAVECTOR Assemble data vector in d = g(m)
%   d is actually a struct containing the data vector and relevant
%   vectors of tags. Based off LoadMSdata_synth.m from Scott Burdick.
%   For Faraday Relative Efficiencies project, 27-Jun-2022 by Noah McLean
%   Revised 26-Jun-2024 by Noah McLean to add BL from ion counter(s)

% Problem: Isotopx method files may list a different set of detectors 
%   than in the data files. This means that the full sequence table columns 
%   from the method file don't necessarily match up with data file columns.
% Solution: determine which Faradays are used in the analysis, harmonize 
%   data/method matrices and indices so that columns line up with columns.
detNamesFromMethod = string(method.OPTable.Properties.VariableNames);
detNamesFromData = data.collectorNames;
[~, detectorIndicesInData] = ...
    intersect(detNamesFromData, detNamesFromMethod, 'stable');
[~, detectorIndicesInMethod] = ...
    intersect(detNamesFromMethod, detNamesFromData, 'stable');

BLTable = method.BLTable{:,detectorIndicesInMethod};
sequenceTableWithSpeciesIndices = method.F_ind(:,detectorIndicesInMethod);
data.BLmatrix = data.BLmatrix(:,detectorIndicesInData);
data.OPmatrix = data.OPmatrix(:,detectorIndicesInData);
data.CollectorTable = data.CollectorTable(detectorIndicesInData,:);


%% determine baseline references

nOPseq = size(sequenceTableWithSpeciesIndices,1);
BLdetIsRef = zeros(size(BLTable)); % true if baseline int is referenced
for iOPseq = 1:nOPseq

    % which baselines are referenced by each sequence?
    refString = string(method.onpeaks(iOPseq).BLReferences);
    %BLrefs is a vector with the BL squence indices referenced in iOPseq
    BLrefs = double(extractAfter(split(refString, ", "), "BL"));

    % make a flag for the BLTable with which intensities are used
    detInSeq = sequenceTableWithSpeciesIndices(iOPseq,:) > 0;
    BLdetIsRef(BLrefs, :) = BLdetIsRef(BLrefs, :) | detInSeq;

end
OPdetIsRef = sequenceTableWithSpeciesIndices > 0; % true if OP int is referenced


%% make masks for BL and OP data used in inversion

BLdataIsRef = logical(BLdetIsRef(data.BLSeqIdx,:));
OPdataIsRef = logical(OPdetIsRef(data.OPSeqIdx,:));


%% calcalute detector-specific parameters indexed to detector


% pull detector properties from CollectorTable at top of data file
isDetectorFaraday = data.CollectorTable(:,2) == "F";
detectorResistance = double(data.CollectorTable(:,3));

% calculate conversion from cps to volts 
ionsPerCoulomb = 6241509074460762607.776;
cpsPerVolt = ionsPerCoulomb ./ detectorResistance;

cpsConversionFactor = 1 .* not(isDetectorFaraday) + ...
    cpsPerVolt .* isDetectorFaraday;


%% assemble baseline data into data vector/struct

% initialize d -- structure with data (intensities)
% and metadata needed to evaluate model to calculate dhat
d.int   = []; % intensity
d.time  = []; % time
d.det   = []; % detector index
d.seq   = []; % sequence index
d.mass  = []; % isotope/species mass
d.iso   = []; % isotope index (for OP)
d.block = []; % block index
d.cycle = [];
d.isOP  = false(0); % is data point an OP measurement?

nDet = size(sequenceTableWithSpeciesIndices,2);

for iDet = 1:nDet

    detRefs = BLdataIsRef(:,iDet); % where is this detector refd in BL?
    nrefs   = sum(detRefs);
    d.int   = [d.int; data.BLmatrix(detRefs, iDet)*cpsConversionFactor(iDet)];
    d.time  = [d.time; data.BLtime(detRefs)];
    d.det   = [d.det; iDet*ones(nrefs,1)];
    BLseqs  = data.BLSeqIdx(detRefs); % BL seqs referenced
    d.seq   = [d.seq; BLseqs];
    d.mass  = [d.mass; BLTable( BLseqs, iDet )];
    d.iso   = [d.iso; zeros(nrefs,1)];
    d.block = [d.block; data.BLserial(detRefs, 1)];
    d.cycle = [d.cycle; zeros(nrefs,1)];
    d.isOP  = [d.isOP; false(nrefs,1)];

end


%% append on-peak data onto d vector/struct

for iDet = 1:nDet

    detRefs = OPdataIsRef(:,iDet); % where is this detector refd in F_ind?
    nrefs   = sum(detRefs);
    d.int   = [d.int; data.OPmatrix(detRefs, iDet)*cpsConversionFactor(iDet)];
    d.time  = [d.time; data.OPtime(detRefs)];
    d.det   = [d.det; iDet*ones(nrefs,1)];
    OPseqs  = data.OPSeqIdx(detRefs); % OP seqs referenced
    d.seq   = [d.seq; OPseqs];
    d.mass  = [d.mass; method.OPMasses{ OPseqs, iDet }];
    d.iso   = [d.iso; sequenceTableWithSpeciesIndices(OPseqs,iDet)];
    d.block = [d.block; data.OPserial(detRefs, 1)];
    d.cycle = [d.cycle; data.OPserial(detRefs, 2)];
    d.isOP  = [d.isOP; true(nrefs,1)];

end



end % function