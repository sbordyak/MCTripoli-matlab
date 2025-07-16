classdef dataVectors
    %DATASET Measured mass spectrometer data
    %   Includes constructor and vectorizing methods
    
    properties
        intensity       (:,1) double  = []
        time            (:,1) double  = []
        integrationTime (:,1) double  = []
        detectorIndex   (:,1) uint8   = []
        sequenceIndex   (:,1) uint8   = []
        mass            (:,1) double  = []
        isotope         (:,1) uint8   = []
        block           (:,1) uint16  = []
        cycle           (:,1) uint32  = []
        isOP            (:,1) logical = []
        isUsed          (:,1) logical = []
        header
        method
    end
    
    methods
        function obj = dataVectors(varargin)

        if nargin > 0 % if constructed with inputs
            dataFolder = varargin(1);
            dataFolder = dataFolder{1};

            % Scott: replicate assembleDataVector(data, method) here
            % to create a dataset


            textFileInfo = dir(dataFolder + "/*.TXT");

            %% parse header block

            % get method row from header, find method name
            opts = delimitedTextImportOptions('NumVariables', 2);
            opts.DataLines = [4, 11];
            methodHeader = readcell(textFileInfo.name, opts);

            obj.header.fileName   = string(methodHeader{1,2});
            obj.header.methodName = string(methodHeader{2,2});
            obj.header.methodPath = string(methodHeader{3,2});
            obj.header.IsoWorksMethod = string(methodHeader{4,2});
            obj.header.FolderPath = string(methodHeader{5,2});
            obj.header.Corrected = string(methodHeader{6,2});
            obj.header.BChannels = string(methodHeader{7,2});
            obj.header.TimeZero = string(methodHeader{8,2});

            
            if ~strcmpi(obj.header.FolderPath,dataFolder)
                obj.header.FolderPath = dataFolder;
            end

        end
        end
           
        
        function obj = loadMethodFile(obj,collectorNames,varargin)
        
        if nargin > 2 % if constructred with inputs

            methodfilename = varargin{1};
        else
            methodfilename = [obj.header.methodPath '\' obj.header.methodName];
        end

        obj.method = parseTIMSAM(methodfilename);
        obj.method = processMethod(obj.method, collectorNames);

        end


        function tmp = loadDataFile(obj)


            textFileInfo = dir(obj.header.FolderPath + "/*.TXT");

            %% new "parse data" required for ATONA A/B channel data

            opts = delimitedTextImportOptions;
            opts.VariableTypes = "string";
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";
            firstColumn = readmatrix(textFileInfo.name, opts);
            % note that default is to skip empty lines
            collectorsStartPosition = find(firstColumn == "#COLLECTORS");
            userTablesStartPosition = find(firstColumn == "#USERTABLES");
            baselinesStartPosition  = find(firstColumn == "#BASELINES");
            onPeakStartPosition     = find(firstColumn == "#ONPEAK");
            endPosition             = find(firstColumn == "#END");

            collectorsEndPosition = min([userTablesStartPosition baselinesStartPosition onPeakStartPosition]) - 2;
            nCollectors = collectorsEndPosition - collectorsStartPosition - 1;
            tmp.collectorNames = firstColumn(collectorsStartPosition+2:collectorsEndPosition);

            if obj.header.BChannels == "No" % if resistor-based amplifiers, no BChannels
                nDataColumns = 7 + nCollectors;
            elseif obj.header.BChannels == "Yes" % if ATONAs
                nDataColumns = 7 + 2*nCollectors - 1; % minus one because PM doesn't have BChannel
            else
                disp('unrecognized text file column setup')
            end

            % extract collector block information
            opts.VariableNames = string(1:6);
            opts.VariableTypes = ["string", "string", "string", "string", "string", "string"];
            opts.DataLines = [collectorsStartPosition+1 collectorsEndPosition];
            tmp.Collectors = readmatrix(textFileInfo.name, opts);
            tmp.CollectorTable = tmp.Collectors(2:end,:);

            % extract Faraday gains and resistances -- not currently used
            firstFaradayRow = find(tmp.Collectors(:,2) == "F", 1, "first");
            lastFaradayRow = find(tmp.Collectors(:,2) == "F", 1, "last");
            Frange = firstFaradayRow:lastFaradayRow;
            tmp.FaradayResist = double(tmp.Collectors(Frange,3)); % resistances (ohms)
            tmp.FaradayGains = double(tmp.Collectors(Frange,4));  % gains (relative to Axial)

            % next: line 104 (parse BL and OP data, separate out indices and data)

            % baselines
            if ~isempty(baselinesStartPosition) % if baselines exist
                opts.VariableNames = string(1:nDataColumns);
                opts.VariableTypes = repmat("string", 1, nDataColumns);
                opts.DataLines = [baselinesStartPosition+2 onPeakStartPosition-2];
                tmp.BLall = readmatrix(textFileInfo.name, opts);

                tmp.BLserial = double(tmp.BLall(:,2:4));% [block cycle integration] serially assigned counts

                tmp.BLmatrix = double(tmp.BLall(:,8:end)); % matrix of collector readings
                if obj.header.BChannels == "Yes"
                    tmp.BLmatrixA = tmp.BLmatrix(:,1:nCollectors);
                    tmp.BLmatrixB = tmp.BLmatrix(:,nCollectors+1:end);
                end

                tmp.BLtime   = double(tmp.BLall(:,7)); % time
                tmp.BLID = tmp.BLall(:,1); % baseline ID, eg "BL1", "BL2", etc. 1st column in TXT data file
                tmp.BLSeqIdx = double(extractAfter(tmp.BLall(:,1), "BL"));
            end % if baselines exist

            % on-peaks
            opts.DataLines = [onPeakStartPosition+2 endPosition-2];
            tmp.OPall = readmatrix(textFileInfo.name, opts);

            tmp.OPserial = double(tmp.OPall(:,2:4));% [block cycle integration] serially assigned counts

            tmp.OPmatrix = double(tmp.OPall(:,8:end)); % matrix of all collector readings
            if obj.header.BChannels == "Yes" % if ATONAs
                tmp.OPmatrixA = tmp.OPmatrix(:,1:nCollectors);
                tmp.OPmatrixB = tmp.OPmatrix(:,nCollectors+1:end);
            end

            tmp.OPtime   = double(tmp.OPall(:,7)); % time
            tmp.OPID = tmp.OPall(:,1); % OnPeak ID, eg "OP1", "OP2", etc.  1st column in TXT data file
            tmp.OPSeqIdx = double(extractAfter(tmp.OPall(:,1), "OP"));

        end


            %%%%%%%%%%%%%%%%%%
            function obj = formDataVectors(obj)

            tmp = loadDataFile(obj);


            detNamesFromMethod = string(obj.method.OPTable.Properties.VariableNames);
            detNamesFromData = tmp.collectorNames;
            [~, detectorIndicesInData] = ...
                intersect(detNamesFromData, detNamesFromMethod, 'stable');
            [~, detectorIndicesInMethod] = ...
                intersect(detNamesFromMethod, detNamesFromData, 'stable');

            BLTable = obj.method.BLTable{:,detectorIndicesInMethod};
            sequenceTableWithSpeciesIndices = obj.method.F_ind(:,detectorIndicesInMethod);
            tmp.BLmatrix = tmp.BLmatrix(:,detectorIndicesInData);
            tmp.OPmatrix = tmp.OPmatrix(:,detectorIndicesInData);
            tmp.CollectorTable = tmp.CollectorTable(detectorIndicesInData,:);


            %% determine baseline references

            nOPseq = size(sequenceTableWithSpeciesIndices,1);
            BLdetIsRef = zeros(size(BLTable)); % true if baseline int is referenced
            for iOPseq = 1:nOPseq

                % which baselines are referenced by each sequence?
                refString = string(obj.method.onpeaks(iOPseq).BLReferences);
                %BLrefs is a vector with the BL squence indices referenced in iOPseq
                BLrefs = double(extractAfter(split(refString, ", "), "BL"));

                % make a flag for the BLTable with which intensities are used
                detInSeq = sequenceTableWithSpeciesIndices(iOPseq,:) > 0;
                BLdetIsRef(BLrefs, :) = BLdetIsRef(BLrefs, :) | detInSeq;

            end
            OPdetIsRef = sequenceTableWithSpeciesIndices > 0; % true if OP int is referenced


            %% make masks for BL and OP data used in inversion

            BLdataIsRef = logical(BLdetIsRef(tmp.BLSeqIdx,:));
            OPdataIsRef = logical(OPdetIsRef(tmp.OPSeqIdx,:));


            %% calcalute detector-specific parameters indexed to detector


            % pull detector properties from CollectorTable at top of data file
            isDetectorFaraday = tmp.CollectorTable(:,2) == "F";
            detectorResistance = double(tmp.CollectorTable(:,3));

            % calculate conversion from cps to volts
            ionsPerCoulomb = 6241509074460762607.776;
            cpsPerVolt = ionsPerCoulomb ./ detectorResistance;

            cpsConversionFactor = 1 .* not(isDetectorFaraday) + ...
                cpsPerVolt .* isDetectorFaraday;


            %% assemble baseline data into data vector/struct

            % % initialize d -- structure with data (intensities)
            % % and metadata needed to evaluate model to calculate dhat
            % d.int   = []; % intensity
            % d.time  = []; % time
            % d.det   = []; % detector index
            % d.seq   = []; % sequence index
            % d.mass  = []; % isotope/species mass
            % d.iso   = []; % isotope index (for OP)
            % d.block = []; % block index
            % d.cycle = [];
            % d.isOP  = false(0); % is data point an OP measurement?

            nDet = size(sequenceTableWithSpeciesIndices,2);

            for iDet = 1:nDet

                detRefs = BLdataIsRef(:,iDet); % where is this detector refd in BL?
                nrefs   = sum(detRefs);
                obj.intensity   = [obj.intensity; tmp.BLmatrix(detRefs, iDet)*cpsConversionFactor(iDet)];
                obj.time  = [obj.time; tmp.BLtime(detRefs)];
                obj.detectorIndex   = [obj.detectorIndex; iDet*ones(nrefs,1)];
                BLseqs  = tmp.BLSeqIdx(detRefs); % BL seqs referenced
                obj.sequenceIndex   = [obj.sequenceIndex; BLseqs];
                obj.mass  = [obj.mass; BLTable( BLseqs, iDet )];
                obj.isotope   = [obj.isotope; zeros(nrefs,1)];
                obj.block = [obj.block; tmp.BLserial(detRefs, 1)];
                obj.cycle = [obj.cycle; zeros(nrefs,1)];
                obj.isOP  = [obj.isOP; false(nrefs,1)];
                obj.isUsed  = [obj.isUsed; true(nrefs,1)];
                obj.integrationTime = [obj.integrationTime; obj.method.BLIntegrationTimes.integrationPeriod(1)*ones(nrefs,1)];

            end


            %% append on-peak data onto d vector/struct

            for iDet = 1:nDet

                detRefs = OPdataIsRef(:,iDet); % where is this detector refd in F_ind?
                nrefs   = sum(detRefs);
                obj.intensity   = [obj.intensity; tmp.OPmatrix(detRefs, iDet)*cpsConversionFactor(iDet)];
                obj.time  = [obj.time; tmp.OPtime(detRefs)];
                obj.detectorIndex   = [obj.detectorIndex; iDet*ones(nrefs,1)];
                OPseqs  = tmp.OPSeqIdx(detRefs); % OP seqs referenced
                obj.sequenceIndex   = [obj.sequenceIndex; OPseqs];
                obj.mass  = [obj.mass; obj.method.OPMasses{ OPseqs, iDet }];
                obj.isotope   = [obj.isotope; sequenceTableWithSpeciesIndices(OPseqs,iDet)];
                obj.block = [obj.block; tmp.OPserial(detRefs, 1)];
                obj.cycle = [obj.cycle; tmp.OPserial(detRefs, 2)];
                obj.isOP  = [obj.isOP; true(nrefs,1)];
                obj.isUsed  = [obj.isUsed; true(nrefs,1)];
                obj.integrationTime = [obj.integrationTime; obj.method.OPIntegrationTiming.integrationPeriod(OPseqs)];

            end

    



        % end % if constructed from a data input
        % 
         end
       


    end
end

