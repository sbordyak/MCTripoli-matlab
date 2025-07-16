%clear

%addpath(genpath("./data"))
addpath(genpath('./matcodes'))
addpath(genpath("~/Documents/GitHub/TripoliTestData/"))

%addpath(genpath('~/Documents/GitHub/MassSpecMCMC-matlab/'))

%%
% 1. input a filename from data folder


%dataFolder = "~/Documents/GitHub/TIMSLAB/CollectorRelativeEfficiencies/data/Pb/C054_Pb.RAW";
dataFolder = "~/Documents/GitHub/TripoliTestData/IsotopxPhoenixTIMS/KU_IGL/IsolinxVersion2/NBS981 230024b.RAW";

% 2. parse the data file


data = parseTXT(dataFolder);
if data.header.BChannels == "Yes" % if ATONAs
    [data, dataB] = parseBChannels(data);
end

% 3. grab the corresponding methods file, make run tables for OP and BL

method = parseTIMSAM(data.header.methodName);

DetNames = ["PM" "L5", "L4", "L3", "L2", "Ax", "H1", "H2", "H3", "H4"];
method = processMethod(method, DetNames);


% 4. assemble data vector and tags

[d00, data] = assembleDataVector(data, method);
% also trims data matrices to collectors from method

return

d00.det = d00.det-1;
d00.det(d00.det==0) = 10;

DetNames = DetNames([2:end 1]);
DetNamesUsed = DetNames(unique(d00.det));

Iso_Name =method.MassIDs;




% 5. Create output folder
runname = 'test_block1';
iset = 1;


foldername = ['./results/' runname];

outfolder = ['./results/' runname '/output/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end



%%



for iblock = 11


    d0 = getblock_MCT(iblock,d00);
    
    %d0 = LoadMSdata_synth(datafile,Isotopes,F_ind);

    % Matrix to project spline knots for intensity function to time samples
    %InterpMat = d0.InterpMat;





    %% Set Up Adaptive MCMC INVERSION
    % Hardcoded here, presumably set by user in file or GUI?

    % MCMC Parameters
    maxcnt = 5000;  % Maximum number of models to save
    hier = 0;  % Hierachical?
    datsav = 10;  % Save model every this many steps

    burn = 1;  % Burn-in, start doing stats after this many saved models

    % Create tempering vector - start high, cool down to 1 then stay there
    temp=1; % Use tempering?
    Ntemp = 10000; % Cool search over this number of steps
    TT = ones(maxcnt*datsav,1);TT(1:Ntemp) = linspace(1,1,Ntemp)';

    % Baseline multiplier - weight Daly more strongly (I think)
    blmult = ones(size(d0.data));
    blmult(d0.axflag)=0.1;

    %sb629 Commented out unused variables
    %Ndata=d0.Ndata; % Number of picks
    %Nsig = d0.Nsig; % Number of noise variables


    % Range for ratios and intensity parameters
    %sb629  Changed priors to infinite where appropropriate
    prior.BL = [-inf inf];  % Faraday baseline
    prior.BLdaly = [0 0];   % Daly baseline (no baseline uncertainty)
    prior.lograt = [-100 100]; % Log ratio
    prior.I = [-inf inf];  % Intensity
    prior.DFgain = [-inf inf];  % Daly-Faraday gain

    prior.sig = [0 1e6];  % Noise hyperparameter for Faraday
    prior.sigdaly = [0 0]; % Gaussian noise on Daly
    prior.sigpois = [0 10]; % Poisson noise on Daly


    psig = []; %






    %% CREATE INITIAL MODEL AND MODEL DATA
    %[x0,d,Intensity] = InitializeModel_synth(d0);

    % Hardcoded variables for whatever I'm testing at the moment
    %d0.ReportInterval = 0.1;
    user_DFGain = 0.9;

    % Hardcoded to test removing first cycle %sb726
    %d0.Include = ~(d0.cycle == 1 & d0.block == 1);
    %d0.IncludeMat{1} = ~(d0.CycleMat{1} == 1);

    [x0,C0,d] = InitializeModelCov_synth(d0,prior,user_DFGain);


x0.Dsig(d0.axflag) = x0.Dsig(d0.axflag)+(2e6);
%x0.Dsig = x0.Dsig;

    %C0(4:end-5,:)=C0(4:end-5,:);
    
    Nmod = size(C0,1);  %sb629 Slightly shorter way to define size of model.

    Dsig = x0.Dsig; % Using fixed noise variance values

    % New function to compute and plot ratios by cycle %sb629
    ViewCycleRatios(x0,d0,Iso_Name)

    % Assign initial values for model x
    x=x0;

    PlotData_MS


    %% Initialize Convergence Criteria Variables
    beta = 0.05;

    % Modified Gelman-Rubin Convergence
    alpha=0.025; % Confidence interval we want to be accurate (0.05 = 95% CI)
    epsilon=0.025; % Relative confidence in mean compared to std dev estimator(?)
    EffectSamp = 10000;%2^(2/Nmod)*pi/(Nmod*gamma(Nmod/2))^(2/Nmod)*chi2inv(1-alpha,Nmod)/epsilon^2;
    Mchain = 1; % Number of Chains
    ExitCrit = sqrt(1+Mchain/EffectSamp); % Exit when G-R criterium less than this


    % Initialize data residual vectors
    restmp=zeros(size(Dsig));
    restmp2=zeros(size(Dsig));

    % Calculate data residuals from starting model
    restmp = (d0.data-d).^2;
    restmp_weight = restmp.*blmult./Dsig;

    % Calculate error function %sb726
    E=sum(restmp_weight(d0.Include)/TT(1));  % Weighted by noise variance (for acceptance)
    E0=sum(restmp(d0.Include));  % Unweighted (for tracking convergence)


    %% Initialize MCMC loop variables
    cnt=0; % Counter
    cnt2=0;
    kept=zeros(5,4); % For displaying how many updates are accepted

    clear ens*
    ensemble=[]; % Make sure to start with new ensemble



    % Data and data covariance vectors
    xmean = zeros(Nmod,1);
    xcov = zeros(Nmod,Nmod);

    % Adaptive MCMC proposal term
    delx_adapt=zeros(Nmod,datsav);

    %%
    %d0.iso_vec(d0.iso_vec==0)=d0.Niso; %Set BL to denominator iso

% SB

    %% MCMC Iterations
    tic


    for m = 1:maxcnt*datsav
        %%
        % Choose an operation for updating model
        oper = RandomOperMS(hier);


        clear delx

        adaptflag = 1;
        allflag = 1;
        temp = 1;

        if m<=2*Nmod   % Use initial covariance until 2*N
            C = C0;
        else % After that begin updating based on model covariance

            % Next proposal based initial variance and iterative covariance
            C = beta*C0 + (1-beta)*2.38^2*Nmod^-1*xcov;
            C=(C'+C)/2; % Make sure it's symmetrical
        end

        diagC = diag(C);
        Cscale = C./sqrt(diagC*diagC');

        % Draw random numbers based on covariance for next proposal
        delx_adapt = mvnrnd(zeros(Nmod,1),Cscale)';

        delx_adapt = delx_adapt.*diagC;

        %delx_adapt(20:end)=0;

        % Update model and save proposed update values (delx)
        [x2,delx] = UpdateMSv2(oper,x,psig,prior,ensemble,xcov,delx_adapt,adaptflag,allflag);



        d2 = ModelMSData(x2,d0);

        % New data covariance vector
        %Dsig2 = x2.sig(d0.det_vec).^2 + x2.sig(d0.iso_vec+d0.Ndet).*dnobl2;
        Dsig2 = Dsig;

        % Calculate residuals for current and new model
        restmp = (d0.data-d).^2;
        restmp2 = (d0.data-d2).^2;

        restmp_weight = restmp.*blmult./Dsig;
        restmp2_weight = restmp2.*blmult./Dsig2;

        %sb726
        E02=sum(restmp2(d0.Include));  % Unweighted error func (for visualization)

        %sb726
        if strcmp(oper,'noise')
            % If noise operation
            E=sum(restmp_weight(d0.Include));
            E2=sum(restmp2_weight(d0.Include));
            dE=E2-E; % Change in misfit
        else
            % If any other model update
            E=sum(restmp_weight(d0.Include)/TT(m));
            E2=sum(restmp2_weight(d0.Include)/TT(m));
            dE=E2-E; % Change in misfit
        end



        %%
        % Decide whether to accept or reject model
        keep = AcceptItMS(oper,dE,psig,delx,prior,Dsig,Dsig2,d0);

        % Update kept variables for display
        kept(OpNumMS(oper),2) = kept(OpNumMS(oper),2)+1;
        kept(OpNumMS(oper),4) = kept(OpNumMS(oper),4)+1;

        % If we accept the new model update values
        if keep>=rand(1)

            E=E2; % Misfit
            E0=E02; % Unweighted misfit
            d=d2; % Data
            x=x2; % Model
            Dsig=Dsig2;  % Model variance
            %        dnobl=dnobl2;  % Data without baseline
            %        Intensity=Intensity2;  % Intensity

            % Display info
            kept(OpNumMS(oper),1) = kept(OpNumMS(oper),1)+1;
            kept(OpNumMS(oper),3) = kept(OpNumMS(oper),3)+1;

        end


        [xmean,xcov] = UpdateMeanCovMS(x,xmean,xcov,m);

        %figure(99);imagesc(xcov);colorbar;drawnow

        % Save model values to ensemble at regular intervals
        if  mod(m,datsav)==0

            cnt=cnt+1; % Increment counter

            ensemble(cnt).lograt=x.lograt; % Log ratios
            for mm=1:d0.Nblock
                ensemble(cnt).I{mm}=x.I{mm}; % Intensity by block
            end
            ensemble(cnt).BL=x.BL;  % Baselines
            ensemble(cnt).DFgain=x.DFgain;  %Daly-Faraday gain
            ensemble(cnt).sig= 1; %x.sig;  % Noise hyperparameter (Calc at beginning now)
            ensemble(cnt).E=E;  % Misfit
            ensemble(cnt).E0=E0; % Unweighted misfit
            ensemble(cnt).C = C;
            ensemble(cnt).delx=delx_adapt;
            ensemble(cnt).dd = d2-d;


            % Display update info to screen
            if  mod(m,100*datsav)==0
                DisplayInfo
                tic
                kept(:,1:2) = 0;
            end


            % If you want, calculate stats and plot regularly during iterations
            if  mod(m,100*datsav)==0
                %%
                %burn = min(1000,cnt-50);
                ens_rat =[ensemble.lograt];
                ens_sig =[ensemble.sig];
                ens_DF =[ensemble.DFgain];
                ens_BL =[ensemble.BL];
                ens_E =[ensemble.E];
                ens_E0 =[ensemble.E0];
                ratmean = mean(ens_rat(:,burn:cnt),2);
                ratstd = std(ens_rat(:,burn:cnt),[],2);


                BLmean = mean(ens_BL(:,burn:cnt),2);
                BLstd = std(ens_BL(:,burn:cnt),[],2);

                sigmean = mean(ens_sig(:,burn:cnt),2);
                sigstd = std(ens_sig(:,burn:cnt),[],2);

                DFmean = mean(ens_DF(:,burn:cnt),2);
                DFstd = std(ens_DF(:,burn:cnt),[],2);

                for m=1:d0.Nblock
                    for n = 1:cnt
                        ens_I{m}(:,n) =[ensemble(n).I{m}];
                    end
                    Imean{m} = mean(ens_I{m}(:,burn:cnt),2);
                    Istd{m} = std(ens_I{m}(:,burn:cnt),[],2);
                end

                %PlotEnsemble_synth;  % Uncomment for plotting
                %PlotEnsemble_MS;
                %PlotData_MS
            end


        end % End saving/displaying/plotting



        % If number of iterations is square number, larger than effective
        % sample size, test for convergence
        if mod(sqrt(cnt),1)==0 && cnt >= EffectSamp/datsav

            cnt2 = cnt2+1;
            Rexit = GRConverge(x,ensemble);  %Gelman-Rubin multivariate criterium

            rrr(cnt2) = Rexit; %debug

            if Rexit<=ExitCrit
                fprintf('MCMC exiting after %d iters with R of %0.6f\n',m,Rexit)
                break
            end


        end



    end  % End of MCMC iterations




    %% Analysis and Plotting

    burn = 1; % Number of models to discard
    ens_rat =[ensemble.lograt];
    ens_sig =[ensemble.sig];
    ens_DF =[ensemble.DFgain];
    ens_BL =[ensemble.BL];
    ens_E =[ensemble.E];
    ens_E0 =[ensemble.E0];


    for m=1:d0.Nblock
        for n = 1:cnt
            ens_I{m}(:,n) =[ensemble(n).I{m}];
        end
        Imean{m} = mean(ens_I{m}(:,burn:cnt),2);  % Mean intensity knots
        Istd{m} = std(ens_I{m}(:,burn:cnt),[],2); % Std dev intensity knots
    end



    % Calculate mean and st dev of ratios after burn in time
    ratmean = mean(ens_rat(:,burn:cnt),2);  % Log ratios
    ratstd = std(ens_rat(:,burn:cnt),[],2);

    BLmean = mean(ens_BL(:,burn:cnt),2);  % Baselines
    BLstd = std(ens_BL(:,burn:cnt),[],2);

    sigmean = mean(ens_sig(:,burn:cnt),2);   % Noise hyperparams
    sigstd = std(ens_sig(:,burn:cnt),[],2);

    DFmean = mean(ens_DF(:,burn:cnt),2);   % Daly-Far gain
    DFstd = std(ens_DF(:,burn:cnt),[],2);




%%

    ensname = sprintf('%s/Ensemble_run%02d.mat',outfolder,iset);
    save(ensname,'ensemble','d0','dataFolder','d00','data','Iso_Name','prior',...
        'psig','maxcnt','datsav','burn','cnt','method','*mean','*std','x','iblock')

     
    PlotEnsemble_MS

    PlotProgress_MS

    PlotData_MS





end


