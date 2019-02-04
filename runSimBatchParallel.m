%just to be sure
%function runSimBatchParallel(strConnFile,strStimFile)
clearvars;
hTic = tic;

%% set log path
global strLogDir;
global boolClust;
boolClust=false;

%add directories to paths
strFile = mfilename;
strHome = mfilename('fullpath');
strHome = strHome(1:(end-length(mfilename)));
strLogDir = [strHome 'logs' filesep];
strOutputDir = 'D:\Data\Processed\V1_LIFmodel\';
strStimDriveDir = [strHome 'StimDrives' filesep];
strConnDir = [strHome 'Connectivity' filesep];
strStimDir = [strHome 'Stimulation' filesep];
strSpec = '';

%% start parallel pool
hPool=gcp('nocreate');
if isempty(hPool)
	hPool = parpool('CPU multicore',14);
end

%% set parameters
sParams = struct;
sParams.strHome = strHome;
sParams.strOutputDir = strOutputDir;
sParams.strStimDriveDir = strStimDriveDir;
sParams.dblDeltaT = 0.5/1000; %seconds
sParams.dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
strStimFile = 'sStim_SquareGrating_x7500R2500_2017-04-20.mat';
intMaxReps = 2;
if isempty(strStimFile)
	sStimParams = struct;
	sStimParams.dblDeltaT = sParams.dblDeltaT;
	sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
	sStimParams.dblStartFirstTrialSecs = 1;
	sStimParams.dblPreStimBlankDur = 0.1;%0.25
	sStimParams.dblStimDur = 0.2;%0.5
	sStimParams.dblPostStimBlankDur = 0.1;%0.25
	sStimParams.vecOrientations = [42.5 47.5];%[0:(180/12):179];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
	sStimParams.vecSpatialFrequencies = 0.5;%2.^[-4:1];
	sStimParams.vecTemporalFrequencies = 0;
	sStimParams.intReps = 2; %number of repetitions
	sStimParams.dblStimSizeRetDeg = 16;
	sStimParams.vecScrPixWidthHeight = [128 128];
	sStimParams.vecScrDegWidthHeight = [25.6 25.6];
	sStimParams.strStimDriveDir = strStimDriveDir;
	
	% build stimuli
	sStimInputs = loadSimStim(sStimParams);
	
	%add
	sParams.sStimParams = sStimParams;
	sParams.sStimInputs = sStimInputs;
else
	[sStimParams,sStimInputs] = loadStimulation(strStimDir,strStimFile);
	sStimParams.intReps = intMaxReps;
	intStimTypes = numel(unique(sStimInputs.vecTrialStimType));
	intNumTrials = intMaxReps*intStimTypes;
	sStimParams.strStimType = 'Uniform';
	
	%remove excess trials
	sStimInputs.vecTrialOris = sStimInputs.vecTrialOris(1:intNumTrials);
	sStimInputs.vecTrialOriIdx = sStimInputs.vecTrialOriIdx(1:intNumTrials);
	sStimInputs.vecTrialSFs = sStimInputs.vecTrialSFs(1:intNumTrials);
	sStimInputs.vecTrialSFIdx = sStimInputs.vecTrialSFIdx(1:intNumTrials);
	sStimInputs.vecTrialTFs = sStimInputs.vecTrialTFs(1:intNumTrials);
	sStimInputs.vecTrialTFIdx = sStimInputs.vecTrialTFIdx(1:intNumTrials);
	sStimInputs.vecTrialContrasts = sStimInputs.vecTrialContrasts(1:intNumTrials);
	sStimInputs.vecTrialContrastIdx = sStimInputs.vecTrialContrastIdx(1:intNumTrials);
	sStimInputs.vecTrialStimType = sStimInputs.vecTrialStimType(1:intNumTrials);
	sStimInputs.vecTrialStimRep = sStimInputs.vecTrialStimRep(1:intNumTrials);
	sStimInputs.vecTrialStartSecs = sStimInputs.vecTrialStartSecs(1:intNumTrials);
	sStimInputs.vecTrialEndSecs = sStimInputs.vecTrialEndSecs(1:intNumTrials);
	sStimInputs.vecStimStartSecs = sStimInputs.vecStimStartSecs(1:intNumTrials);
	sStimInputs.vecStimStopSecs = sStimInputs.vecStimStopSecs(1:intNumTrials);

	%add
	sParams.sStimParams = sStimParams;
	sParams.sStimInputs = sStimInputs;
	
	% msg
	fprintf('Loaded stimulation file <%s> [%s]\n',strStimFile,getTime);
end
strSpec = 'xAreaTransferFeedback3';

%% load connectivity
if ~exist('strConnFile','var')
	%strConnFile = 'sConn_Col48N2160S637056_2017-03-31.mat';
	strConnFile = 'sConn_Col48N2160S637056_2017-04-10.mat';
	%strConnFile = 'sConn_Ret1N1440S403200_2017-04-10.mat';
	%strConnFile = 'sConn_Col48N2160S637056_2017-04-11.mat';
	%strConnFile = 'sConnInhFB_Col48N2160S671616_2017-04-13.mat';
	%strConnFile = 'sConn_Col32N1920S481920_2017-04-13.mat';
end
[sConnParams,sConnectivity] = loadConnectivity_xArea(strConnDir,strConnFile);

%assign
sParams.sConnParams = sConnParams;
sParams.sConnectivity = sConnectivity;
sParams.strConnFile = strConnFile;

%% msg
printf('Preparations complete; will be using connectivity matrix <%s> [%s]\n',strConnFile,getTime);

%% set variables for parallel runs
intSplitType = 1; %split by repetitions (1) or by trials (2)
if intSplitType == 1
	%get data
	vecRepetitions = unique(sParams.sStimInputs.vecTrialStimRep);
	indTrialNewRep = logical([1 diff(sParams.sStimInputs.vecTrialStimRep)]);
	intCortexCells = sConnectivity.intCellsV2 + sConnectivity.intCellsV1;
	
	%get starts/stops of repetition runs
	vecStartRuns = sParams.sStimInputs.vecTrialStartSecs(indTrialNewRep);
	vecStopRuns = [vecStartRuns(2:end) sParams.sStimInputs.vecTrialEndSecs(end)]-sParams.dblDeltaT;
	
	%build matrix
	matOverallT=[];
	for intRun=1:numel(vecStartRuns)
		matOverallT = [matOverallT; vecStartRuns(intRun):sParams.dblDeltaT:vecStopRuns(intRun)];
	end
	
	%assign start trials
	vecPrevTrial = find(indTrialNewRep)-1;
	vecTrialT = zeros(size(vecPrevTrial));
elseif intSplitType == 2
	%get data
	vecRepetitions = unique(sParams.sStimInputs.vecTrialStimRep);
	indTrialNewRep = logical([1 diff(sParams.sStimInputs.vecTrialStimRep)]);
	
	%get starts/stops of repetition runs
	vecStartRuns = sParams.sStimInputs.vecTrialStartSecs;
	vecStopRuns = [vecStartRuns(2:end) sParams.sStimInputs.vecTrialEndSecs(end)]-sParams.dblDeltaT;
	
	%build matrix
	matOverallT=[];
	for intRun=1:numel(vecStartRuns)
		matOverallT = [matOverallT; vecStartRuns(intRun):sParams.dblDeltaT:vecStopRuns(intRun)];
	end
	
	%assign start trials
	vecPrevTrial = 0:(numel(vecStartRuns)-1);
	vecTrialT = zeros(size(vecPrevTrial));
end

%put in structure
sParallel = struct;
sParallel.vecRepetitions = vecRepetitions;
sParallel.indTrialNewRep = indTrialNewRep;
sParallel.vecStartRuns = vecStartRuns;
sParallel.vecStopRuns = vecStopRuns;
sParallel.matOverallT = matOverallT;
sParallel.vecPrevTrial = vecPrevTrial;
sParallel.vecTrialT = vecTrialT;

%% run simulation
intRuns = numel(vecRepetitions);
sSimRun = struct;
%try
parfor intRun=1:size(matOverallT,1)
	printf('Processing repetition %d/%d [%s]\n',intRun,intRuns,getTime);
	vecOverallT = matOverallT(intRun,:);
	intPrevTrial = vecPrevTrial(intRun);
	intTrialT = vecTrialT(intRun);
	
	%run sim
	sData = getSimLIF_xArea(sParams,vecOverallT,intPrevTrial,intTrialT);
	
	%save
	sSimRun(intRun).cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	sSimRun(intRun).cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	sSimRun(intRun).cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	sSimRun(intRun).vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	sSimRun(intRun).vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	sSimRun(intRun).vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
end

%save data
strDataFile = [strOutputDir 'Simulation_PP'  strSpec getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('_%s_%s.mat',sprintf('%f',now),getDate)];
printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
save(strDataFile,'sSimRun','sParallel','sParams');
printf(' .. Data saved! [%s]\n',getTime);

%catch ME
%	printf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%	fprintf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%end
