%just to be sure
%function runSimBatcSerial(strConnFile,strStimFile)
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

%% set parameters
sParams = struct;
sParams.strHome = strHome;
sParams.strOutputDir = strOutputDir;
sParams.strStimDriveDir = strStimDriveDir;
sParams.dblDeltaT = 0.5/1000; %seconds
sParams.dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
strStimFile = 'sStim_xAreaTransfer_x251R1_2017-04-12.mat';
intMaxReps = 1;
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
strSpec = 'xAreaTransferFeedback';

%% load connectivity
if ~exist('strConnFile','var')
	%strConnFile = 'sConn_Col48N2160S637056_2017-03-31.mat';
	%strConnFile = 'sConn_Col48N2160S637056_2017-04-10.mat';
	%strConnFile = 'sConn_Ret1N1440S403200_2017-04-10.mat';
	%strConnFile = 'sConn_Col48N2160S637056_2017-04-11.mat';
	strConnFile = 'sConn_Col48N2160S671616_2017-04-12.mat'; %with feedback V2=> V1
end
[sConnParams,sConnectivity] = loadConnectivity_xArea(strConnDir,strConnFile);

%assign
sParams.sConnParams = sConnParams;
sParams.sConnectivity = sConnectivity;
sParams.strConnFile = strConnFile;

%% msg
fprintf('Preparations complete; will be using connectivity matrix <%s> [%s]\n',strConnFile,getTime);

%% prepare variables for simulation
% run stitched together trials
dblDeltaT = sStimParams.dblDeltaT;
dblSimDur = max(sStimInputs.vecTrialEndSecs); %seconds
intCortexCells = sConnectivity.intCortexCells;
vecThisV = ones(intCortexCells,1)*-65;
vecPrevV = vecThisV;
boolStimPresent = false;
intPrevTrial = 0;
intTrialT = 0;
vecOverallT = dblDeltaT:dblDeltaT:dblSimDur;
fprintf('\n   >>> Starting simulation run [%s]\n',getTime);
intIter = 0;

intLGN_CellsPerStream = size(sStimInputs.cellLGN_ON{1},1)*size(sStimInputs.cellLGN_ON{1},2);
cellSpikeTimesLGN_ON = cell(intLGN_CellsPerStream,1);
cellSpikeTimesLGN_OFF = cell(intLGN_CellsPerStream,1);
cellSpikeTimesCortex = cell(intCortexCells,1);

vecSpikeCounterLGN_ON = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterLGN_OFF = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterCortex = zeros(intCortexCells,1);

vecSpikeCounterPreAllocatedLGN_ON = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterPreAllocatedLGN_OFF = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterPreAllocatedCortex = zeros(intCortexCells,1);

intPreAllocationSize = 1000;

%% assign to structure
vecTauPeakByType = [1/1000 2/1000];
sData.vecOverallT = vecOverallT;
sData.dblDeltaT = dblDeltaT;
sData.intCellsV1 = sConnectivity.intCellsV1;
sData.intCellsV2 = sConnectivity.intCellsV2;

sData.matSynFromTo = sConnectivity.matSynFromTo;
sData.dblSynSpikeMem = 0.200;
sData.vecSynExcInh = sConnectivity.vecSynExcInh;
sData.intCortexCells = sConnectivity.intCortexCells;
sData.vecSynDelay = sConnectivity.vecSynDelay;
sData.vecSynConductance = sConnectivity.vecSynConductance;
sData.vecSynWeight = sConnectivity.vecSynWeight;
sData.vecSynType = sConnectivity.vecSynType;

sData.vecCellThresh = sConnectivity.vecCellThresh;
sData.vecTauPeakByType = sConnectivity.vecTauPeakByType;
sData.vecCellV_E = sConnectivity.vecCellV_E;
sData.vecCellV_I = sConnectivity.vecCellV_I;
sData.vecCellV_AHP = sConnectivity.vecCellV_AHP;
sData.vecCellV_Leak = sConnectivity.vecCellV_Leak;
sData.vecCellCm = sConnectivity.vecCellCm;
sData.vecCellG_Leak = sConnectivity.vecCellG_Leak;
sData.vecCellG_AHP = sConnectivity.vecCellG_AHP;
sData.vecSynConductanceON_to_Cort = sConnectivity.vecSynConductanceON_to_Cort;
sData.vecSynConductanceOFF_to_Cort = sConnectivity.vecSynConductanceOFF_to_Cort;
sData.vecSynWeightON_to_Cort = sConnectivity.vecSynWeightON_to_Cort;
sData.vecSynWeightOFF_to_Cort = sConnectivity.vecSynWeightOFF_to_Cort;
sData.vecSynDelayON_to_Cort = sConnectivity.vecSynDelayON_to_Cort;
sData.vecSynDelayOFF_to_Cort = sConnectivity.vecSynDelayOFF_to_Cort;
sData.matSynConnON_to_Cort = sConnectivity.matSynConnON_to_Cort;
sData.matSynConnOFF_to_Cort = sConnectivity.matSynConnOFF_to_Cort;

sData.matBlankLGN_ON = sStimInputs.matBlankLGN_ON;
sData.matBlankLGN_OFF = sStimInputs.matBlankLGN_OFF;
sData.cellLGN_ON = sStimInputs.cellLGN_ON;
sData.cellLGN_OFF = sStimInputs.cellLGN_OFF;
sData.vecTrialOris = sStimInputs.vecTrialOris;
sData.vecTrialOriIdx = sStimInputs.vecTrialOriIdx;
sData.vecTrialStimType = sStimInputs.vecTrialStimType;
sData.vecTrialStartSecs = sStimInputs.vecTrialStartSecs;
sData.vecStimStartSecs = sStimInputs.vecStimStartSecs;
sData.vecStimStopSecs = sStimInputs.vecStimStopSecs;
sData.vecTrialEndSecs = sStimInputs.vecTrialEndSecs;
sData.vecThisV = vecThisV;
sData.boolStimPresent = double(boolStimPresent);
sData.intPrevTrial = intPrevTrial;
sData.intTrialT = intTrialT;
sData.intIter = intIter;

sData.cellSpikeTimesLGN_ON = cellSpikeTimesLGN_ON;
sData.cellSpikeTimesLGN_OFF = cellSpikeTimesLGN_OFF;
sData.cellSpikeTimesCortex = cellSpikeTimesCortex;
sData.vecSpikeCounterLGN_ON = vecSpikeCounterLGN_ON;
sData.vecSpikeCounterLGN_OFF = vecSpikeCounterLGN_OFF;
sData.vecSpikeCounterCortex = vecSpikeCounterCortex;
sData.intPreAllocationSize = intPreAllocationSize;
sData.dblStimDur = sStimParams.dblStimDur;

%cell preferences
sData.vecCellTypes = sConnectivity.vecCellTypes;
sData.vecPrefPsi = sConnectivity.vecPrefPsi; %cell's phase offset
sData.vecPrefOri = sConnectivity.vecPrefOri; %cell's preferred orientation
sData.vecPrefSF = sConnectivity.vecPrefSF;
sData.vecPrefRF_X = sConnectivity.vecPrefRF_X;
sData.vecPrefRF_Y = sConnectivity.vecPrefRF_Y;
sData.dblSigmaX = sConnectivity.dblSigmaX; %length of gabor response
sData.dblSigmaY = sConnectivity.dblSigmaY; %width of gabor response
sData.matPrefGabors = sConnectivity.matPrefGabors;

%optional fields
if isfield(sStimInputs,'vecTrialInputStrength'),sData.vecTrialInputStrength=sStimInputs.vecTrialInputStrength;end
if isfield(sStimInputs,'vecTrialSpikeProbLGN'),sData.vecTrialSpikeProbLGN=sStimInputs.vecTrialSpikeProbLGN;end

%% run simulation
sData = getSimulationRun_xArea(sData);

%% remove placeholder spike entries
cellSpikeTimesTemp = sData.cellSpikeTimesCortex;
for intN=1:numel(cellSpikeTimesTemp)
	if numel(cellSpikeTimesTemp{intN}) > 1
		cellSpikeTimesTemp{intN} = cellSpikeTimesTemp{intN}(~isnan(cellSpikeTimesTemp{intN}));
	else
		cellSpikeTimesTemp{intN} = [];
	end
end
sData.cellSpikeTimesCortex = cellSpikeTimesTemp;

cellSpikeTimesTemp = sData.cellSpikeTimesLGN_ON;
for intN=1:numel(cellSpikeTimesTemp)
	if numel(cellSpikeTimesTemp{intN}) > 1
		cellSpikeTimesTemp{intN} = cellSpikeTimesTemp{intN}(~isnan(cellSpikeTimesTemp{intN}));
	else
		cellSpikeTimesTemp{intN} = [];
	end
end
sData.cellSpikeTimesLGN_ON = cellSpikeTimesTemp;

cellSpikeTimesTemp = sData.cellSpikeTimesLGN_OFF;
for intN=1:numel(cellSpikeTimesTemp)
	if numel(cellSpikeTimesTemp{intN}) > 1
		cellSpikeTimesTemp{intN} = cellSpikeTimesTemp{intN}(~isnan(cellSpikeTimesTemp{intN}));
	else
		cellSpikeTimesTemp{intN} = [];
	end
end
sData.cellSpikeTimesLGN_OFF = cellSpikeTimesTemp;

%% save data
strDataFile = [strOutputDir 'Simulation_'  strSpec getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('_%s_%s.mat',sprintf('%f',now),getDate)];
fprintf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
save(strDataFile,'sData','sParams');
fprintf(' .. Data saved! [%s]\n',getTime);

%catch ME
%	printf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%	fprintf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%end
