%just to be sure
%function runSimBatch
clear all;

%% set log path
global strLogDir;
global boolClust;

%% running on cluster?
boolClust=false;
if boolClust
	
	%add directories to paths
	strFile = mfilename;
	strHome = mfilename('fullpath');
	strHome = strHome(1:(end-length(strFile)));
	addpath([strHome 'Code']);
	strLogDir = [strHome 'logs' filesep];
	strOutputDir = [strHome 'SimResults' filesep];
	strStimDriveDir = [strHome 'StimDrives' filesep];
	%add recursive subdirectories
	cellPaths = getSubDirs(strHome,inf,{'old','backup'});
	printf('path 1: %s; nr of paths=%d\n',cellPaths{1},numel(cellPaths));
	for intPath=1:numel(cellPaths)
		addpath(cellPaths{intPath});
	end
	
	%start parallel pool
	intCPUs = feature('numThreads');
	printf('Number of available CPUs: %d\n',intCPUs);
	parpool(intCPUs);
else
	%add directories to paths
	strFile = mfilename;
	strHome = mfilename('fullpath');
	strHome = strHome(1:(end-length(mfilename)));
	strLogDir = [strHome 'logs' filesep];
	strOutputDir = 'D:\Data\Processed\V1_LIFmodel\';
	strStimDriveDir = [strHome 'StimDrives' filesep];
end

%% set parameters
sParams = struct;
sParams.strHome = strHome;
sParams.strOutputDir = strOutputDir;
sParams.strStimDriveDir = strStimDriveDir;
sParams.dblDeltaT = 0.5/1000; %seconds
sParams.dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
sStimParams = struct;
sStimParams.dblDeltaT = sParams.dblDeltaT;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0.1;
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 0.3;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [42.5 47.5];%[0:(180/12):179];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
sStimParams.vecSpatialFrequencies = 0.5;%2.^[-4:1];
sStimParams.vecTemporalFrequencies = 3;
sStimParams.intReps = 20; %number of repetitions
sStimParams.dblStimSizeRetDeg = 16;
sStimParams.vecScrPixWidthHeight = [128 128];
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
sStimParams.strStimDriveDir = strStimDriveDir;
sParams.sStimParams = sStimParams;

%% build stimuli
sParams.sStimInputs = loadSimStim(sStimParams);

%% set connectivity parameters, or load pre-existing connectivity structure
intColumns = 48; %48 / 252
sConnParams = struct;
sConnParams.dblVisSpacing = 0.2;
sConnParams.vecSizeInput = [64 64];

%connection definition LGN
sConnParams.vecConnsPerTypeON = [24 16]; %[pyramid interneuron]
sConnParams.vecConnsPerTypeOFF = [24 16]; %[pyramid interneuron]

sConnParams.dblSigmaX = 0.7; %length of gabor response
sConnParams.dblSigmaY = 0.47; %width of gabor response
sConnParams.vecConductance_FromLGN_ToCort = [5.5 6]*0.3; %to [pyramid interneuron]
sConnParams.vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
sConnParams.vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

%V1 def
sConnParams.vecDefinitionV1PrefOri = 0:pi/intColumns:(pi-pi/intColumns);
sConnParams.vecDefinitionV1SpatFreq = 2.^[-4:1:1];
sConnParams.vecDefinitionV1CellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
sConnParams.intCellsV1 = numel(sConnParams.vecDefinitionV1PrefOri) * numel(sConnParams.vecDefinitionV1SpatFreq) * numel(sConnParams.vecDefinitionV1CellTypes);
	
%cortical connectivity
%number of connections
dblScalingFactor = 4; %4
sConnParams.matConnCortFromTo(1,:) = [40 40]*dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConnCortFromTo(2,:) = [30 30]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
sConnParams.matConductancesFromTo(1,:) = [1.1 1.6]/dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromTo(2,:) = [1.5 1.0]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanCortToCort = 3/1000; %in ms
sConnParams.dblDelaySDCortToCort = 1/1000; %in ms

%connection probability excitatory cells
sConnParams.vecConnProbSD = ang2rad([7.5 80]); %degrees difference in pref ori, for [pyramid interneuron]

%V2 params
sConnParams.dblSpatialDropoffV1V2 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
sConnParams.dblSpatialDropoffInterneuronsV2 = 3; %for interneurons
sConnParams.intCellsV2 = round(sConnParams.intCellsV1/2);

%create cell-based parameters
sConnParams.vecDefinitionV2CellTypes = [1 1 1 1 2];
sConnParams.vecCellTypesV2 = repmat(sConnParams.vecDefinitionV2CellTypes,[1 ceil(sConnParams.intCellsV2/numel(sConnParams.vecDefinitionV2CellTypes))]);
sConnParams.vecCellTypesV2(sConnParams.intCellsV2+1:end) = [];
sConnParams.vecCellFractionsV2 = [0.8 0.2]; %[pyramid interneuron]
	
%V1=>V2
dblInterArealFactor = 2.6;
sConnParams.vecConnsPerTypeV1V2 = [48 32];%[48 32]; %pyr/int
sConnParams.matConductancesFromToV1V2(1,:) = [1.1 1.6]*dblInterArealFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromToV1V2(2,:) = [1.5 1.0]*dblInterArealFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanV1ToV2 = 4/1000; %in ms
sConnParams.dblDelaySDV1ToV2 = 1/1000; %in ms

%assign
sParams.sConnParams = sConnParams;

%% build connectivity
sParams.sConnectivity = buildConnectivity_xArea(sConnParams);

%% set variables for parallel runs
%get data
vecRepetitions = unique(sParams.sStimInputs.vecTrialStimRep);
indTrialNewRep = logical([1 diff(sParams.sStimInputs.vecTrialStimRep)]);
intCortexCells = sConnParams.intCellsV2 + sConnParams.intCellsV1;

%get starts/stops of repetition runs
vecStartRuns = sParams.sStimInputs.vecTrialStartSecs(indTrialNewRep);
vecStopRuns = [vecStartRuns(2:end) sParams.sStimInputs.vecTrialEndSecs(end)]-sParams.dblDeltaT;

%build matrix
matOverallT=[];
for intRep=vecRepetitions
	matOverallT = [matOverallT; vecStartRuns(intRep):sParams.dblDeltaT:vecStopRuns(intRep)];
end

%assign start trials
vecPrevTrial = find(indTrialNewRep)-1;
vecTrialT = find(indTrialNewRep)-1;

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
intReps = numel(vecRepetitions);
sSimRun = struct;
%try
parfor intRep=vecRepetitions
	printf('Processing repetition %d/%d [%s]\n',intRep,intReps,getTime);
	vecOverallT = matOverallT(intRep,:);
	intPrevTrial = vecPrevTrial(intRep);
	intTrialT = vecTrialT(intRep);
	
	%run sim
	sData = getSimLIF_xArea(sParams,vecOverallT,intPrevTrial,intTrialT);
	
	%save
	sSimRun(intRep).cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	sSimRun(intRep).cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	sSimRun(intRep).cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	sSimRun(intRep).vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	sSimRun(intRep).vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	sSimRun(intRep).vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
end
boolRemNans = true;
cellFields = {'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF','cellSpikeTimesCortex'};
sSimRun = catfields(sSimRun,cellFields,boolRemNans);

%save data
strDataFile = [strOutputDir 'Simulation_'  getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('_%f_%s.mat',now,getDate)];
printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
save(strDataFile,'sSimRun','sParallel','sParams');
printf(' .. Data saved! [%s]\n',getTime);

%catch ME
%	printf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%	fprintf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);	
%end
