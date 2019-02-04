%function runBuildConnectivityFromFile

%% set target file
clearvars;
strSimulation = 'SquareGrating2017-01-25';
load(['D:\Data\Processed\V1_LIFmodel\Simulation_' strSimulation '.mat']);

%% retrieve fields
cellFields=fieldnames(sData);
for intField = 1:numel(cellFields)
	strField = cellFields{intField};
	eval([strField ' = sData.' cellFields{intField} ';']);
end

%% get defaults
if ~exist('intCellsV2','var'),intCellsV2=0;end
if ~exist('intCellsV1','var'),intCellsV1=intCortexCells;end
if exist('matCortConn','var'),
	matSynFromTo=matCortConn;
	vecSynExcInh=vecCortSynType;
	vecSynDelay = vecCortDelay;
	vecSynWeight = vecCortDelay;
	vecSynConductance = vecCortConductance;
	vecSynType = ones(size(vecSynExcInh));
	
	vecCellTauPeak = vecTauPeakByType(vecCellTypes);
	vecCellArea = ones(size(vecCellTypes));
	
	dblSigmaX = 0.7; %length of gabor response
	dblSigmaY = 0.47; %width of gabor response
	matPrefGabors = [];
	matFieldsV2 = [];
end

%% combine data
sConnParams = struct;
sConnectivity = struct;

%summary variables
sConnectivity.intCellsV2 = intCellsV2;
sConnectivity.intCellsV1 = intCellsV1;
sConnectivity.intCortexCells = intCellsV1+intCellsV2;

%LGN->cortex connectivity
sConnectivity.vecSynConductanceON_to_Cort = vecSynConductanceON_to_Cort;
sConnectivity.vecSynConductanceOFF_to_Cort = vecSynConductanceOFF_to_Cort;
sConnectivity.vecSynWeightON_to_Cort = vecSynWeightON_to_Cort;
sConnectivity.vecSynWeightOFF_to_Cort = vecSynWeightOFF_to_Cort;
sConnectivity.vecSynDelayON_to_Cort = vecSynDelayON_to_Cort;
sConnectivity.vecSynDelayOFF_to_Cort = vecSynDelayOFF_to_Cort;
sConnectivity.matSynConnON_to_Cort = matSynConnON_to_Cort;
sConnectivity.matSynConnOFF_to_Cort = matSynConnOFF_to_Cort;

%cortical connectivity
sConnectivity.matSynFromTo = matSynFromTo;
sConnectivity.vecSynExcInh = vecSynExcInh;
sConnectivity.vecSynDelay = vecSynDelay;
sConnectivity.vecSynWeight = vecSynWeight;
sConnectivity.vecSynConductance = vecSynConductance;
sConnectivity.vecSynType = vecSynType;

%cortical cell parameters
sConnectivity.vecCellThresh = vecCellThresh;
sConnectivity.vecTauPeakByType = vecTauPeakByType;
sConnectivity.vecCellTauPeak = vecCellTauPeak;
sConnectivity.vecCellV_E = vecCellV_E;
sConnectivity.vecCellV_I = vecCellV_I;
sConnectivity.vecCellV_AHP = vecCellV_AHP;
sConnectivity.vecCellV_Leak = vecCellV_Leak;
sConnectivity.vecCellCm = vecCellCm;
sConnectivity.vecCellG_Leak = vecCellG_Leak;
sConnectivity.vecCellG_AHP = vecCellG_AHP;

%cell preferences
sConnectivity.vecCellArea = vecCellArea;
sConnectivity.vecCellTypes = vecCellTypes;
sConnectivity.vecPrefPsi = vecPrefPsi; %cell's phase offset
sConnectivity.vecPrefOri=vecPrefOri; %cell's preferred orientation
sConnectivity.vecPrefSF = vecPrefSF;
sConnectivity.vecPrefRF_X = vecPrefRF_X;
sConnectivity.vecPrefRF_Y = vecPrefRF_Y;
sConnectivity.dblSigmaX = dblSigmaX; %length of gabor response
sConnectivity.dblSigmaY = dblSigmaY; %width of gabor response
sConnectivity.matPrefGabors = matPrefGabors;
sConnectivity.matFieldsV2 = matFieldsV2;


%% save
strConnFile = sprintf('sConn_Ret1N%dS%d_%s.mat',sConnectivity.intCortexCells,numel(sConnectivity.vecSynExcInh),getDate);
strConnDir = 'D:\Simulations\Connectivity\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strConnFile,strConnDir,getTime);
save([strConnDir strConnFile],'sConnParams','sConnectivity');

%% stimulus
%stim inputs
sStimInputs.dblTrialDur = sStimInputs.dblTrialDur;
sStimInputs.dblVisSpacing = sStimInputs.dblVisSpacing;
sStimInputs.varDeltaSyn = sStimInputs.varDeltaSyn;
sStimInputs.matBlankLGN_ON = sStimInputs.matBlankLGN_ON;
sStimInputs.matBlankLGN_OFF = sStimInputs.matBlankLGN_OFF;

sStimInputs.cellLGN_ON = sStimInputs.cellLGN_ON;
sStimInputs.cellLGN_OFF = sStimInputs.cellLGN_OFF;
sStimInputs.cellContrast = cellfill(6,size(sStimInputs.cellLGN_ON));
sStimInputs.vecTrialOris = sStimInputs.vecTrialOris;
sStimInputs.vecTrialOriIdx = sStimInputs.vecTrialOriIdx;
sStimInputs.vecTrialSFs = sStimInputs.vecTrialSFs;
sStimInputs.vecTrialSFIdx = sStimInputs.vecTrialSFIdx;
sStimInputs.vecTrialTFs = sStimInputs.vecTrialTFs;
sStimInputs.vecTrialTFIdx = sStimInputs.vecTrialTFIdx;
sStimInputs.vecTrialStimType = sStimInputs.vecTrialStimType;
sStimInputs.vecTrialStimRep = sort(repmat(1:sStimParams.intReps,[1 numel(unique(sStimInputs.vecTrialStimType))]),'ascend');
sStimInputs.vecTrialStartSecs = sStimInputs.vecTrialEndSecs - sStimInputs.dblTrialDur; %% !!!!!!! sStimInputs.vecTrialStartSecs = sData.vecTrialStartSecs;
sStimInputs.vecTrialEndSecs = sStimInputs.vecTrialEndSecs;
sStimInputs.vecStimStartSecs = sStimInputs.vecStimStartSecs;
sStimInputs.vecStimStopSecs = sStimInputs.vecStimStartSecs+sData.dblStimDur;
sStimInputs.matStimTypeCombos = sStimInputs.matStimTypeCombos;
sStimInputs.vecStimTypeOris = sStimInputs.vecStimTypeOris;
sStimInputs.vecStimTypeSFs = sStimInputs.vecStimTypeSFs;
sStimInputs.vecStimTypeTFs = sStimInputs.vecStimTypeTFs;

%% stim params
sStimParams.dblDeltaT = sStimParams.dblDeltaT;
sStimParams.strStimType = sStimParams.strStimType;
sStimParams.dblStartFirstTrialSecs = sStimParams.dblStartFirstTrialSecs;
sStimParams.dblPreStimBlankDur = sStimParams.dblPreStimBlankDur ;
sStimParams.dblStimDur = sStimParams.dblStimDur;
sStimParams.dblPostStimBlankDur = sStimParams.dblPostStimBlankDur;
sStimParams.vecOrientations = sStimParams.vecOrientations;
sStimParams.vecSpatialFrequencies = sStimParams.vecSpatialFrequencies;
sStimParams.vecTemporalFrequencies = sStimParams.vecTemporalFrequencies;
sStimParams.intReps = sStimParams.intReps;
sStimParams.dblStimSizeRetDeg = sStimParams.dblStimSizeRetDeg;
sStimParams.vecScrPixWidthHeight = sStimParams.vecScrPixWidthHeight;
sStimParams.vecScrDegWidthHeight = sStimParams.vecScrDegWidthHeight;

%% save
strStimFile = sprintf('sStim_Ret1FixedTrials%dTypes%dReps%d_%s.mat',numel(sStimInputs.vecStimStartSecs),numel(unique(sStimInputs.vecTrialStimType)),sStimParams.intReps,getDate);
strStimDir = 'D:\Simulations\Stimulation\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strStimFile,strStimDir,getTime);
save([strStimDir strStimFile],'sStimInputs','sStimParams');
