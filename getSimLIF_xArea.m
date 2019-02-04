function sData = getSimLIF_xArea(sParams,vecOverallT,intPrevTrial,intTrialT)

%% set parameters
dblDeltaT = sParams.dblDeltaT; %seconds
dblSynSpikeMem = sParams.dblSynSpikeMem; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% check supplied paths
if nargin == 0 || isempty(sParams)
	strHome = 'D:\Simulations';
	strOutputDir = 'D:\Data\Processed\V1_LIFmodel\';
	strStimDriveDir = [strHome 'StimDrives' filesep];
	boolVerbose = true;
else
	strHome = sParams.strHome;
	strOutputDir = sParams.strOutputDir;
	strStimDriveDir = sParams.strStimDriveDir;
	boolVerbose = false;
end

%% prepare stimulus list
if ~isfield(sParams,'sStimParams')
	sStimParams = struct;
	sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
	sStimParams.dblStartFirstTrialSecs = 0.1;
	sStimParams.dblPreStimBlankDur = 0.1;%0.25
	sStimParams.dblStimDur = 0.3;%0.5
	sStimParams.dblPostStimBlankDur = 0.1;%0.25
	sStimParams.vecOrientations = [42.5 47.5];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
	sStimParams.vecSpatialFrequencies = 0.5;%2.^[-4:1];
	sStimParams.vecTemporalFrequencies = 3;
	sStimParams.intReps = 2; %number of repetitions
	sStimParams.dblStimSizeRetDeg = 16;
	sStimParams.vecScrPixWidthHeight = [128 128];
	sStimParams.vecScrDegWidthHeight = [25.6 25.6];
else
	sStimParams = sParams.sStimParams;
end
sStimParams.dblDeltaT = dblDeltaT;
sStimParams.strStimDriveDir = strStimDriveDir;

%get stimuli
if ~isfield(sParams,'sStimInputs')
	sStimInputs = loadSimStim(sStimParams);
else
	sStimInputs = sParams.sStimInputs;
end

if nargin == 1
	figure
	colormap(grey)
	matPlotOFF = sStimInputs.cellLGN_OFF{end};
	matPlotON = sStimInputs.cellLGN_ON{end};
	matPlot = matPlotON - matPlotOFF;
	dblMin = min(matPlot(:));
	dblMax = max(matPlot(:));
	for i=1:10:size(matPlot,3)
		imagesc(matPlot(:,:,i),[dblMin dblMax]);
		title(sprintf('%d',i));
		drawnow
	end
	close;
end

%extract output
dblStimDur = sStimParams.dblStimDur;
dblTrialDur = sStimInputs.dblTrialDur;
varDeltaSyn = sStimInputs.varDeltaSyn;
dblVisSpacing = sStimInputs.dblVisSpacing;
matBlankLGN_ON = sStimInputs.matBlankLGN_ON;
matBlankLGN_OFF = sStimInputs.matBlankLGN_OFF;
cellLGN_ON = sStimInputs.cellLGN_ON;
cellLGN_OFF = sStimInputs.cellLGN_OFF;
vecTrialOris = sStimInputs.vecTrialOris;
vecTrialOriIdx = sStimInputs.vecTrialOriIdx;
vecTrialSFs = sStimInputs.vecTrialSFs;
vecTrialSFIdx = sStimInputs.vecTrialSFIdx;
vecTrialStimType = sStimInputs.vecTrialStimType;
vecTrialStimRep = sStimInputs.vecTrialStimRep;
vecStimStartSecs = sStimInputs.vecStimStartSecs;
vecStimStopSecs = sStimInputs.vecStimStartSecs+sStimParams.dblStimDur;
vecTrialStartSecs = vecStimStartSecs-sStimParams.dblPreStimBlankDur;
vecTrialEndSecs = sStimInputs.vecTrialEndSecs;
vecStimTypeSFs = sStimInputs.vecStimTypeSFs;
vecStimTypeOris = sStimInputs.vecStimTypeOris;

dblSimDur = max(vecTrialEndSecs); %seconds
if boolVerbose,printf(' .. Parameters: time step %.1f ms for %.3f seconds\n',dblDeltaT*1000,dblSimDur);end

%% build connectivity
if ~isfield(sParams,'sConnectivity')
	sData = buildConnectivity_xArea(sParams.sConnParams);
else
	sData = sParams.sConnectivity;
end

%% prepare variables for simulation
% run stitched together trials
intCortexCells = sData.intCortexCells;
vecThisV = gaussrnd(-58,1,[intCortexCells,1]);
boolStimPresent = false;
if nargin == 1
	intPrevTrial = 0;
	intTrialT = 0;
	vecOverallT = dblDeltaT:dblDeltaT:dblSimDur;
end
if boolVerbose,printf('\n   >>> Starting simulation run [%s]\n',getTime);end
intIter = 0;

intLGN_CellsPerStream = size(cellLGN_ON{1},1)*size(cellLGN_ON{1},2);
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
%simulation run parameters
sData.vecOverallT = vecOverallT;
sData.dblDeltaT = dblDeltaT;
sData.dblSynSpikeMem = dblSynSpikeMem;

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

%stim params
sData.matBlankLGN_ON = matBlankLGN_ON;
sData.matBlankLGN_OFF = matBlankLGN_OFF;
sData.cellLGN_ON = cellLGN_ON;
sData.cellLGN_OFF = cellLGN_OFF;

sData.sStimParams = sStimParams;
sData.dblStimDur = dblStimDur;

sData.dblTrialDur = dblTrialDur;
sData.varDeltaSyn = varDeltaSyn;
sData.dblVisSpacing = dblVisSpacing;

sData.vecTrialOris = vecTrialOris;
sData.vecTrialOriIdx = vecTrialOriIdx;
sData.vecStimTypeOris = vecStimTypeOris;

sData.vecTrialSFs = vecTrialSFs;
sData.vecTrialSFIdx = vecTrialSFIdx;
sData.vecStimTypeSFs = vecStimTypeSFs;

sData.vecTrialStimRep = vecTrialStimRep;
sData.vecTrialStimType = vecTrialStimType;
sData.vecTrialStartSecs = vecTrialStartSecs;
sData.vecStimStartSecs = vecStimStartSecs;
sData.vecStimStopSecs = vecStimStopSecs;
sData.vecTrialEndSecs = vecTrialEndSecs;


%% implicitly included by connectivity structure
%{
%overall variables
sData.intCellsV2 = intCellsV2;
sData.intCellsV1 = intCellsV1;
sData.intCortexCells = intCortexCells;

%LGN->cortex connectivity
sData.vecSynConductanceON_to_Cort = vecSynConductanceON_to_Cort;
sData.vecSynConductanceOFF_to_Cort = vecSynConductanceOFF_to_Cort;
sData.vecSynWeightON_to_Cort = vecSynWeightON_to_Cort;
sData.vecSynWeightOFF_to_Cort = vecSynWeightOFF_to_Cort;
sData.vecSynDelayON_to_Cort = vecSynDelayON_to_Cort;
sData.vecSynDelayOFF_to_Cort = vecSynDelayOFF_to_Cort;
sData.matSynConnON_to_Cort = matSynConnON_to_Cort;
sData.matSynConnOFF_to_Cort = matSynConnOFF_to_Cort;

%cortical connectivity
sData.matSynFromTo = matSynFromTo;
sData.vecSynExcInh = vecSynExcInh;
sData.vecSynDelay = vecSynDelay;
sData.vecSynWeight = vecSynWeight;
sData.vecSynConductance = vecSynConductance;
sData.vecSynType = vecSynType;

%cortical cell parameters
sData.vecCellThresh = vecCellThresh;
sData.vecTauPeakByType = vecTauPeakByType;
sData.vecCellV_E = vecCellV_E;
sData.vecCellV_I = vecCellV_I;
sData.vecCellV_AHP = vecCellV_AHP;
sData.vecCellV_Leak = vecCellV_Leak;
sData.vecCellCm = vecCellCm;
sData.vecCellG_Leak = vecCellG_Leak;
sData.vecCellG_AHP = vecCellG_AHP;

%cell preferences
sData.vecCellArea = vecCellArea;
sData.vecCellTypes = vecCellTypes;
sData.vecPrefPsi = vecPrefPsi; %cell's phase offset
sData.vecPrefOri=vecPrefOri; %cell's preferred orientation
sData.vecPrefSF = vecPrefSF;
sData.vecPrefRF_X = vecPrefRF_X;
sData.vecPrefRF_Y = vecPrefRF_Y;
sData.dblSigmaX = dblSigmaX; %length of gabor response
sData.dblSigmaY = dblSigmaY; %width of gabor response
sData.matPrefGabors = matPrefGabors;
sData.matFieldsV2 = matFieldsV2;
%}

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
if nargout == 0
	strDataFile = [strOutputDir 'Simulation_'  getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('%s.mat',getDate)];
	printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
	save(strDataFile,'sData','sStimInputs','sStimParams');
	printf(' .. Data saved! [%s]\n',getTime);
end

%end
%{
%plot spike times
figure
hold on
for intN=1:intCortexCells
	if numel(cellSpikeTimesCortex{intN}) > 0
		line([cellSpikeTimesCortex{intN};cellSpikeTimesCortex{intN}],repmat([intN-0.5;intN+0.5],[1 numel(cellSpikeTimesCortex{intN})]),'Color','k')
	end
end
%}
