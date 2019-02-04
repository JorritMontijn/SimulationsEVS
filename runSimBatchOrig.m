%function runSimLIF
%% starting message
clearvars;
clc;
fprintf(' \nStarting Leaky Integrate and Fire simulation of Retina, LGN and V1 [%s]\n\n',getTime);
%try objParPool = parpool('CPU multicore',10);end%delete(gcp('nocreate'))

%% set parameters
dblDeltaT = 0.5/1000; %seconds
dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
sStimParams = struct;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0.1;
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 0.2;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [42.5 47.5];%[0:(360/16):359]; [42.5 47.5];
sStimParams.vecSpatialFrequencies = 0.5;
sStimParams.vecTemporalFrequencies = 0;
sStimParams.intReps = 1500; %number of repetitions
sStimParams.dblDeltaT = dblDeltaT;
sStimParams.dblStimSizeRetDeg = 16;
sStimParams.vecScrPixWidthHeight = [128 128];
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
	
%get stimuli
sStimInputs = loadSimStim(sStimParams);

figure
colormap(grey)
dblMin = min(sStimInputs.cellLGN_ON{end}(:));
dblMax = max(sStimInputs.cellLGN_ON{end}(:));
for i=1:10:size(sStimInputs.cellLGN_ON{end},3)
	imagesc(sStimInputs.cellLGN_ON{end}(:,:,i),[dblMin dblMax]);
	title(sprintf('%d',i));
	drawnow
end
close;

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
vecTrialStartSecs = sStimInputs.vecTrialStartSecs;
vecStimStartSecs = sStimInputs.vecStimStartSecs;
vecTrialEndSecs = sStimInputs.vecTrialEndSecs;
vecStimTypeSFs = sStimInputs.vecStimTypeSFs;
vecStimTypeOris = sStimInputs.vecStimTypeOris;

dblSimDur = max(vecTrialEndSecs); %seconds
fprintf(' .. Parameters: time step %.1f ms for %.3f seconds\n',dblDeltaT*1000,dblSimDur);

%% get LGN to cortex parameters
%get spatial receptive fields
intImX = size(varDeltaSyn,1);
intImY = size(varDeltaSyn,2);

vecSpaceX = dblVisSpacing*((-(intImX - 1)/2):(intImX - 1)/2);
vecSpaceY = dblVisSpacing*((-(intImY - 1)/2):(intImY - 1)/2);
[matMeshX,matMeshY] = meshgrid(vecSpaceX,vecSpaceY);

%parameters for all cells
intColumns = 48; %252
vecDefinitionPrefOri = 0:pi/intColumns:(pi-pi/intColumns);
vecDefinitionSpatFreq = 2.^[-4:1:1];
vecDefinitionCellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
intCells = numel(vecDefinitionPrefOri) * numel(vecDefinitionSpatFreq) * numel(vecDefinitionCellTypes);

%distribute RF center locations
intCellsPerSpatialDim = ceil(sqrt(intCells));
vecPrefDistroX = linspace(min(vecSpaceX)/2,max(vecSpaceX)/2,intCellsPerSpatialDim);
vecPrefDistroY = linspace(min(vecSpaceY)/2,max(vecSpaceY)/2,intCellsPerSpatialDim);
[matPrefRF_X,matPrefRF_Y] = meshgrid(vecPrefDistroX,vecPrefDistroY);
vecRandRF = randperm(intCellsPerSpatialDim.^2,intCells);
dblSpatialStep = min(diff(vecPrefDistroX));
vecPrefRF_X = matPrefRF_X(vecRandRF)+(dblSpatialStep)*(rand(1,intCells)-0.5);
vecPrefRF_Y = matPrefRF_Y(vecRandRF)+(dblSpatialStep)*(rand(1,intCells)-0.5);
matPrefGabors = nan(numel(vecSpaceX),numel(vecSpaceY),intCells);

%create cell-based parameters
vecPrefOri = nan(1,intCells);
vecPrefSF = nan(1,intCells);
vecCellTypes = nan(1,intCells);
intCounter=0;
intIncrement = numel(vecDefinitionCellTypes);
for dblOri=vecDefinitionPrefOri
	for dblSF = vecDefinitionSpatFreq
		vecPrefOri((intCounter+1):(intCounter+intIncrement)) = dblOri;
		vecPrefSF((intCounter+1):(intCounter+intIncrement)) = dblSF;
		vecCellTypes((intCounter+1):(intCounter+intIncrement)) = vecDefinitionCellTypes;
		intCounter = intCounter + intIncrement;
	end
end

%connection definition LGN
vecCellTypeLabels = unique(vecCellTypes);
vecConnsPerTypeON = [24 16]; %[pyramid interneuron]
vecConnsPerTypeOFF = [24 16]; %[pyramid interneuron]

dblSigmaX = 0.7;
dblSigmaY = 0.47;
intCortexCells = length(vecCellTypes);
vecPrefPsi = 2*pi*rand(1,intCortexCells);%phase offset
intConnsLGN_to_CortOFF = sum(vecConnsPerTypeOFF(vecDefinitionCellTypes))*(numel(vecDefinitionPrefOri) * numel(vecDefinitionSpatFreq));
intConnsLGN_to_CortON = sum(vecConnsPerTypeON(vecDefinitionCellTypes))*(numel(vecDefinitionPrefOri) * numel(vecDefinitionSpatFreq));
	
%create connection parameters from LGN to cortex
vecConductance_FromLGN_ToCort = [5.5 6]*0.15; %to [pyramid interneuron]
vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

matSynConnON_to_Cort = nan(intConnsLGN_to_CortON,2);
vecSynWeightON_to_Cort = nan(intConnsLGN_to_CortON,1);
vecSynConductanceON_to_Cort = nan(intConnsLGN_to_CortON,1);
vecSynDelayON_to_Cort = nan(intConnsLGN_to_CortON,1);

matSynConnOFF_to_Cort = nan(intConnsLGN_to_CortOFF,2);
vecSynWeightOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);
vecSynConductanceOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);
vecSynDelayOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);

intC_ON = 1;
intC_OFF = 1;
for intNeuron=1:intCortexCells
	%get cell parameters
	dblPsi = vecPrefPsi(intNeuron); %cell's phase offset
	dblTheta=vecPrefOri(intNeuron); %cell's preferred orientation
	dblPrefSF = vecPrefSF(intNeuron);
	dblPrefRF_X = vecPrefRF_X(intNeuron);
	dblPrefRF_Y = vecPrefRF_Y(intNeuron);
	
	intCellType = vecCellTypes(intNeuron); %pyramidal cell or interneuron
	intConnsON = vecConnsPerTypeON(intCellType);
	intConnsOFF = vecConnsPerTypeOFF(intCellType);
	
	%get assignment vectors
	vecConnsAssignON = intC_ON:(intC_ON+intConnsON-1);
	vecConnsAssignOFF = intC_OFF:(intC_OFF+intConnsOFF-1);
	intC_ON = intC_ON + intConnsON;
	intC_OFF = intC_OFF + intConnsOFF;
	
	%% create Gabor
	% get rotation matrices
	matX_theta=(matMeshX*cos(dblTheta)+matMeshY*sin(dblTheta))+dblPrefRF_X;
	matY_theta=(-matMeshX*sin(dblTheta)+matMeshY*cos(dblTheta))+dblPrefRF_Y;
	
	%get gabor
	matG = exp(-.5*(matX_theta.^2/dblSigmaX^2+matY_theta.^2/dblSigmaY^2)).*cos(2*pi*dblPrefSF*matX_theta+dblPsi);
	matPrefGabors(:,:,intNeuron) = matG;
	
	%get connection probabilities
	matProbConn_ON = max(0,matG) / sum(sum(max(0,matG)));
	matProbConn_OFF = max(0,-matG) / sum(sum(max(0,-matG)));
	
	%% ON fields
	%ON; linearize, get random connections, and reshape to matrix
	vecProbLinON = matProbConn_ON(:);
	vecProbCumSumON = cumsum(vecProbLinON);
	%matConnON = false(size(matProbConn_ON));
	vecConnsON = [];
	while numel(vecConnsON) < intConnsON
		vecConnsON = unique([vecConnsON sum(bsxfun(@gt,rand(1,intConnsON),vecProbCumSumON),1)+1]);
		if numel(vecConnsON) > intConnsON
			vecConnsON = vecConnsON(randperm(numel(vecConnsON),intConnsON));
		end
	end
	%[I,J] = ind2sub([intImX intImY],vecConnsON);
	%matConnON = logical(getFillGrid(zeros(intImX,intImY),I,J,true));
	
	%put in 2D connection matrix
	matConnections2D = [vecConnsON' repmat(intNeuron,[length(vecConnsON) 1])];
	matSynConnON_to_Cort(vecConnsAssignON,:) = matConnections2D;
	
	%put synaptic properties in vector
	vecSynWeightON_to_Cort(vecConnsAssignON) = vecProbLinON(vecConnsON);
	vecSynConductanceON_to_Cort(vecConnsAssignON) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignON));
	vecSynDelayON_to_Cort(vecConnsAssignON) = max(0,gaussrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignON)));
	
	%% OFF fields
	%OFF; linearize, get random connections
	vecProbLinOFF = matProbConn_OFF(:);
	vecProbCumSumOFF = cumsum(vecProbLinOFF);
	%matConnOFF = false(size(matProbConn_OFF));
	vecConnsOFF = [];
	while numel(vecConnsOFF) < intConnsOFF
		vecConnsOFF = unique([vecConnsOFF sum(bsxfun(@gt,rand(1,intConnsOFF),vecProbCumSumOFF),1)+1]);
		if numel(vecConnsOFF) > intConnsOFF
			vecConnsOFF = vecConnsOFF(randperm(numel(vecConnsOFF),intConnsOFF));
		end
	end
	%[I,J] = ind2sub([intImX intImY],vecConnsOFF);
	%matConnOFF = getFillGrid(zeros(intImX,intImY),I,J,true);
	
	%put in 2D connection matrix
	matConnections2D = [vecConnsOFF' repmat(intNeuron,[length(vecConnsOFF) 1])];
	matSynConnOFF_to_Cort(vecConnsAssignOFF,:) = matConnections2D;
	
	%put synaptic properties in vector
	vecSynWeightOFF_to_Cort(vecConnsAssignOFF) = vecProbLinOFF(vecConnsOFF);
	vecSynConductanceOFF_to_Cort(vecConnsAssignOFF) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignOFF));
	vecSynDelayOFF_to_Cort(vecConnsAssignOFF) = max(0,gaussrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignOFF)));
	
end

fprintf(' .. Created LGN to Cortex connectivity; %d active OFF field synapses, and %d active ON field synapses [%s]\n',intConnsLGN_to_CortOFF,intConnsLGN_to_CortON,getTime);

%% define LIF model parameters
%% neuron parameters
%pyramids
sParamCells{1}.dblThresh = -55; %mV
sParamCells{1}.dblTau_Peak = 1/1000; %ms
sParamCells{1}.dblV_E = 0; %excitatory reversal potential
sParamCells{1}.dblV_I = -70; %inhibitory reversal potential
sParamCells{1}.dblV_AHP = -90;%after-hyperpol reversal potential
sParamCells{1}.dblV_Leak = -65; %resting potential
sParamCells{1}.dblCm = 0.5; %capacitance
sParamCells{1}.dblG_Leak = 25; %nS
sParamCells{1}.dblG_AHP = 40; %nS

%interneurons
sParamCells{2}.dblThresh = -55; %mV
sParamCells{2}.dblTau_Peak = 2/1000; %ms
sParamCells{2}.dblV_E = 0; %excitatory reversal potential
sParamCells{2}.dblV_I = -70; %inhibitory reversal potential
sParamCells{2}.dblV_AHP = -90;%after-hyperpol reversal potential
sParamCells{2}.dblV_Leak = -65; %resting potential
sParamCells{2}.dblCm = 0.2; %capacitance
sParamCells{2}.dblG_Leak = 20; %nS
sParamCells{2}.dblG_AHP = 20; %nS

%define tau by cell type
vecTauPeakByType = [sParamCells{1}.dblTau_Peak sParamCells{2}.dblTau_Peak];

%% connection parameters
%number of connections
dblScalingFactor = 4;
matConnCortFromTo(1,:) = [40 40]*dblScalingFactor; %from pyramid to [pyr inter]
matConnCortFromTo(2,:) = [30 30]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
matConductancesFromTo(1,:) = [1.1 1.6]/dblScalingFactor; %from pyramid to [pyr inter]
matConductancesFromTo(2,:) = [1.5 1.0]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
dblDelayMeanCortToCort = 3/1000; %in ms
dblDelaySDCortToCort = 1/1000; %in ms

%msg
fprintf(' .. Created LIF model parameters [%s]\n',getTime);

%% create cortical connectivity and parameter vectors
%connection probability excitatory cells
vecConnProbSD = ang2rad([7.5 80]); %degrees difference in pref ori, for [pyramid interneuron]

%pre-allocate cortical synaptic connection matrix
vecTotalConnectionsPerType = sum(matConnCortFromTo,1);
intTotalCortConnections = sum(vecTotalConnectionsPerType(vecCellTypes));
matSynFromTo = nan(intTotalCortConnections,2);
vecSynConductance = nan(intTotalCortConnections,1);
vecSynDelay = nan(intTotalCortConnections,1);
vecSynExcInh = nan(intTotalCortConnections,1);

%pre-allocate cell parameters
vecCellThresh = nan(intCortexCells,1);
vecCellTauPeak = nan(intCortexCells,1);
vecCellV_E = nan(intCortexCells,1);
vecCellV_I = nan(intCortexCells,1);
vecCellV_AHP = nan(intCortexCells,1);
vecCellV_Leak = nan(intCortexCells,1);
vecCellCm =	nan(intCortexCells,1);
vecCellG_Leak = nan(intCortexCells,1);
vecCellG_AHP = nan(intCortexCells,1);

%get random connections
intC = 1;
for intNeuron=1:intCortexCells
	dblPrefOri = vecPrefOri(intNeuron);
	vecOriDiff = abs(circ_dist(vecPrefOri*2,dblPrefOri*2)/2)';
		
	%get this cell's type
	intTargetCellType = vecCellTypes(intNeuron); %pyramidal cell or interneuron
	vecConnsFromTypesToThisType = matConnCortFromTo(:,intTargetCellType); %from [pyr inter]
	
	%assign parameter values
	vecCellThresh(intNeuron) = sParamCells{intTargetCellType}.dblThresh;
	vecCellTauPeak(intNeuron) = sParamCells{intTargetCellType}.dblTau_Peak; %ms
	vecCellV_E(intNeuron) = sParamCells{intTargetCellType}.dblV_E; %excitatory reversal potential
	vecCellV_I(intNeuron) = sParamCells{intTargetCellType}.dblV_I; %inhibitory reversal potential
	vecCellV_AHP(intNeuron) = sParamCells{intTargetCellType}.dblV_AHP;%after-hyperpol reversal potential
	vecCellV_Leak(intNeuron) = sParamCells{intTargetCellType}.dblV_Leak; %resting potential
	vecCellCm(intNeuron) = sParamCells{intTargetCellType}.dblCm; %capacitance
	vecCellG_Leak(intNeuron) = sParamCells{intTargetCellType}.dblG_Leak; %nS
	vecCellG_AHP(intNeuron) = sParamCells{intTargetCellType}.dblG_AHP; %nS
	
	%loop through subtypes to make connections
	for intSubIdx=1:length(vecCellTypeLabels)
		intSourceCellType = vecCellTypeLabels(intSubIdx);
		indThisType = vecCellTypes==intSourceCellType;
		intConnsFromThisType = vecConnsFromTypesToThisType(intSourceCellType);
		dblConnProbOriSD = vecConnProbSD(intSourceCellType);
		vecConnProb = gausspdf(vecOriDiff,0,dblConnProbOriSD);
		vecConnProb(~indThisType) = 0; %remove other types
		vecConnProb(intNeuron) = 0; %remove the neuron itself
		vecConnProb = vecConnProb/sum(vecConnProb); %normalize so integral is 1
		vecProbCumSum = cumsum(vecConnProb);
		
		%get assignment vectors
		vecConnsAssign = intC:(intC+intConnsFromThisType-1);
		intC = intC + intConnsFromThisType;
		
		%make connections
		vecConns = [];
		while numel(vecConns) < intConnsFromThisType
			vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsFromThisType),vecProbCumSum),1)+1]);
			if numel(vecConns) > intConnsFromThisType
				vecConns = vecConns(randperm(numel(vecConns),intConnsFromThisType));
			end
		end
		
		%assign to connection matrix
		dblG = matConductancesFromTo(intSourceCellType,intTargetCellType);
		matSynFromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
		vecSynConductance(vecConnsAssign) = dblG;
		vecSynDelay(vecConnsAssign) = max(0,gaussrnd(dblDelayMeanCortToCort,dblDelaySDCortToCort,size(vecConnsAssign)));
		vecSynExcInh(vecConnsAssign) = intSourceCellType;
	end
end

fprintf(' .. Created cortical connectivity; %d active corticocortical synapses between %d cells [%s]\n',intTotalCortConnections,intCortexCells,getTime);

%% prepare variables for simulation
% run stitched together trials
vecThisV = ones(intCortexCells,1)*-65;
vecPrevV = vecThisV;
boolStimPresent = false;
intPrevTrial = 0;
intTrialT = 0;
vecOverallT = dblDeltaT:dblDeltaT:dblSimDur;
fprintf('\n   >>> Starting simulation run [%s]\n',getTime);
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
vecTauPeakByType = [1/1000 2/1000];
sData.vecOverallT = vecOverallT;
sData.dblDeltaT = dblDeltaT;
sData.matSynFromTo = matSynFromTo;
sData.dblSynSpikeMem = dblSynSpikeMem;
sData.vecSynExcInh = vecSynExcInh;
sData.intCortexCells = intCortexCells;
sData.vecSynDelay = vecSynDelay;
sData.vecSynConductance = vecSynConductance;
sData.vecCellThresh = vecCellThresh;
sData.vecTauPeakByType = vecTauPeakByType;
sData.vecCellV_E = vecCellV_E;
sData.vecCellV_I = vecCellV_I;
sData.vecCellV_AHP = vecCellV_AHP;
sData.vecCellV_Leak = vecCellV_Leak;
sData.vecCellCm = vecCellCm;
sData.vecCellG_Leak = vecCellG_Leak;
sData.vecCellG_AHP = vecCellG_AHP;
sData.vecSynConductanceON_to_Cort = vecSynConductanceON_to_Cort;
sData.vecSynConductanceOFF_to_Cort = vecSynConductanceOFF_to_Cort;
sData.vecSynWeightON_to_Cort = vecSynWeightON_to_Cort;
sData.vecSynWeightOFF_to_Cort = vecSynWeightOFF_to_Cort;
sData.vecSynDelayON_to_Cort = vecSynDelayON_to_Cort;
sData.vecSynDelayOFF_to_Cort = vecSynDelayOFF_to_Cort;
sData.matSynConnON_to_Cort = matSynConnON_to_Cort;
sData.matSynConnOFF_to_Cort = matSynConnOFF_to_Cort;
sData.matBlankLGN_ON = matBlankLGN_ON;
sData.matBlankLGN_OFF = matBlankLGN_OFF;
sData.cellLGN_ON = cellLGN_ON;
sData.cellLGN_OFF = cellLGN_OFF;
sData.vecTrialOris = vecTrialOris;
sData.vecTrialOriIdx = vecTrialOriIdx;
sData.vecTrialStimType = vecTrialStimType;
sData.vecTrialStartSecs = vecTrialStartSecs;
sData.vecStimStartSecs = vecStimStartSecs;
sData.vecTrialEndSecs = vecTrialEndSecs;
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
sData.dblStimDur = dblStimDur;

%cell preferences
sData.vecCellTypes = vecCellTypes;
sData.vecPrefPsi = vecPrefPsi; %cell's phase offset
sData.vecPrefOri = vecPrefOri; %cell's preferred orientation
sData.vecPrefSF = vecPrefSF;
sData.vecPrefRF_X = vecPrefRF_X;
sData.vecPrefRF_Y = vecPrefRF_Y;
sData.dblSigmaX = dblSigmaX; %length of gabor response
sData.dblSigmaY = dblSigmaY; %width of gabor response
sData.matPrefGabors = matPrefGabors;

%% run simulation
sData = getSimulationRun(sData);

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
strDataFile = ['D:\Data\Processed\V1_LIFmodel\Simulation_Orig' sStimParams.strStimType sprintf('%s.mat',getDate)];
fprintf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
save(strDataFile,'sData','sStimInputs','sStimParams');
printf(' .. Data saved! [%s]\n',getTime);

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
