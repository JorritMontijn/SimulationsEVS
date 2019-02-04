%{
function [dblCurT,... %only output
		vecThisV,boolStimPresent,intPrevTrial,intTrialT,intIter,... %input & output; line 1
		cellSpikeTimesLGN_ON,cellSpikeTimesLGN_OFF,cellSpikeTimesCortex,... %input & output; line 2
		vecSpikeCounterLGN_ON,vecSpikeCounterLGN_OFF,vecSpikeCounterCortex,... %input & output; line 3
		vecSpikeCounterPreAllocatedLGN_ON,vecSpikeCounterPreAllocatedLGN_OFF,vecSpikeCounterPreAllocatedCortex,intPreAllocationSize]... %input & output; line 4
		= getSimulationRun(vecOverallT,dblDeltaT,matCortConn,dblSynSpikeMem,vecCortSynType,intCortexCells,vecCortDelay,vecCortConductance,...%only input
		vecCellThresh,vecTauPeakByType,vecCellV_E,vecCellV_I,vecCellV_AHP,vecCellV_Leak,vecCellCm,vecCellG_Leak,vecCellG_AHP,...
		vecSynConductanceON_to_Cort,vecSynConductanceOFF_to_Cort,vecSynWeightON_to_Cort,vecSynWeightOFF_to_Cort,vecSynDelayON_to_Cort,vecSynDelayOFF_to_Cort,...%only input
		matSynConnON_to_Cort,matSynConnOFF_to_Cort,matBlankLGN_ON,matBlankLGN_OFF,cellLGN_ON,cellLGN_OFF,... %only input
		vecTrialOris,vecTrialOriIdx,vecStimStartSecs,vecTrialEndSecs,...  %only input
		vecThisV,boolStimPresent,intPrevTrial,intTrialT,intIter,...  %input & output; line 1
		cellSpikeTimesLGN_ON,cellSpikeTimesLGN_OFF,cellSpikeTimesCortex,... %input & output; line 2
		vecSpikeCounterLGN_ON,vecSpikeCounterLGN_OFF,vecSpikeCounterCortex,... %input & output; line 3
		vecSpikeCounterPreAllocatedLGN_ON,vecSpikeCounterPreAllocatedLGN_OFF,vecSpikeCounterPreAllocatedCortex,intPreAllocationSize) %input & output; line 4
%}

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
function sData = getDynSimRunIndInp(sData,dblMaxRunningTime,intWorker)
	%% prep
	%check if we're on cluster
	global boolClust;
	global boolSaveVm;
	if isempty(boolClust),boolClust=false;end
	if isempty(boolSaveVm),boolSaveVm=false;end
	hTic = tic;
	
	%unpack
	intCortexCells = sData.intCortexCells;
	intCellsV1 = sData.intCellsV1;
	intCellsV2 = sData.intCellsV2;
	dblDeltaT = sData.dblDeltaT;
	
	matSynFromTo = sData.matSynFromTo;
	dblSynSpikeMem = sData.dblSynSpikeMem;
	vecSynExcInh = sData.vecSynExcInh;
	vecSynDelay = sData.vecSynDelay;
	vecSynConductance = sData.vecSynConductance;
	vecSynWeight = sData.vecSynWeight;
	vecSynType = sData.vecSynType;
	
	vecCellTypes = sData.vecCellTypes';
	vecCellThresh = sData.vecCellThresh;
	vecCellTauPeak = sData.vecCellTauPeak;
	vecCellV_E = sData.vecCellV_E;
	vecCellV_I = sData.vecCellV_I;
	vecCellV_AHP = sData.vecCellV_AHP;
	vecCellV_Leak = sData.vecCellV_Leak;
	vecCellCm = sData.vecCellCm;
	vecCellG_Leak = sData.vecCellG_Leak;
	vecCellG_AHP = sData.vecCellG_AHP;
	%vecCellRefracT = max(dblDeltaT,sData.vecCellRefracT);
	if isfield(sData,'vecTauPeakByType')
		vecTauPeakByType = sData.vecTauPeakByType;
	else
		vecTauPeakByType = unique(vecCellTauPeak);
	end
	if isfield(sData,'dblAttention'),dblAttention = sData.dblAttention;else dblAttention = 0;end
	if isfield(sData,'intAttArea'),intAttArea = sData.intAttArea;else intAttArea = 0;end
	
	vecSynConductanceON_to_Cort = sData.vecSynConductanceON_to_Cort;
	vecSynConductanceOFF_to_Cort = sData.vecSynConductanceOFF_to_Cort;
	vecSynWeightON_to_Cort = sData.vecSynWeightON_to_Cort;
	vecSynWeightOFF_to_Cort = sData.vecSynWeightOFF_to_Cort;
	vecSynDelayON_to_Cort = sData.vecSynDelayON_to_Cort;
	vecSynDelayOFF_to_Cort = sData.vecSynDelayOFF_to_Cort;
	matSynConnON_to_Cort = sData.matSynConnON_to_Cort;
	matSynConnOFF_to_Cort = sData.matSynConnOFF_to_Cort;
	matBlankLGN_ON = sData.matBlankLGN_ON;
	matBlankLGN_OFF = sData.matBlankLGN_OFF;
	
	vecThisV = sData.vecThisV;
	boolStimPresent = sData.boolStimPresent;
	intPrevTrial = sData.intPrevTrial;
	intTrialT = sData.intTrialT;
	cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	vecSpikeCounterLGN_ON = sData.vecSpikeCounterLGN_ON;
	vecSpikeCounterLGN_OFF = sData.vecSpikeCounterLGN_OFF;
	vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
	intPreAllocationSize = sData.intPreAllocationSize;
	
	%trial timing
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	
	%trial stim type
	vecTrialStimType = sData.vecTrialStimType;
	vecTrialStimRep = sData.vecTrialStimRep;
	
	%stim type props
	vecStimTypeOris = sData.vecStimTypeOris;
	vecStimTypeOriNoise = sData.vecStimTypeOriNoise;
	vecStimTypeSFs = sData.vecStimTypeSFs;
	vecStimTypeSFNoise = sData.vecStimTypeSFNoise;
	vecStimTypeTFs = sData.vecStimTypeTFs;
	vecStimTypeTFNoise = sData.vecStimTypeTFNoise;
	vecStimTypeContrasts = sData.vecStimTypeContrasts;
	vecStimTypeContrastNoise = sData.vecStimTypeContrastNoise;
	vecStimTypeLuminances = sData.vecStimTypeLuminances;
	vecStimTypeLuminanceNoise = sData.vecStimTypeLuminanceNoise;
	vecStimTypePhase = sData.vecStimTypePhase;
	vecStimTypePhaseNoise = sData.vecStimTypePhaseNoise;
	vecStimTypeGain = sData.vecStimTypeGain;
	vecStimTypeGainNoise = sData.vecStimTypeGainNoise;
	
	%get stimulus parameters
	sStimParams = sData.sStimParams;
	
	%define vars
	vecSpikeCounterPreAllocatedLGN_ON = vecSpikeCounterLGN_ON;
	vecSpikeCounterPreAllocatedLGN_OFF = vecSpikeCounterLGN_OFF;
	vecSpikeCounterPreAllocatedCortex = vecSpikeCounterCortex;
	vecCortTauSyn = vecTauPeakByType(vecSynExcInh);
	
	%synapse function
	fPSP = @(dblT,vecSpikeTimes,dblTauPeak) sum(max(0,dblT - vecSpikeTimes).*(1/dblTauPeak).*exp(1-(dblT - vecSpikeTimes)/dblTauPeak));
	fPSP2 = @(matT,dblTauPeak) sum(max(0,matT).*(1/dblTauPeak).*exp(1-matT/dblTauPeak),2);
	
	%create synaptic target cell array
	cellSynTargets = cell(1,intCortexCells);
	for intNeuron=1:intCortexCells
		cellSynTargets{intNeuron} = find(matSynFromTo(:,1) == intNeuron);
	end
	
	%create LGN source cell array for V1
	cellSynONSourcesV1 = cell(1,intCellsV1);
	for intNeuron=1:intCellsV1
		cellSynONSourcesV1{intNeuron} = find(matSynConnON_to_Cort(:,2) == intNeuron);
	end
	cellSynOFFSourcesV1 = cell(1,intCellsV1);
	for intNeuron=1:intCellsV1
		cellSynOFFSourcesV1{intNeuron} = find(matSynConnOFF_to_Cort(:,2) == intNeuron);
	end
	
	%% attention & pixel noise
	indAttentionTargetCells = vecCellTypes==2;
	if intAttArea==0,indAttentionTargetAreas=isnumeric(sData.vecCellArea');
	else indAttentionTargetAreas = sData.vecCellArea'==intAttArea;end
	vecAttention = dblAttention*indAttentionTargetAreas.*indAttentionTargetCells;
	dblPixNoise = sData.dblPixNoise;
	
	%% create recent spike array Cortex
	%build cylindrical recent spike storage
	intMaxRecentSpikes = 100;
	matRecentSpikes = nan(intCortexCells,intMaxRecentSpikes);
	vecRSCounter = zeros(1,intCortexCells);
	vecRSPos = ones(1,intCortexCells);
	dblCurT = 0;
	
	%get recent spikes
	for intN = 1:intCortexCells
		vecRecentSpikes = cellSpikeTimesCortex{intN}(((cellSpikeTimesCortex{intN} - dblCurT) + dblSynSpikeMem) > 0);
		if ~isempty(vecRecentSpikes)
			intNumRecentSpikes = numel(vecRecentSpikes);
			matRecentSpikes(intN,vecRSPos(intN):(vecRSPos(intN)+intNumRecentSpikes-1)) = vecRecentSpikes;
			vecRSPos(intN) = vecRSPos(intN) + intNumRecentSpikes;
			vecRSCounter(intN) = vecRSCounter(intN) + intNumRecentSpikes;
		end
	end
	
	%% copy trial based vars
	vecTrialStimRepIn = vecTrialStimRep;
	vecTrialStimTypeIn = vecTrialStimType;
	dblOffset = vecTrialStartSecs(1);
	vecTrialStartSecsIn = vecTrialStartSecs-dblOffset;
	vecStimStartSecsIn = vecStimStartSecs-dblOffset;
	vecStimStopSecsIn = vecStimStopSecs-dblOffset;
	vecTrialEndSecsIn = vecTrialEndSecs-dblOffset;
	dblRepDur = vecTrialEndSecsIn(end);
	
	%% randomize first run
	if boolSaveVm,matVm = nan(intCortexCells,10000,7);end
	%trial stim type
	vecTrialStimType = randperm(numel(vecTrialStimTypeIn));
	vecTrialStimRep = vecTrialStimRepIn;
	
	%stim type props
	vecTrialOris = vecStimTypeOris(vecTrialStimType);
	vecTrialOriNoise = vecStimTypeOriNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialSFs = vecStimTypeSFs(vecTrialStimType);
	vecTrialSFNoise = vecStimTypeSFNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialTFs = vecStimTypeTFs(vecTrialStimType);
	vecTrialTFNoise = vecStimTypeTFNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialContrasts = vecStimTypeContrasts(vecTrialStimType);
	vecTrialContrastNoise = vecStimTypeContrastNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialLuminances = vecStimTypeLuminances(vecTrialStimType);
	vecTrialLuminanceNoise = vecStimTypeLuminanceNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialPhases = vecStimTypePhase(vecTrialStimType);
	vecTrialPhaseNoise = vecStimTypePhaseNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	vecTrialGains = vecStimTypeGain(vecTrialStimType);
	vecTrialGainNoise = vecStimTypeGainNoise(vecTrialStimType).*randn(size(vecTrialStimType));
	
	%% check running type
	if isa(dblMaxRunningTime,'uint64')
		intMaxRep = dblMaxRunningTime;
		dblMaxRunningTime = inf;
	else
		intMaxRep = inf;
	end
	
	%% run
	%vecTimeSinceSpike = (dblCurT + max(vecCellRefracT(:)))*ones(size(vecThisV));
	intRepCounter = 0;
	intRepIter = 0;
	intIter = 0;
	vecOverallT = [];
	intTimeStepCounter = 0;
	while toc(hTic) < dblMaxRunningTime && intRepCounter < intMaxRep
		
		%increment repetition counters
		intIter = intIter + intRepIter;
		if intIter > 0
			intRepCounter = intRepCounter + 1;
			intRepIter = 0;
			vecOverallT = cat(2,vecOverallT,vecOverallRepT);
			
			%build stimulus vectors for next repetition and prep random noise
			vecRandStimTypes = randperm(numel(vecTrialStimTypeIn));
			
			%concatenate vectors
			vecTrialStimType = cat(2,vecTrialStimType,vecRandStimTypes);
			vecTrialStimRep = cat(2,vecTrialStimRep,vecTrialStimRepIn*(intRepCounter+1));
			
			vecTrialStartSecs = cat(2,vecTrialStartSecs,vecTrialStartSecsIn+dblRepDur*intRepCounter);
			vecStimStartSecs = cat(2,vecStimStartSecs,vecStimStartSecsIn+dblRepDur*intRepCounter);
			vecStimStopSecs = cat(2,vecStimStopSecs,vecStimStopSecsIn+dblRepDur*intRepCounter);
			vecTrialEndSecs = cat(2,vecTrialEndSecs,vecTrialEndSecsIn+dblRepDur*intRepCounter);
			
			%stim type props
			vecTrialOris = cat(2,vecTrialOris,vecStimTypeOris(vecRandStimTypes));
			vecTrialOriNoise = cat(2,vecTrialOriNoise,vecStimTypeOriNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
			vecTrialSFs = cat(2,vecTrialSFs,vecStimTypeSFs(vecRandStimTypes));
			vecTrialSFNoise = cat(2,vecTrialSFNoise,vecStimTypeSFNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
			vecTrialTFs = cat(2,vecTrialTFs,vecStimTypeTFs(vecRandStimTypes));
			vecTrialTFNoise = cat(2,vecTrialTFNoise,vecStimTypeTFNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
			vecTrialContrasts = cat(2,vecTrialContrasts,vecStimTypeContrasts(vecRandStimTypes));
			vecTrialContrastNoise = cat(2,vecTrialContrastNoise,vecStimTypeContrastNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
			vecTrialLuminances = cat(2,vecTrialLuminances,vecStimTypeLuminances(vecRandStimTypes));
			vecTrialLuminanceNoise = cat(2,vecTrialLuminanceNoise,vecStimTypeLuminanceNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
			vecTrialPhases = cat(2,vecTrialPhases,vecStimTypePhase(vecRandStimTypes));
			vecTrialPhaseNoise =  cat(2,vecTrialPhaseNoise,vecStimTypePhaseNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
            vecTrialGains = cat(2,vecTrialGains,vecStimTypeGain(vecRandStimTypes));
			vecTrialGainNoise =  cat(2,vecTrialGainNoise,vecStimTypeGainNoise(vecRandStimTypes).*randn(size(vecRandStimTypes)));
	
		end
		%msg
		printf('   > Worker %d: Starting repetition %d/%d; %.3fs remaining... [%s]\n',intWorker,intRepCounter,intMaxRep,dblMaxRunningTime-toc(hTic),getTime);
		
		%build timing for next stretch
		dblNextRepStartT = dblCurT + dblDeltaT;
		dblRepEnd = dblNextRepStartT + dblRepDur - dblDeltaT;
		vecOverallRepT = dblNextRepStartT:dblDeltaT:dblRepEnd;
		hTimer=tic;
		%%
		for dblCurT = vecOverallRepT
			%% check which trial
			intTimeStepCounter = intTimeStepCounter + 1;
			intRepIter = intRepIter + 1;
			intThisTrial = sum(dblCurT>vecTrialStartSecs);
			
			%% prepare trial
			%start new trial
			if intThisTrial > intPrevTrial
				%% get trial data
				intPrevTrial = intThisTrial;
				
				%set counter for this trial
				intTrialT = 0;
				boolStimPresent = true;
				
				%% generate stimulus
				%generate static sStimParams
				sSP = struct;
				sSP.dblPixNoise = dblPixNoise;
				sSP.intFrameRate=sStimParams.intFrameRate;
				sSP.dblAngleInDeg=vecTrialOris(intThisTrial) + vecTrialOriNoise(intThisTrial);
				sSP.dblSF=vecTrialSFs(intThisTrial) + vecTrialSFNoise(intThisTrial);
				sSP.dblTF=vecTrialTFs(intThisTrial) + vecTrialTFNoise(intThisTrial);
				sSP.dblStimSizeRetDeg=sStimParams.dblStimSizeRetDeg;
				sSP.vecScrPixWidthHeight=sStimParams.vecScrPixWidthHeight;
				sSP.vecScrDegWidthHeight=sStimParams.vecScrDegWidthHeight;
				sSP.dblDeltaT=sStimParams.dblDeltaT;
				sSP.dblStimDur=sStimParams.dblStimDur;
				sSP.dblBlankDur=sStimParams.dblBlankDur;
				sSP.varDeltaSyn=sStimParams.varDeltaSyn;
				sSP.strStimType = sStimParams.strStimType;
				sSP.dblContrast=vecTrialContrasts(intThisTrial) + vecTrialContrastNoise(intThisTrial);
				sSP.dblLuminance=vecTrialLuminances(intThisTrial) + vecTrialLuminanceNoise(intThisTrial);
				sSP.dblPhase=vecTrialPhases(intThisTrial) + vecTrialPhaseNoise(intThisTrial);
				sSP.dblGain=vecTrialGains(intThisTrial) + vecTrialGainNoise(intThisTrial);
				
				%msg
				printf('   > Worker %d: Generating orientation %.1f degs, %.0fms, SF=%.3f, TF=%.3f, PN=%.3f, Gain=%.3f; %.3fs remaining... [%s]\n',...
					intWorker,sSP.dblAngleInDeg,sSP.dblStimDur*1000,sSP.dblSF,sSP.dblTF,sSP.dblPixNoise,sSP.dblGain,...
					dblMaxRunningTime-toc(hTic),getTime);
		
				%get stim drive
				sStimDrive = getDynBottomUpInputs(sSP);
				
				%assign
				matLGN_ON = sStimDrive.matLGN_ON;
				matLGN_OFF = sStimDrive.matLGN_OFF;
			end
			
			%check to stop trial
			intStopTrial = sum(dblCurT>vecTrialEndSecs);
			if intStopTrial >= intThisTrial
				boolStimPresent = false;
			end
			
			% check if stimulus should be presented
			intTrialT = intTrialT + 1;
			if boolStimPresent && intTrialT < size(matLGN_ON,3)
				%create stochastic spiking
				matThisIterON = matLGN_ON(:,:,intTrialT);
				matThisIterOFF = matLGN_OFF(:,:,intTrialT);
			else
				matThisIterON = matBlankLGN_ON;
				matThisIterOFF = matBlankLGN_OFF;
			end
			
			
			%% calculate all excitatory currents due to spikes from LGN
			intSubType = 1;
			dblTauPeak = vecTauPeakByType(intSubType);
			intStepsT = dblSynSpikeMem/dblDeltaT;
			vecDelayT = dblDeltaT:dblDeltaT:dblSynSpikeMem;
			vecLGNCurrents = zeros(size(vecThisV));
			for intField=[1 2]
				if intField == 1
					vecSW = vecSynWeightON_to_Cort; %#ok<NASGU>
					vecSC = vecSynConductanceON_to_Cort;
					matC = matSynConnON_to_Cort;
					cellSynSource = cellSynONSourcesV1;
					matThisIterDrive = matThisIterON;
				else
					vecSW = vecSynWeightOFF_to_Cort; %#ok<NASGU>
					vecSC = vecSynConductanceOFF_to_Cort;
					matC = matSynConnOFF_to_Cort;
					cellSynSource = cellSynOFFSourcesV1;
					matThisIterDrive = matThisIterOFF;
				end
				
				%loop
				for intTargetNeuron=1:intCellsV1
					%get cells
					vecSynapses = cellSynSource{intTargetNeuron};
					vecSourceCells = matC(vecSynapses,1);
					
					%generate spikes
					vecDrive = matThisIterDrive(vecSourceCells);
					matSpikeTimes = bsxfun(@times,vecDelayT,bsxfun(@gt,vecDrive,rand(numel(vecDrive),intStepsT)));
					
					%get currents per synapse
					vecPSP = fPSP2(matSpikeTimes,dblTauPeak);
					
					%get neuron based properties
					dblTotalCurrent = sum(vecSW(vecSynapses) .* vecSC(vecSynapses) .* vecPSP);
					
					%calculate total current
					vecLGNCurrents(intTargetNeuron) = vecLGNCurrents(intTargetNeuron) + dblTotalCurrent;
				end
			end
			
			%%
			% remove old spikes
			for intN = find(vecRSCounter>0)
				intPos = getCycPos(vecRSPos(intN),vecRSCounter(intN),intMaxRecentSpikes);
				if dblCurT - matRecentSpikes(intN,intPos) > dblSynSpikeMem
					matRecentSpikes(intN,intPos) = nan;
					vecRSCounter(intN) = vecRSCounter(intN) - 1;
				end
			end
			
			%%
			%get PSPs for all spiking neurons
			vecPSPs = zeros(size(matSynFromTo,1),1);
			for intNeuron=find(vecRSCounter>0)
				%get recent spikes
				vecPos = getCycPosVec(vecRSPos(intNeuron),vecRSCounter(intNeuron),intMaxRecentSpikes);
				vecSpikeTimes = matRecentSpikes(intNeuron,vecPos);
				vecTargetSynapses = cellSynTargets{intNeuron};
				if isempty(vecTargetSynapses),continue;end
				%get neuron based properties
				dblTauPeak = vecCortTauSyn(vecTargetSynapses(1));
				
				%get synapse data
				%dblWeight = 1;
				
				%calculate PSPs
				vecPSPs(vecTargetSynapses) = vecSynWeight(vecTargetSynapses) .* vecSynConductance(vecTargetSynapses) .* fPSP2(bsxfun(@minus,dblCurT-vecSynDelay(vecTargetSynapses),vecSpikeTimes),dblTauPeak);
			end
			
			
			% sum all transmissions incoming to each neuron
			indNonZero = vecPSPs~=0;
			indExcSynapses = vecSynExcInh==1;
			indInhSynapses = vecSynExcInh==2;
			vecExcCurrents = accumarray(matSynFromTo(indExcSynapses&indNonZero,2),vecPSPs(indExcSynapses&indNonZero),[intCortexCells 1],@sum,0);
			vecInhCurrents = accumarray(matSynFromTo(indInhSynapses&indNonZero,2),vecPSPs(indInhSynapses&indNonZero),[intCortexCells 1],@sum,0);
			
			
			% calculate after hyper polarizations
			intSubType = 2;
			dblThisTau =  vecTauPeakByType(intSubType);
			vecAHPs = zeros(intCortexCells,1);
			for intNeuron=1:intCortexCells
				%get recent spikes
				vecPos = getCycPosVec(vecRSPos(intNeuron),vecRSCounter(intNeuron),intMaxRecentSpikes);
				vecSpikeTimes = matRecentSpikes(intNeuron,vecPos);
				if ~isempty(vecSpikeTimes)
					%get synapse data
					dblWeight = 1;%vecSW(intSynapse);
					dblDelay = 0;
					%dblTauPeak = dblThisTau;
					dblCond = vecCellG_AHP(intNeuron);
					
					%calculate PSPs
					vecAHPs(intNeuron) = max(0,dblWeight * dblCond * fPSP(dblCurT-dblDelay,vecSpikeTimes,dblThisTau));
				end
			end
			
			% integrate all inputs
			vecPrevV = vecThisV;
			vecThisV = vecPrevV + dblDeltaT*(...
				-vecLGNCurrents.*(vecPrevV-vecCellV_E) ...
				-vecExcCurrents.*(vecPrevV-vecCellV_E) ...
				-vecAttention.*(vecPrevV-vecCellV_E) ...
				-vecInhCurrents.*(vecPrevV-vecCellV_I) ...
				-vecCellG_Leak.*(vecPrevV-vecCellV_Leak) ...
				-vecAHPs.*(vecPrevV-vecCellV_AHP) ...
				) ./ vecCellCm;
			if any(isnan(vecThisV) | isinf(vecThisV));
				vecErrN = find(isnan(vecThisV) | isinf(vecThisV));
				intErrs = numel(vecErrN);
				intErrN = vecErrN(randi(intErrs));
				if boolClust
					printf('\n.. [%d] ERRORS DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrs,intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
				else
					fprintf('\n.. [%d] ERRORS DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrs,intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
					warning([mfilename ':VmOutOfBounds'],'Membrane voltage out of bounds!');
				end
				break;
			end
			
			if boolSaveVm
				if size(matVm,2) < intTimeStepCounter
					matVm(:,(end+1):(end+10000),:) = nan;
				end
				matVm(:,intTimeStepCounter,1) = vecThisV;
				matVm(:,intTimeStepCounter,2) = vecLGNCurrents.*(vecPrevV-vecCellV_E);
				matVm(:,intTimeStepCounter,3) = vecExcCurrents.*(vecPrevV-vecCellV_E);
				matVm(:,intTimeStepCounter,4) = vecInhCurrents.*(vecPrevV-vecCellV_I);
				matVm(:,intTimeStepCounter,5) = vecCellG_Leak.*(vecPrevV-vecCellV_Leak);
				matVm(:,intTimeStepCounter,6) = vecAHPs.*(vecPrevV-vecCellV_AHP);
				matVm(:,intTimeStepCounter,7) = vecAttention.*(vecPrevV-vecCellV_E);
			end
			%refractory period
			%indRefrac = (vecTimeSinceSpike < vecCellRefracT);
			%if any(indRefrac)
			%	vecThisV(indRefrac) = min(vecCellV_Leak(indRefrac),vecThisV(indRefrac));
			%end
			
			% spiking
			%vecTimeSinceSpike = vecTimeSinceSpike + dblDeltaT;
			indSpiking = vecThisV>=vecCellThresh;
			if any(indSpiking)
				for intN = find(indSpiking)'
					%check number of spikes and pre-allocated size
					vecSpikeCounterCortex(intN) = vecSpikeCounterCortex(intN) + 1;
					if vecSpikeCounterPreAllocatedCortex(intN) < vecSpikeCounterCortex(intN);
						vecSpikeCounterPreAllocatedCortex(intN) = vecSpikeCounterPreAllocatedCortex(intN) + intPreAllocationSize;
						cellSpikeTimesCortex{intN} = [cellSpikeTimesCortex{intN} nan(1,intPreAllocationSize)];
					end
					
					% add spike time
					cellSpikeTimesCortex{intN}(vecSpikeCounterCortex(intN)) = dblCurT;
					vecThisV(intN) = vecCellV_Leak(intN);
					%vecTimeSinceSpike(intN) = 0;
					
					%add recent spikes
					matRecentSpikes(intN,vecRSPos(intN)) = dblCurT;
					vecRSPos(intN) = mod(vecRSPos(intN) + 1,intMaxRecentSpikes);
					vecRSCounter(intN) = vecRSCounter(intN) + 1;
				end
			end
			vecRSPos(vecRSPos==0)=intMaxRecentSpikes; %mod(y,y) is 0, but we want it to be y
			
			%% check running time
			if toc(hTic) > dblMaxRunningTime || intRepCounter > intMaxRep
				%remove current unfinished repetition
				%{
				indRemTrials = vecTrialStimRep==max(vecTrialStimRep);
				vecTrialStimType(indRemTrials) = [];
				vecTrialStimRep(indRemTrials) = [];
				vecTrialStartSecs(indRemTrials) = [];
				vecStimStartSecs(indRemTrials) = [];
				vecStimStopSecs(indRemTrials) = [];
				vecTrialEndSecs(indRemTrials) = [];
				%}
				%break
				break;
			end
			if ispc && toc(hTimer) > 5
				hTimer = tic;
				vecSpikeCounts = cellfun(@numel,cellSpikeTimesCortex)-cellfun(@sum,cellfun(@isnan,cellSpikeTimesCortex,'UniformOutput',false));
				printf('Elapsed: %.1fs; now at t=%.3fs; mean rate (Hz): %.3f (V1 Pyr); %.3f (V1 Int); %.3f (V2 Pyr); %.3f (V2 Int) [%s]\n',...
					toc(hTic),dblCurT,...
					mean(vecSpikeCounts(vecCellTypes==1 & sData.vecCellArea'==1))/dblCurT,...
					mean(vecSpikeCounts(vecCellTypes==2 & sData.vecCellArea'==1))/dblCurT,...
					mean(vecSpikeCounts(vecCellTypes==1 & sData.vecCellArea'==2))/dblCurT,...
					mean(vecSpikeCounts(vecCellTypes==2 & sData.vecCellArea'==2))/dblCurT,...
					getTime);
				cla;
				indPlotCells = vecCellTypes==1 & sData.vecCellArea'==1;
				scatter(1:sum(indPlotCells),vecSpikeCounts(indPlotCells)/dblCurT);
				ylim([0 max(get(gca,'ylim'))]);drawnow;
				
			end
		end
	end
	
	%end msg
	hToc = toc(hTic);
	if boolClust
		printf('Done! Total run took %f seconds [%s]\n',hToc,getTime);
	else
		fprintf('Done! Total run took %f seconds [%s]\n',hToc,getTime);
	end
	
	%% repack
	sData.intCellsV1 = intCellsV1;
	sData.intCellsV2 = intCellsV2;
	sData.intCortexCells = intCortexCells;
	sData.vecOverallT = vecOverallT;
	sData.dblDeltaT = dblDeltaT;
	
	sData.matSynFromTo = matSynFromTo;
	sData.dblSynSpikeMem = dblSynSpikeMem;
	sData.vecSynExcInh = vecSynExcInh;
	sData.vecSynDelay = vecSynDelay;
	sData.vecSynWeight = vecSynWeight;
	sData.vecSynConductance = vecSynConductance;
	sData.vecSynType = vecSynType;
	
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
	
	sData.vecThisV = vecThisV;
	sData.boolStimPresent = double(boolStimPresent);
	sData.intPrevTrial = intPrevTrial;
	sData.intTrialT = intTrialT;
	sData.intIter = intRepIter;
	sData.cellSpikeTimesLGN_ON = cellSpikeTimesLGN_ON;
	sData.cellSpikeTimesLGN_OFF = cellSpikeTimesLGN_OFF;
	sData.cellSpikeTimesCortex = cellSpikeTimesCortex;
	sData.vecSpikeCounterLGN_ON = vecSpikeCounterLGN_ON;
	sData.vecSpikeCounterLGN_OFF = vecSpikeCounterLGN_OFF;
	sData.vecSpikeCounterCortex = vecSpikeCounterCortex;
	sData.intPreAllocationSize = intPreAllocationSize;
	
	%actual trial data
	sData.vecTrialStimRep = vecTrialStimRep;
	sData.vecTrialStimType = vecTrialStimType;
	sData.vecTrialStartSecs = vecTrialStartSecs;
	sData.vecStimStartSecs = vecStimStartSecs;
	sData.vecStimStopSecs = vecStimStopSecs;
	sData.vecTrialEndSecs = vecTrialEndSecs;
	
	sData.vecTrialOris = vecTrialOris;
	sData.vecTrialOriNoise = vecTrialOriNoise;
	sData.vecTrialSFs = vecTrialSFs;
	sData.vecTrialSFNoise = vecTrialSFNoise;
	sData.vecTrialTFs = vecTrialTFs;
	sData.vecTrialTFNoise = vecTrialTFNoise;
	sData.vecTrialContrasts = vecTrialContrasts;
	sData.vecTrialContrastNoise = vecTrialContrastNoise;
	sData.vecTrialLuminances = vecTrialLuminances;
	sData.vecTrialLuminanceNoise = vecTrialLuminanceNoise;
    sData.vecTrialPhases = vecTrialPhases;
	sData.vecTrialPhaseNoise = vecTrialPhaseNoise;
    sData.vecTrialGains = vecTrialGains;
	sData.vecTrialGainNoise = vecTrialGainNoise;
	
	if boolSaveVm
		matVm(:,isnan(matVm(1,:,1)),:) = [];
		sData.matVm = matVm;
	end
end
