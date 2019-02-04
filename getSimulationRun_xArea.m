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
function sData = getSimulationRun_xArea(sData)
	%check if we're on cluster
	global boolClust;
	if isempty(boolClust),boolClust=false;end
	
	%unpack
	intCortexCells = sData.intCortexCells;
	intCellsV1 = sData.intCellsV1;
	intCellsV2 = sData.intCellsV2;
	vecOverallT = sData.vecOverallT;
	dblDeltaT = sData.dblDeltaT;
	
	matSynFromTo = sData.matSynFromTo;
	dblSynSpikeMem = sData.dblSynSpikeMem;
	vecSynExcInh = sData.vecSynExcInh;
	vecSynDelay = sData.vecSynDelay;
	vecSynConductance = sData.vecSynConductance;
	vecSynWeight = sData.vecSynWeight;
	vecSynType = sData.vecSynType;
	
	vecCellThresh = sData.vecCellThresh;
	vecTauPeakByType = sData.vecTauPeakByType;
	vecCellV_E = sData.vecCellV_E;
	vecCellV_I = sData.vecCellV_I;
	vecCellV_AHP = sData.vecCellV_AHP;
	vecCellV_Leak = sData.vecCellV_Leak;
	vecCellCm = sData.vecCellCm;
	vecCellG_Leak = sData.vecCellG_Leak;
	vecCellG_AHP = sData.vecCellG_AHP;
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
	cellLGN_ON = sData.cellLGN_ON;
	cellLGN_OFF = sData.cellLGN_OFF;
	vecTrialOris = sData.vecTrialOris;
	vecTrialOriIdx = sData.vecTrialOriIdx;
	vecTrialStimType = sData.vecTrialStimType;
	vecStimStartSecs = sData.vecStimStartSecs;
	vecStimStopSecs = sData.vecStimStopSecs;
	vecTrialStartSecs = sData.vecTrialStartSecs;
	vecTrialEndSecs = sData.vecTrialEndSecs;
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
	
	%create synaptic target cell array for LGN
	intCellsLGN = numel(cellSpikeTimesLGN_ON);
	cellSynTargetsLGN_ON = cell(1,intCellsLGN);
	for intNeuron=1:intCellsLGN
		cellSynTargetsLGN_ON{intNeuron} = find(matSynConnON_to_Cort(:,1) == intNeuron);
	end
	indMootCells_ON = cellfun(@numel,cellSynTargetsLGN_ON)==0;
	cellSynTargetsLGN_OFF = cell(1,intCellsLGN);
	for intNeuron=1:intCellsLGN
		cellSynTargetsLGN_OFF{intNeuron} = find(matSynConnOFF_to_Cort(:,1) == intNeuron);
	end
	indMootCells_OFF = cellfun(@numel,cellSynTargetsLGN_OFF)==0;
	
	%% create recent spike array Cortex
	%build cylindrical recent spike storage
	intMaxRecentSpikes = 100;
	matRecentSpikes = nan(intCortexCells,intMaxRecentSpikes);
	vecRSCounter = zeros(1,intCortexCells);
	vecRSPos = ones(1,intCortexCells);
	dblCurT = vecOverallT(end);%vecOverallT(1);
	
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
	
	
	%% create recent spike arrays LGN ON
	%build cylindrical recent spike storage
	matRecentSpikes_ON = nan(intCellsLGN,intMaxRecentSpikes);
	vecRSCounter_ON = zeros(1,intCellsLGN);
	vecRSPos_ON = ones(1,intCellsLGN);
	
	%get recent spikes
	for intN = 1:intCellsLGN
		vecRecentSpikes_ON = cellSpikeTimesLGN_ON{intN}(((cellSpikeTimesLGN_ON{intN} - dblCurT) + dblSynSpikeMem) > 0);
		if ~isempty(vecRecentSpikes_ON)
			intNumRecentSpikes = numel(vecRecentSpikes_ON);
			matRecentSpikes_ON(intN,vecRSPos_ON(intN):(vecRSPos_ON(intN)+intNumRecentSpikes-1)) = vecRecentSpikes_ON;
			vecRSPos_ON(intN) = vecRSPos_ON(intN) + intNumRecentSpikes;
			vecRSCounter_ON(intN) = vecRSCounter_ON(intN) + intNumRecentSpikes;
		end
	end
	
	%% create recent spike arrays LGN OFF
	%build cylindrical recent spike storage
	matRecentSpikes_OFF = nan(intCellsLGN,intMaxRecentSpikes);
	vecRSCounter_OFF = zeros(1,intCellsLGN);
	vecRSPos_OFF = ones(1,intCellsLGN);
	
	%get recent spikes
	for intN = 1:intCellsLGN
		vecRecentSpikes_OFF = cellSpikeTimesLGN_OFF{intN}(((cellSpikeTimesLGN_OFF{intN} - dblCurT) + dblSynSpikeMem) > 0);
		if ~isempty(vecRecentSpikes_OFF)
			intNumRecentSpikes = numel(vecRecentSpikes_OFF);
			matRecentSpikes_OFF(intN,vecRSPos_OFF(intN):(vecRSPos_OFF(intN)+intNumRecentSpikes-1)) = vecRecentSpikes_OFF;
			vecRSPos_OFF(intN) = vecRSPos_OFF(intN) + intNumRecentSpikes;
			vecRSCounter_OFF(intN) = vecRSCounter_OFF(intN) + intNumRecentSpikes;
		end
	end
	
	% run
	intIter = 0;
	hTic = tic;
	for dblCurT = vecOverallT%(end):dblDeltaT:(vecOverallT(end)+dblDeltaT*100)
		%% check which trial
		intIter = intIter + 1;
		intThisTrial = sum(dblCurT>vecTrialStartSecs);
		
		%% prepare trial
		%start new trial
		if intThisTrial > intPrevTrial
			%% get trial data
			intPrevTrial = intThisTrial;
			intStimType = vecTrialStimType(intThisTrial);
			
			%set counter for this trial
			intTrialT = 0;
			boolStimPresent = true;
		end
		
		%check to stop trial
		intStopTrial = sum(dblCurT>vecStimStopSecs);
		if intStopTrial >= intThisTrial
			boolStimPresent = false;
		end
		
		% check if stimulus should be presented
		intTrialT = intTrialT + 1;
		matRandON = rand(size(matBlankLGN_ON));
		matRandOFF = rand(size(matBlankLGN_ON));
		if boolStimPresent && ((intTrialT <= size(cellLGN_ON{intStimType},3)) || (size(cellLGN_ON{intStimType},3) == 1))
			%rotate template LGN response
			if (size(cellLGN_ON{intStimType},3) == 1),intTrialT = 1;end
			matLGN_ON = cellLGN_ON{intStimType}(:,:,intTrialT);
			matLGN_OFF = cellLGN_OFF{intStimType}(:,:,intTrialT);
			
			%create stochastic spiking
			matThisIterON = matLGN_ON>matRandON;
			matThisIterOFF = matLGN_OFF>matRandOFF;
		else
			matThisIterON = matBlankLGN_ON>matRandON;
			matThisIterOFF = matBlankLGN_OFF>matRandOFF;
		end
		%{
		if boolStimPresent
		subplot(2,2,1)
		imagesc(matLGN_ON-matLGN_OFF);
		subplot(2,2,2)
		imagesc(matRandON);
		subplot(2,2,3)
		imagesc(matThisIterON-matThisIterOFF)
		title(sprintf('Time: %.3fs',dblCurT));
		drawnow;
		pause
		end
		%}
		% assign new LGN spikes to cell arrays
		indSpikingON = matThisIterON(:);
		if any(indSpikingON)
			for intN = find(indSpikingON)'
				%check number of spikes and pre-allocated size
				vecSpikeCounterLGN_ON(intN) = vecSpikeCounterLGN_ON(intN) + 1;
				if vecSpikeCounterPreAllocatedLGN_ON(intN) < vecSpikeCounterLGN_ON(intN)
					vecSpikeCounterPreAllocatedLGN_ON(intN) = vecSpikeCounterPreAllocatedLGN_ON(intN) + intPreAllocationSize;
					cellSpikeTimesLGN_ON{intN} = [cellSpikeTimesLGN_ON{intN} nan(1,intPreAllocationSize)];
				end
				
				% add spike time
				cellSpikeTimesLGN_ON{intN}(vecSpikeCounterLGN_ON(intN)) = dblCurT;
				
				%add recent spikes
				matRecentSpikes_ON(intN,vecRSPos_ON(intN)) = dblCurT;
				vecRSPos_ON(intN) = mod(vecRSPos_ON(intN) + 1,intMaxRecentSpikes);
				vecRSCounter_ON(intN) = vecRSCounter_ON(intN) + 1;
			end
		end
		vecRSPos_ON(vecRSPos_ON==0)=intMaxRecentSpikes; %mod(y,y) is 0, but we want it to be y
		indSpikingOFF = matThisIterOFF(:);
		if any(indSpikingOFF)
			for intN = find(indSpikingOFF)'
				%check number of spikes and pre-allocated size
				vecSpikeCounterLGN_OFF(intN) = vecSpikeCounterLGN_OFF(intN) + 1;
				if vecSpikeCounterPreAllocatedLGN_OFF(intN) < vecSpikeCounterLGN_OFF(intN);
					vecSpikeCounterPreAllocatedLGN_OFF(intN) = vecSpikeCounterPreAllocatedLGN_OFF(intN) + intPreAllocationSize;
					cellSpikeTimesLGN_OFF{intN} = [cellSpikeTimesLGN_OFF{intN} nan(1,intPreAllocationSize)];
				end
				
				% add spike time
				cellSpikeTimesLGN_OFF{intN}(vecSpikeCounterLGN_OFF(intN)) = dblCurT;
				
				%add recent spikes
				matRecentSpikes_OFF(intN,vecRSPos_OFF(intN)) = dblCurT;
				vecRSPos_OFF(intN) = mod(vecRSPos_OFF(intN) + 1,intMaxRecentSpikes);
				vecRSCounter_OFF(intN) = vecRSCounter_OFF(intN) + 1;
			end
		end
		vecRSPos_OFF(vecRSPos_OFF==0)=intMaxRecentSpikes; %mod(y,y) is 0, but we want it to be y
		
		%%% remove old spikes
		for intN = find(vecRSCounter_OFF>0)
			intPos = getCycPos(vecRSPos_OFF(intN),vecRSCounter_OFF(intN),intMaxRecentSpikes);
			if dblCurT - matRecentSpikes_OFF(intN,intPos) > dblSynSpikeMem
				matRecentSpikes_OFF(intN,intPos) = nan;
				vecRSCounter_OFF(intN) = vecRSCounter_OFF(intN) - 1;
			end
		end
		for intN = find(vecRSCounter_ON>0)
			intPos = getCycPos(vecRSPos_ON(intN),vecRSCounter_ON(intN),intMaxRecentSpikes);
			if dblCurT - matRecentSpikes_ON(intN,intPos) > dblSynSpikeMem
				matRecentSpikes_ON(intN,intPos) = nan;
				vecRSCounter_ON(intN) = vecRSCounter_ON(intN) - 1;
			end
		end
		
		%% calculate all excitatory currents due to spikes from LGN
		intSubType = 1;
		vecLGNCurrents = zeros(size(vecThisV));
		for intField=[1 2]
			if intField == 1
				vecSW = vecSynWeightON_to_Cort; %#ok<NASGU>
				vecSC = vecSynConductanceON_to_Cort;
				vecSD = vecSynDelayON_to_Cort;
				matC = matSynConnON_to_Cort;
				cellSynTarg = cellSynTargetsLGN_ON;
				indMoot = indMootCells_ON;
				matRS = matRecentSpikes_ON;
				vecRSPosLGN = vecRSPos_ON;
				vecRSCounterLGN = vecRSCounter_ON;
			else
				vecSW = vecSynWeightOFF_to_Cort; %#ok<NASGU>
				vecSC = vecSynConductanceOFF_to_Cort;
				vecSD = vecSynDelayOFF_to_Cort;
				matC = matSynConnOFF_to_Cort;
				cellSynTarg = cellSynTargetsLGN_OFF;
				indMoot = indMootCells_OFF;
				matRS = matRecentSpikes_OFF;
				vecRSPosLGN = vecRSPos_OFF;
				vecRSCounterLGN = vecRSCounter_OFF;
			end
			
			%loop
			vecPSPs = zeros(size(matC,1),1);
			for intNeuron=find(vecRSCounterLGN>0 & ~indMoot)
				%get recent spikes
				vecPos = getCycPosVec(vecRSPosLGN(intNeuron),vecRSCounterLGN(intNeuron),intMaxRecentSpikes);
				vecSpikeTimes = matRS(intNeuron,vecPos);
				vecTargetSynapses = cellSynTarg{intNeuron};
				
				%get neuron based properties
				dblDelay =  vecSD(vecTargetSynapses(1));
				dblTauPeak = vecTauPeakByType(intSubType);
				
				%get synapse data
				dblWeight = 1;
				
				%calculate PSPs
				vecPSPs(vecTargetSynapses) = dblWeight .* vecSC(vecTargetSynapses) .* fPSP(dblCurT-dblDelay,vecSpikeTimes,dblTauPeak);
			end
			
			% sum all transmissions incoming to each neuron
			indNonZero = vecPSPs~=0;
			if any(indNonZero)
				vecLGNCurrents = vecLGNCurrents + accumarray(matC(indNonZero,2),vecPSPs(indNonZero),[intCortexCells 1],@sum,0);
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
			%{
			vecAllTargetSynapses = cellSynTargets{intNeuron};
			vecTargetSynTypes = vecSynType(vecAllTargetSynapses);
			for intSynSubType = unique(vecTargetSynTypes(:))'
				%get subtypes
				vecTargetSynapses = vecAllTargetSynapses(vecTargetSynTypes==intSynSubType);
				
				%get neuron based properties
				dblDelay =  vecSynDelay(vecTargetSynapses(1));
				dblTauPeak = vecCortTauSyn(vecTargetSynapses(1));
				
				%get synapse data
				dblWeight = 1;
				
				%calculate PSPs
				vecPSPs(vecTargetSynapses) = dblWeight .* vecSynConductance(vecTargetSynapses) .* fPSP(dblCurT-dblDelay,vecSpikeTimes,dblTauPeak);
				%calculate PSPs
				%vecPSPs(vecTargetSynapses) = dblWeight .* vecSynConductance(vecTargetSynapses) .* fPSP2(bsxfun(@minus,dblCurT-vecSynDelay(vecTargetSynapses),vecSpikeTimes),dblTauPeak);
			end
			%}
			vecTargetSynapses = cellSynTargets{intNeuron};
			if isempty(vecTargetSynapses),continue;end
			%get neuron based properties
			dblTauPeak = vecCortTauSyn(vecTargetSynapses(1));
			
			%get synapse data
			dblWeight = 1;
			
			%calculate PSPs
			%vecPSPs(vecTargetSynapses) = vecSynWeight(vecTargetSynapses) .* vecSynConductance(vecTargetSynapses) .* fPSP2(bsxfun(@minus,dblCurT-vecSynDelay(vecTargetSynapses),vecSpikeTimes),dblTauPeak);
			vecPSPs(vecTargetSynapses) = dblWeight .* vecSynConductance(vecTargetSynapses) .* fPSP2(bsxfun(@minus,dblCurT-vecSynDelay(vecTargetSynapses),vecSpikeTimes),dblTauPeak);
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
			-vecInhCurrents.*(vecPrevV-vecCellV_I) ...
			-vecCellG_Leak.*(vecPrevV-vecCellV_Leak) ...
			-vecAHPs.*(vecPrevV-vecCellV_AHP) ...
			) ./ vecCellCm;
		if any(isnan(vecThisV) | isinf(vecThisV));
			intErrN = find(isnan(vecThisV) | isinf(vecThisV),1);
			if boolClust
				printf('\n.. ERROR DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
			else
				fprintf('\n.. ERROR DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
				warning([mfilename ':VmOutOfBounds'],'Membrane voltage out of bounds!');
			end
			break;
		end
		%{
		subplot(2,2,1)
		plot(vecThisV);
		title('Membrane voltage')
		subplot(2,2,2)
		plot(vecSpikeCounterCortex)
		title('Spikes Cortex')
		subplot(2,2,3)
		plot(vecSpikeCounterLGN_ON)
		title('Spikes LGN ON')
		subplot(2,2,4)
		plot(vecSpikeCounterLGN_OFF)
		title('Spikes LGN OFF')
		drawnow;
		%}
		
		% spiking
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
				
				%add recent spikes
				matRecentSpikes(intN,vecRSPos(intN)) = dblCurT;
				vecRSPos(intN) = mod(vecRSPos(intN) + 1,intMaxRecentSpikes);
				vecRSCounter(intN) = vecRSCounter(intN) + 1;
			end
		end
		vecRSPos(vecRSPos==0)=intMaxRecentSpikes; %mod(y,y) is 0, but we want it to be y
		
		%% message
		if mod(intIter,200)==0
			dblFracDone = intIter/numel(vecOverallT);
			dblTimeElapsed = toc(hTic);
			dblTimeRemaining = (dblTimeElapsed/dblFracDone)-dblTimeElapsed;
			strTimeDone = datestr(datenum(datetime('now'))+dblTimeRemaining/(24*60*60));
			
			vecTime = fix(clock);
			strTime = sprintf('%02d:%02d:%02d',vecTime(4:6));
			if boolClust
				printf('Now at t=%.1fs / %.1fs [%s]; estimated completion: [%s]\n',dblCurT,vecOverallT(end),strTime,strTimeDone);
			else
				fprintf('Now at t=%.1fs / %.1fs [%s]; estimated completion: [%s]\n',dblCurT,vecOverallT(end),strTime,strTimeDone);
			end
		end
	end
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
	sData.cellLGN_ON = cellLGN_ON;
	sData.cellLGN_OFF = cellLGN_OFF;
	sData.vecTrialOris = vecTrialOris;
	sData.vecTrialOriIdx = vecTrialStimType;
	sData.vecStimStartSecs = vecStimStartSecs;
	sData.vecStimStopSecs = vecStimStopSecs;
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
end
