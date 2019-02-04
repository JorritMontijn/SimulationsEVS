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
function sData = getSimRunNoStim(sData)
	%check if we're on cluster
	global boolClust;
	if isempty(boolClust),boolClust=false;end
	hTic = tic;
	
	%unpack
	intCortexCells = sData.intCortexCells;
	intCellsV1 = sData.intCellsV1;
	intCellsV2 = sData.intCellsV2;
	dblDeltaT = sData.dblDeltaT;
	vecOverallT = sData.vecOverallT;
	
	matSynFromTo = sData.matSynFromTo;
	dblSynSpikeMem = sData.dblSynSpikeMem;
	vecSynExcInh = sData.vecSynExcInh;
	vecSynDelay = sData.vecSynDelay;
	vecSynConductance = sData.vecSynConductance;
	vecSynWeight = sData.vecSynWeight;
	vecSynType = sData.vecSynType;
	
	vecCellThresh = sData.vecCellThresh;
	vecCellTypes = sData.vecCellTypes;
	vecCellTauPeak = sData.vecCellTauPeak;
	vecTauPeakByType=nan(1,numel(unique(vecCellTypes)));
	vecTauPeakByType(vecCellTypes) = vecCellTauPeak;
	vecCellV_E = sData.vecCellV_E;
	vecCellV_I = sData.vecCellV_I;
	vecCellV_AHP = sData.vecCellV_AHP;
	vecCellV_Leak = sData.vecCellV_Leak;
	vecCellCm = sData.vecCellCm;
	vecCellG_Leak = sData.vecCellG_Leak;
	vecCellG_AHP = sData.vecCellG_AHP;
	%vecCellRefracT = max(dblDeltaT,sData.vecCellRefracT);
	
	vecThisV = sData.vecThisV;
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	vecSpikeCounterCortex = sData.vecSpikeCounterCortex;
	intPreAllocationSize = sData.intPreAllocationSize;
	if isfield(sData,'boolSaveVm'),boolSaveVm = sData.boolSaveVm;else boolSaveVm = false;end
	
	%get input
	vecInputIdx = sData.vecInputIdx; %[1 x T] with M index values
	matInput = sData.matInput; %[N x M]; neurons by input indices
	
	%define vars
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
	%vecTimeSinceSpike = (vecOverallT(1) + 10)*ones(size(vecThisV));
	
	%% pre-allocate Vm matrix
	if boolSaveVm
		matVm = nan(intCortexCells,numel(vecOverallT));
	end
	
	% run
	intIter = 0;
	hTic = tic;
	%fprintf('\n   >>> Starting simulation run [%s]\n',getTime);
	dblLastMsg = toc(hTic);
	for dblCurT = vecOverallT%(end):dblDeltaT:(vecOverallT(end)+dblDeltaT*100)
		intIter = intIter + 1;
		
		%% get inputs
		intInputIdx = vecInputIdx(intIter); %[1 x T] with M index values
		vecInputG = matInput(:,intInputIdx); %[N x M]; neurons by input indices
	
		
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
			dblWeight = 1;
			
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
			vecInputG ...
			-vecExcCurrents.*(vecPrevV-vecCellV_E) ...
			-vecInhCurrents.*(vecPrevV-vecCellV_I) ...
			-vecCellG_Leak.*(vecPrevV-vecCellV_Leak) ...
			-vecAHPs.*(vecPrevV-vecCellV_AHP) ...
			) ./ vecCellCm;
		
		if any(isnan(vecThisV) | isinf(vecThisV));
			intErrN = find(isnan(vecThisV) | isinf(vecThisV));
			if boolClust
				printf('\n.. ERROR DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
			else
				fprintf('\n.. ERROR DETECTED! Please check; neuron %d; Currents; LGN: %f, E: %f, I: %f, L: %f, H: %f [%s]\n',intErrN,vecLGNCurrents(intErrN),vecExcCurrents(intErrN),vecInhCurrents(intErrN),vecCellG_Leak(intErrN),vecAHPs(intErrN),getTime);
				warning([mfilename ':VmOutOfBounds'],'Membrane voltage out of bounds!');
			end
			break;
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
		
		
		%save Vm
		if boolSaveVm 
			matVm(:,intIter) = vecThisV;
		end
		
		%% message
		if (toc(hTic) - dblLastMsg) > 5
			dblLastMsg = toc(hTic);
			dblFracDone = intIter/numel(vecOverallT);
			dblTimeElapsed = toc(hTic);
			dblTimeRemaining = (dblTimeElapsed/dblFracDone)-dblTimeElapsed;
			strTimeDone = datestr(datenum(datetime('now'))+dblTimeRemaining/(24*60*60));
			
			vecTime = fix(clock);
			strTime = sprintf('%02d:%02d:%02d',vecTime(4:6));
			
			%fprintf('Now at t=%.1fs / %.1fs [%s]; estimated completion: [%s]\n',dblCurT,vecOverallT(end),strTime,strTimeDone);
		end
	end
	
	%end msg
	hToc = toc(hTic);
	if boolClust
		printf('Done! Total run took %f seconds [%s]\n',hToc,getTime);
	else
		%fprintf('Done! Total run took %f seconds [%s]\n',hToc,getTime);
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
	
	sData.vecThisV = vecThisV;
	
	sData.intIter = intIter;
	
	sData.cellSpikeTimesCortex = cellSpikeTimesCortex;
	
	sData.vecSpikeCounterCortex = vecSpikeCounterCortex;
	sData.intPreAllocationSize = intPreAllocationSize;
	if boolSaveVm,sData.matVm = matVm;end
end
