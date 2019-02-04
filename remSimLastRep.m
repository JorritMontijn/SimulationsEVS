function sSimRun = remSimLastRep(sSimRun)
	%get max time
	dblMaxTimeSpikes = max(cellfun(@max, sSimRun.cellSpikeTimesCortex));
	dblMaxTimeSim = max(sSimRun.vecOverallT);
	dblMaxTime = min(dblMaxTimeSim,dblMaxTimeSpikes);
	
	%get rem trials
	indRemTrials = sSimRun.vecTrialStimRep >= min(sSimRun.vecTrialStimRep(~(sSimRun.vecTrialEndSecs <= dblMaxTime)));
	if any(~indRemTrials)
		dblMaxTime = sSimRun.vecTrialEndSecs(find(~indRemTrials,1,'last'));
	else
		dblMaxTime = 0;
	end
	
	%remove stim times
	sSimRun.vecStimStartSecs(indRemTrials) = [];
	sSimRun.vecStimStopSecs(indRemTrials) = [];
	
	%trial + time vars
	sSimRun.vecOverallT(sSimRun.vecOverallT>dblMaxTime) = [];
	cellFields = fieldnames(sSimRun);
	for intField=1:numel(cellFields)
		strField = cellFields{intField};
		if ~isempty(strfind(strField,'vecTrial'))
			sSimRun.(strField)(indRemTrials) = [];
		end
	end
	
	
	%remove spikes
	%cortex
	cellSpikeTimes = sSimRun.cellSpikeTimesCortex;
	for i=1:numel(cellSpikeTimes)
		vecSpikeTimes = cellSpikeTimes{i};
		vecSpikeTimes(vecSpikeTimes > dblMaxTime) = [];
		cellSpikeTimes{i} = vecSpikeTimes;
	end
	sSimRun.cellSpikeTimesCortex = cellSpikeTimes;
	
	%LGN ON
	cellSpikeTimes = sSimRun.cellSpikeTimesLGN_ON;
	for i=1:numel(cellSpikeTimes)
		vecSpikeTimes = cellSpikeTimes{i};
		vecSpikeTimes(vecSpikeTimes > dblMaxTime) = [];
		cellSpikeTimes{i} = vecSpikeTimes;
	end
	sSimRun.cellSpikeTimesLGN_ON = cellSpikeTimes;
	
	%LGN OFF
	cellSpikeTimes = sSimRun.cellSpikeTimesLGN_OFF;
	for i=1:numel(cellSpikeTimes)
		vecSpikeTimes = cellSpikeTimes{i};
		vecSpikeTimes(vecSpikeTimes > dblMaxTime) = [];
		cellSpikeTimes{i} = vecSpikeTimes;
	end
	sSimRun.cellSpikeTimesLGN_OFF = cellSpikeTimes;
end
