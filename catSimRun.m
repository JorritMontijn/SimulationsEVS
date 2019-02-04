function sSimRunAgg = catSimRun(sSimRunAgg,sSimRun,dblLastT,dblLastRep,boolIncludeLGN)
	if nargin < 5
		boolIncludeLGN = true;
	end
	if isempty(sSimRunAgg)
		if nargin < 3
			sSimRunAgg = sSimRun;
		else
			%timing variables
			sSimRunAgg.vecOverallT = sSimRun.vecOverallT + dblLastT;
			sSimRunAgg.vecTrialStartSecs = sSimRun.vecTrialStartSecs + dblLastT;
			sSimRunAgg.vecStimStartSecs = sSimRun.vecStimStartSecs + dblLastT;
			sSimRunAgg.vecStimStopSecs = sSimRun.vecStimStopSecs + dblLastT;
			sSimRunAgg.vecTrialEndSecs = sSimRun.vecTrialEndSecs + dblLastT;
			sSimRunAgg.vecTrialStimRep = sSimRun.vecTrialStimRep + dblLastRep;
			
			%other variables
			cellIgnore = {'vecTrialStartSecs','vecStimStartSecs','vecStimStopSecs','vecTrialEndSecs','vecTrialStimRep'};
			cellFields = fieldnames(sSimRun);
			for intField=1:numel(cellFields)
				strField = cellFields{intField};
				if ~isempty(strfind(strField,'vecTrial')) && ~ismember(strField,cellIgnore)
					sSimRunAgg.(strField) = sSimRun.(strField);
				end
			end
			
			%cortex
			for i=1:numel(sSimRun.cellSpikeTimesCortex)
				sSimRunAgg.cellSpikeTimesCortex{i} = sSimRun.cellSpikeTimesCortex{i} + dblLastT;
			end
			if boolIncludeLGN
				%LGN ON
				for i=1:numel(sSimRun.cellSpikeTimesLGN_ON)
					sSimRunAgg.cellSpikeTimesLGN_ON{i} = sSimRun.cellSpikeTimesLGN_ON{i} + dblLastT;
				end
				
				%LGN OFF
				for i=1:numel(sSimRun.cellSpikeTimesLGN_OFF)
					sSimRunAgg.cellSpikeTimesLGN_OFF{i} = sSimRun.cellSpikeTimesLGN_OFF{i} + dblLastT;
				end
			end
		end
	else
		if nargin < 3 || isempty(dblLastT)
			dblLastT = sSimRunAgg.vecOverallT(end);
		end
		if nargin < 4 || isempty(dblLastRep)
			dblLastRep = max(sSimRunAgg.vecTrialStimRep);
		end
		
		%timing variables
		sSimRunAgg.vecOverallT = cat(2,sSimRunAgg.vecOverallT,sSimRun.vecOverallT + dblLastT);
		sSimRunAgg.vecTrialStartSecs = cat(2,sSimRunAgg.vecTrialStartSecs,sSimRun.vecTrialStartSecs + dblLastT);
		sSimRunAgg.vecStimStartSecs = cat(2,sSimRunAgg.vecStimStartSecs,sSimRun.vecStimStartSecs + dblLastT);
		sSimRunAgg.vecStimStopSecs = cat(2,sSimRunAgg.vecStimStopSecs,sSimRun.vecStimStopSecs + dblLastT);
		sSimRunAgg.vecTrialEndSecs = cat(2,sSimRunAgg.vecTrialEndSecs,sSimRun.vecTrialEndSecs + dblLastT);
		sSimRunAgg.vecTrialStimRep = cat(2,sSimRunAgg.vecTrialStimRep,sSimRun.vecTrialStimRep+dblLastRep);
		
		%other variables
		cellIgnore = {'vecTrialStartSecs','vecStimStartSecs','vecStimStopSecs','vecTrialEndSecs','vecTrialStimRep'};
		cellFields = fieldnames(sSimRun);
		for intField=1:numel(cellFields)
			strField = cellFields{intField};
			if ~isempty(strfind(strField,'vecTrial')) && ~ismember(strField,cellIgnore)
				sSimRunAgg.(strField) = cat(2,sSimRunAgg.(strField),sSimRun.(strField));
			end
		end
		
		
		%cortex
		for i=1:numel(sSimRun.cellSpikeTimesCortex)
			sSimRunAgg.cellSpikeTimesCortex{i} = cat(2,sSimRunAgg.cellSpikeTimesCortex{i},sSimRun.cellSpikeTimesCortex{i} + dblLastT);
		end
		
		if boolIncludeLGN
			%LGN ON
			for i=1:numel(sSimRun.cellSpikeTimesLGN_ON)
				sSimRunAgg.cellSpikeTimesLGN_ON{i} = cat(2,sSimRunAgg.cellSpikeTimesLGN_ON{i},sSimRun.cellSpikeTimesLGN_ON{i} + dblLastT);
			end
			
			%LGN OFF
			for i=1:numel(sSimRun.cellSpikeTimesLGN_OFF)
				sSimRunAgg.cellSpikeTimesLGN_OFF{i} = cat(2,sSimRunAgg.cellSpikeTimesLGN_OFF{i},sSimRun.cellSpikeTimesLGN_OFF{i} + dblLastT);
			end
		end
	end
end
