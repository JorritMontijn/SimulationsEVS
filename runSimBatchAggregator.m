%function runSimBatchAggregator

%% define directories
clearvars;
strSourceDir = 'D:\Simulations\SimResults\';
strTargetDir = 'D:\Data\Processed\V1_LIFmodel\';

%% load directory contents
sFiles = dir(strSourceDir);
intFiles = numel(sFiles);
strProcessingID = '';
vecProcessed = false(1,intFiles);
%% loop through files
for intFile = 1:intFiles
	strFile = sFiles(intFile).name;
	if length(strFile) > 4 && strcmpi(strFile(end-3:end),'.mat') && ~vecProcessed(intFile)
		fprintf('Checking file %s (%d/%d)\n',strFile,intFile,intFiles);
		sLoad = load([strSourceDir strFile]);
		if ((isfield(sLoad,'sSimRun') && ~isempty(strProcessingID)) || (isfield(sLoad,'sParams') && isfield(sLoad.sParams,'strConnFile')))
			if isfield(sLoad,'sParams') && isfield(sLoad.sParams,'strCurrID')
				strCurrID = sLoad.sParams.strConnFile;
			end
			if isempty(strProcessingID)
				%% first file; assign all persistent variables
				sParams = sLoad.sParams;
				sSimRun = sLoad.sSimRun;
				clear sLoad;
				strProcessingID = strCurrID;
				
				%msg
				fprintf('\b; Processing connectivity [%s]; file <%s> [%s]\n',strProcessingID,strFile,getTime);
				
				%rem fields
				sData = rmfield(sParams,{'sStimParams','sStimInputs','sConnParams','sConnectivity'});
				
				%suppress duplicate warning, concat struct, & re-enable
				warning('off','catstruct:DuplicatesFound');
				sData = catstruct(sData,sSimRun(end),sParams.sStimParams,sParams.sStimInputs,sParams.sConnParams,sParams.sConnectivity);
				warning('on','catstruct:DuplicatesFound');
				
				%remove last repetition and add to aggregate
				sSimRunAgg = remSimLastRep(sSimRun(end));
				
				%remove to-be-overwritten fields from sData
				sData = rmfield(sData,fieldnames(sSimRunAgg));
				
				%loop through SimRuns if multiple exist
				for intRun=1:(numel(sSimRun)-1)
					%skip if empty
					if isempty(sSimRun(intRun).vecOverallT);continue;end
					
					% add spiking arrays to existing structure
					fprintf('   Adding %d trials (total: %d); [%d/%d]/<%s> of [%s]; [%s]\n',numel(sSimRun(intRun).vecTrialStimType),numel(sSimRun(intRun).vecTrialStimType)+numel(sSimRunAgg.vecTrialStimType),...
						intRun,numel(sSimRun),strFile,strProcessingID,getTime);
					
					%remove last repetition
					sSimRun(intRun) = remSimLastRep(sSimRun(intRun));
					
					%add spiking arrays to aggregate
					sSimRunAgg = catSimRun(sSimRunAgg,sSimRun(intRun));
				end
				
				%set flag
				vecProcessed(intFile) = true;
			elseif strcmpi(strProcessingID,strCurrID)
				%% extract
				if ~isfield(sLoad,'sSimRun')
					fprintf('   Warning: file <%s> of [%s] has no sSimRun structure; [%s]\n',strFile,strProcessingID,getTime);
					continue;
				end
				clearvars sParams sSimRun;
				sSimRun = sLoad.sSimRun;
				clearvars sLoad;
				
				%loop through SimRuns if multiple exist
				for intRun=1:numel(sSimRun)
					%skip if empty
					if isempty(sSimRun(intRun).vecOverallT);continue;end
					if isempty(sSimRunAgg)
						intNumPrevTrials=0;
					else
						intNumPrevTrials=numel(sSimRunAgg.vecTrialStimType);
						dblLastT = sSimRunAgg.vecOverallT(end);
						dblLastRep = max(sSimRunAgg.vecTrialStimRep);
					end
					% add spiking arrays to existing structure
					fprintf('   Adding %d trials (total: %d); run [%d/%d]/<%s> of [%s]; [%s]\n',numel(sSimRun(intRun).vecTrialStimType),numel(sSimRun(intRun).vecTrialStimType)+intNumPrevTrials,intRun,numel(sSimRun),strFile,strProcessingID,getTime);
					
					%remove last repetition
					sSimRun(intRun) = remSimLastRep(sSimRun(intRun));
					
		
					%add spiking arrays to aggregate
					%sSimRunAgg = catSimRun(sSimRunAgg,sSimRun(intRun));
					sSimRunAgg = catSimRun(sSimRunAgg,sSimRun,dblLastT,dblLastRep);
					
					%sSimRunAgg = catSimRun(sSimRunAgg2,sSimRunAgg,0,0);
				end
				
				%set flag
				vecProcessed(intFile) = true;
			else
				warning([mfilename ':MultipleDataFiles'],'Data from multiple connectivity structures detected; will only process <%s>',strProcessingID);
			end
		end
	end
end

%% transform cellSpikes to uint32 arrays
dblStepSize = (sSimRunAgg.vecOverallT(end)-sSimRunAgg.vecOverallT(1))/(numel(sSimRunAgg.vecOverallT)-1);
vecUniqueSteps = unique(diff(sSimRunAgg.vecOverallT));
vecUniqueOffsets = abs(vecUniqueSteps-dblStepSize);
if any(vecUniqueOffsets > (dblStepSize/1000))
	warning([mfilename ':StepSizeError'],'Step size offset is larger than acceptable tolerance: %f',max(vecUniqueOffsets));
end
fprintf('Transforming Doubles to Int32 timestamps; maximum step-error is %e with step size %e [%s]\n',max(vecUniqueOffsets),dblStepSize,getTime);
sSimRunAgg.cellSpikeTimesCortex = doSpikeDoubleToInt(sSimRunAgg.cellSpikeTimesCortex,dblStepSize,sSimRunAgg.vecOverallT(1));
fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
sSimRunAgg.cellSpikeTimesLGN_ON = doSpikeDoubleToInt(sSimRunAgg.cellSpikeTimesLGN_ON,dblStepSize,sSimRunAgg.vecOverallT(1));
fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
sSimRunAgg.cellSpikeTimesLGN_OFF = doSpikeDoubleToInt(sSimRunAgg.cellSpikeTimesLGN_OFF,dblStepSize,sSimRunAgg.vecOverallT(1));
fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);

%% add sim run aggregate to sData
% add spiking arrays to existing structure
fprintf('Catstructing sData with SimRunAgg [%s]\n',getTime);
sData = catstruct(sData,sSimRunAgg);

%% save data
strOutputFile = sprintf('Simulation_xAreaDistributed_%d_%s.mat',sum(cat(2,cast(getDate,'double')*1000,cast(getTime,'double'))),getDate);
fprintf('Aggregation complete, saving data to <%s%s> [%s]\n',strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7');
fprintf('Done! [%s]\n',getTime);
