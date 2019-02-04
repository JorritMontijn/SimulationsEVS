%function runSimBatchAggregator

%% define directories
clearvars;
strSourceDir = 'A:\SimResults\';
%strSourceDir = 'B:\DataBackup\SimResults\ExcOnly\';
strTargetDir = 'A:\SimAggregates\';

boolIncludeLGN = false;

%% load directory contents and prep vars
sFiles = dir(strSourceDir);
intFiles = numel(sFiles);
strProcessingID = '';
vecProcessed = false(1,intFiles);
intBufferSize = 20;
intCurBuf = 0;
sSimRunAgg = [];
sSimRunBuf = [];
intNumPrevTrials=0;
dblLastT = 0;
dblLastRep = 0;

%% loop through files
for intFile = [intFiles 1:(intFiles-1)]
	strFile = sFiles(intFile).name;
	if length(strFile) > 4 && strcmpi(strFile(end-3:end),'.mat') && ~vecProcessed(intFile)
		fprintf('Checking file %s (%d/%d)\n',strFile,intFile,intFiles);
		sLoad = load([strSourceDir strFile]);
		if ((isfield(sLoad,'sSimRun') && ~isempty(strProcessingID)) || (isfield(sLoad,'sParams') && isfield(sLoad.sParams,'strConnFile')))
			%if isfield(sLoad,'sParams') && isfield(sLoad.sParams,'strConnFile')
				strCurrID = sLoad.sParams.strConnFile;
			%end
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
				
				if isempty(sSimRun(end).vecOverallT)
					%add general info
					sData = catstruct(sData,sParams.sStimInputs,sParams.sConnectivity,sParams.sStimParams,sParams.sConnParams);
				else
					%suppress duplicate warning, concat struct, & re-enable
					warning('off','catstruct:DuplicatesFound');
					sData = catstruct(sData,sSimRun(end),sParams.sStimParams,sParams.sStimInputs,sParams.sConnParams,sParams.sConnectivity);
					warning('on','catstruct:DuplicatesFound');
					
					%remove last repetition and add to aggregate
					sSimRunBuf = remSimLastRep(sSimRun(end));
					intCurBuf = intCurBuf + 1;
					
					%remove to-be-overwritten fields from sData
					sData = rmfield(sData,fieldnames(sSimRunBuf));
					
					%loop through SimRuns if multiple exist
					for intRun=1:(numel(sSimRun)-1)
						%skip if empty
						if isempty(sSimRun(intRun).vecOverallT);continue;end
						
						%remove last repetition
						sSimRun(intRun) = remSimLastRep(sSimRun(intRun));
						
						% add spiking arrays to existing structure
						fprintf('   Adding %d trials (total: %d); [%d/%d]/<%s> of [%s]; [%s]\n',numel(sSimRun(intRun).vecTrialStimType),numel(sSimRun(intRun).vecTrialStimType)+numel(sSimRunBuf.vecTrialStimType),...
							intRun,numel(sSimRun),strFile,strProcessingID,getTime);
						
						%add spiking arrays to aggregate
						sSimRunBuf = catSimRun(sSimRunBuf,sSimRun(intRun),0,0,boolIncludeLGN);
						intCurBuf = intCurBuf + 1;
					end
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
					
					%remove last repetition
					sSimRun(intRun) = remSimLastRep(sSimRun(intRun));
					intCurTrials = numel(sSimRun(intRun).vecTrialStimType);
					dblCurT = sSimRun(intRun).vecOverallT(end);
					dblCurRep = max(sSimRun(intRun).vecTrialStimRep);
					
					% add spiking arrays to existing structure
					fprintf('   Adding %d trials (total: %d); run [%d/%d]/<%s> of [%s]; [%s]\n',numel(sSimRun(intRun).vecTrialStimType),numel(sSimRun(intRun).vecTrialStimType)+intNumPrevTrials,intRun,numel(sSimRun),strFile,strProcessingID,getTime);
					
					%add spiking arrays to aggregate
					sSimRunBuf = catSimRun(sSimRunBuf,sSimRun(intRun),dblLastT,dblLastRep,boolIncludeLGN);
					intCurBuf = intCurBuf + 1;
					
					%increment counters
					intNumPrevTrials=intNumPrevTrials + intCurTrials;
					dblLastT = dblLastT + dblCurT;
					dblLastRep = dblLastRep + dblCurRep;
				end
				
				%%
				%check buffer
				if intCurBuf >= intBufferSize
					%transform buffer
					fprintf('Buffer is at maximum size; adding to aggregate [%s]\n',getTime);
					
					dblStepSize = (sSimRunBuf.vecOverallT(end)-sSimRunBuf.vecOverallT(1))/(numel(sSimRunBuf.vecOverallT)-1);
					vecUniqueSteps = unique(diff(sSimRunBuf.vecOverallT));
					vecUniqueOffsets = abs(vecUniqueSteps-dblStepSize);
					if any(vecUniqueOffsets > (dblStepSize/1000))
						warning([mfilename ':StepSizeError'],'Step size offset is larger than acceptable tolerance: %f',max(vecUniqueOffsets));
					end
					fprintf('Transforming Doubles to Int32 timestamps; maximum step-error is %e with step size %e [%s]\n',max(vecUniqueOffsets),dblStepSize,getTime);
					sSimRunBuf.cellSpikeTimesCortex = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesCortex,dblStepSize,0);
					fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
					if boolIncludeLGN
						sSimRunBuf.cellSpikeTimesLGN_ON = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesLGN_ON,dblStepSize,0);
						fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
						sSimRunBuf.cellSpikeTimesLGN_OFF = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesLGN_OFF,dblStepSize,0);
						fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);
					end
					
					%add buffer to aggregate
					sSimRunAgg = catSimRun(sSimRunAgg,sSimRunBuf,0,0,boolIncludeLGN);
					
					%reset buffer
					intCurBuf = 0;
					sSimRunBuf = [];
				end
				
				%set flag
				vecProcessed(intFile) = true;
			else
				warning([mfilename ':MultipleDataFiles'],'Data from multiple connectivity structures detected; will only process <%s>',strProcessingID);
			end
		end
	end
end
%% add last buffer
if intCurBuf > 0
	fprintf('Adding final buffer to aggregate [%s]\n',getTime);

	dblStepSize = (sSimRunBuf.vecOverallT(end)-sSimRunBuf.vecOverallT(1))/(numel(sSimRunBuf.vecOverallT)-1);
	vecUniqueSteps = unique(diff(sSimRunBuf.vecOverallT));
	vecUniqueOffsets = abs(vecUniqueSteps-dblStepSize);
	if any(vecUniqueOffsets > (dblStepSize/1000))
		warning([mfilename ':StepSizeError'],'Step size offset is larger than acceptable tolerance: %f',max(vecUniqueOffsets));
	end
	fprintf('Transforming Doubles to Int32 timestamps; maximum step-error is %e with step size %e [%s]\n',max(vecUniqueOffsets),dblStepSize,getTime);
	sSimRunBuf.cellSpikeTimesCortex = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesCortex,dblStepSize,0);
	fprintf('   Cortex spikes have been transformed... [%s]\n',getTime);
	
	if boolIncludeLGN
		sSimRunBuf.cellSpikeTimesLGN_ON = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesLGN_ON,dblStepSize,0);
		fprintf('   LGN ON spikes have been transformed... [%s]\n',getTime);
		sSimRunBuf.cellSpikeTimesLGN_OFF = doSpikeDoubleToInt(sSimRunBuf.cellSpikeTimesLGN_OFF,dblStepSize,0);
		fprintf('   LGN OFF spikes have been transformed... [%s]\n',getTime);
	end
	%add buffer to aggregate
	sSimRunAgg = catSimRun(sSimRunAgg,sSimRunBuf,0,0,boolIncludeLGN);

	%reset buffer
	intCurBuf = 0;
	sSimRunBuf = [];
end
%% add sim run aggregate to sData
% add spiking arrays to existing structure
fprintf('Catstructing sData with SimRunAgg [%s]\n',getTime);
sData = catstruct(sData,sSimRunAgg);

%% save data
strOutputFile = sprintf('Simulation_xAreaDistributed_%d_%s.mat',sum(cat(2,cast(getDate,'double')*1000,cast(getTime,'double'))),getDate);
fprintf('Aggregation complete, saving data to <%s%s> [%s]\n',strTargetDir,strOutputFile,getTime);
save([strTargetDir strOutputFile],'-struct','sData','-v7.3');
fprintf('Done! [%s]\n',getTime);
