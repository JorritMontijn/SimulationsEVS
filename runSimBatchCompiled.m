%just to be sure
function sData = runSimBatchCompiled(strInput)
	hTic = tic;
	%% split input
	cellIn = strsplit(strInput,',');
	strRunningTime = cellIn{1};
	strConnFile = cellIn{2};
	strStimFile = cellIn{3};
	intWorker = str2double(cellIn{4});
	intTag = 5;
	dblAttention = 0;
	%% msg
	printf('Starting simulation; input was <%s> (%d args); [%s]\n',strInput,numel(cellIn),getTime);
	%% parse
	if numel(cellIn)>4
		strAtt = cellIn{5};
		cellAtt = strsplit(strAtt,'=');
		if numel(cellAtt) > 1
			boolAtt = false;
			if strcmpi(cellAtt{1},'att')
				cellAttSub = strsplit(cellAtt{2},'_');
				dblAttention = str2double(cellAttSub{1});intAttArea=[];
				if ~isempty(dblAttention),boolAtt=true;end
				if numel(cellAttSub)>1;intAttArea=str2double(cellAttSub{2});end
				if isempty(intAttArea) || ~(intAttArea > 0),intAttArea=0;end
			end
			if boolAtt
				printf('Attention requested; setting factor to <%.2f>, applying to area %d (%s)\n',dblAttention,intAttArea,getTime);
				intTag=intTag+1;
			else
				printf('Supplied argument not recognized! Input was: (%s)\n',strInput);
			end
		end
	end
	if numel(cellIn) >= intTag
		strTag = cellIn{intTag};
	else
		strTag = '';
	end
	
	%% set log path
	global strLogDir;
	global boolClust;
	global boolSaveVm;
	boolSaveVm = false;
	
	%% running on cluster?
	if ispc
		%add directories to paths
		strHome = mfilename('fullpath');
		strHome = strHome(1:(end-length(mfilename)));
		strLogDir = [strHome 'logs' filesep];
		strOutputDir = ['A:\SimResults' filesep];
		strStimDir = [strHome 'Stimulation' filesep];
		strConnDir = [strHome 'Connectivity' filesep];
	else
		%add directories to paths
		strHome = [filesep 'home' filesep 'montijn' filesep];
		strLogDir = [strHome 'logs' filesep];
		strOutputDir = [strHome 'SimResults' filesep];
		strStimDir = [strHome 'Stimulation' filesep];
		strConnDir = [strHome 'Connectivity' filesep];
		boolSaveVm = false; %always disable Vm saving on cluster
	end
	boolClust = false;
	
	%% transform running time to double
	if isnumeric(str2double(strRunningTime)) && ~isnan(str2double(strRunningTime))
		dblMaxRunningTime = uint64(str2double(strRunningTime));
		
		%% start msg
		printf('Running time requested is number of repetitions: %d; [%s]\n',dblMaxRunningTime,getTime);
	else
		if ~ismember('-',strRunningTime),strRunningTime=['0-' strRunningTime];end
		dblDays = str2double(getFlankedBy(strRunningTime,'','-'));
		dblHours = str2double(getFlankedBy(strRunningTime,'-',':'));
		dblMins = str2double(getFlankedBy(strRunningTime,':',':'));
		dblSecs = str2double(getFlankedBy(getFlankedBy(strRunningTime,':',''),':',''));
		vecV = [24*60*60*dblDays 60*60*dblHours 60*dblMins dblSecs];
		vecV(isnan(vecV)) = 0;
		dblRunningTimeInput = sum(vecV);
		dblMaxRunningTime = dblRunningTimeInput - 60*5; %stop 5 minutes before requested end
		
		%% start msg
		printf('Running time requested: %.3fs; so we will be running for %.3fs [%s]\n',dblRunningTimeInput,dblMaxRunningTime,getTime);
	end
	
	%% set parameters
	sParams = struct;
	sParams.strHome = strHome;
	sParams.strOutputDir = strOutputDir;
	sParams.strStimDir = strStimDir;
	sParams.dblDeltaT = 0.5/1000; %seconds
	sParams.dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs
	sParams.intWorker = intWorker; %worker number
	sParams.dblAttention = dblAttention;
	sParams.intAttArea = intAttArea;
	
	%% load stimulus list
	[sStimParams,sStimInputs] = loadStimulation(strStimDir,strStimFile);
	
	%assign
	sParams.sStimParams = sStimParams;
	sParams.sStimInputs = sStimInputs;
	
	%msg
	printf('Loaded stimulation file <%s> [%s]\n',strStimFile,getTime);
	
	%% load connectivity
	[sConnParams,sConnectivity] = loadConnectivity_xArea(strConnDir,strConnFile);
	
	%assign
	sParams.sConnParams = sConnParams;
	sParams.sConnectivity = sConnectivity;
	sParams.strConnFile = strConnFile;
	
	% msg
	printf('Preparations complete; will be using connectivity matrix <%s> [%s]\n',strConnFile,getTime);
	
	%% seed random number generator
	rng('shuffle');
	intRandVals = randi(10^6);
	rng(intWorker);
	vecDummyVals = rand(1,intRandVals);clear vecDummyVals; %#ok<NASGU>
	sOut=rng;
	intSeed = sOut.Seed;
	intState = sOut.State(1);
	sParams.intSeed = intSeed;
	sParams.intState = intState;
	
	%% run simulation
	printf('Starting distributed run on compiled worker %d with random seed/state <%d/%d> [%s]\n',intWorker,intSeed,intState,getTime);
	
	%run
	if ~isa(dblMaxRunningTime,'uint64'),dblMaxRunningTime=dblMaxRunningTime-toc(hTic);end
	sData = getSimDistribLIF_xArea(sParams,dblMaxRunningTime,intWorker);
	
	%trial + time vars
	sSimRun = struct;
	sSimRun.vecOverallT = sData.vecOverallT;
	sSimRun.vecTrialStimType = sData.vecTrialStimType;
	sSimRun.vecTrialStimRep = sData.vecTrialStimRep;
	sSimRun.vecTrialStartSecs = sData.vecTrialStartSecs;
	sSimRun.vecStimStartSecs = sData.vecStimStartSecs;
	sSimRun.vecStimStopSecs = sData.vecStimStopSecs;
	sSimRun.vecTrialEndSecs = sData.vecTrialEndSecs;
	
	%spiking
	sSimRun.cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	sSimRun.cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	sSimRun.cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	
	%save data
	if ispc
		strDataFile = [strOutputDir 'xMaster_'  strTag sStimParams.strStimType sprintf('_%03d_%d_%s.mat',intSeed,intState,getDate)];
		printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
		save(strDataFile,'sSimRun','sParams','-v7.3');
	else
		strDataFile = [strOutputDir 'Simulation_' strTag sStimParams.strStimType sprintf('_%03d_%d_%s.mat',intSeed,intState,getDate)];
		printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
		sParams = rmfield(sParams,{'sConnParams','sConnectivity','sStimParams','sStimInputs'});
		save(strDataFile,'sSimRun','sParams','-v7');
	end
	printf(' .. Data saved! [%s]\n',getTime);
	
	%catch ME
	%	printf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);
	%	fprintf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);
	%end
