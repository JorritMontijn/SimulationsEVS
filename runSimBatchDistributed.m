%just to be sure
function runSimBatchDistributed(strRunningTime,strConnFile)
hTic = tic;

%% set log path
global strLogDir;
global boolClust;

%% running on cluster?
boolClust=false;
if boolClust
	
	%add directories to paths
	strHome = [filesep 'home' filesep 'montijn' filesep];
	if isdir([strHome 'Code']) && ~isdeployed,addpath([strHome 'Code']);end
	strLogDir = [strHome 'logs' filesep];
	strOutputDir = [strHome 'SimResults' filesep];
	strStimDriveDir = [strHome 'StimDrives' filesep];
	strConnDir = [strHome 'Connectivity' filesep];
	
	%add recursive subdirectories
	if ~isdeployed
		cellPaths = getSubDirs(strHome,inf,{'old','backup',' ','.','..','.mcrCache9.1'});
		fprintf('path 1: %s; nr of paths=%d\n',cellPaths{1},numel(cellPaths));
		for intPath=1:numel(cellPaths)
			addpath(cellPaths{intPath});
		end
	end
else
	%add directories to paths
	strFile = mfilename;
	strHome = mfilename('fullpath');
	strHome = strHome(1:(end-length(mfilename)));
	strLogDir = [strHome 'logs' filesep];
	strOutputDir = 'D:\Data\Processed\V1_LIFmodel\';
	strStimDriveDir = [strHome 'StimDrives' filesep];
	strConnDir = [strHome 'Connectivity' filesep];
end
boolClust = false;

%% transform running time to double
if ~ismember('-',strRunningTime),strRunningTime=['0-' strRunningTime];end
dblDays = str2double(getFlankedBy(strRunningTime,'','-'));
dblHours = str2double(getFlankedBy(strRunningTime,'-',':'));
dblMins = str2double(getFlankedBy(strRunningTime,':',':'));
dblSecs = str2double(getFlankedBy(getFlankedBy(strRunningTime,':',''),':',''));
vecV = [24*60*60*dblDays 60*60*dblHours 60*dblMins dblSecs];
vecV(isnan(vecV)) = 0;
dblRunningTimeInput = sum(vecV);
dblRunningTime = dblRunningTimeInput - 60*5; %stop 5 minutes before requested end

%% start msg
printf('Running time requested: %.3fs; so we will be running for %.3fs on each worker [%s]\n',dblRunningTimeInput,dblRunningTime,getTime);

%% start parallel pool
intCPUs = feature('numThreads');
printf('Number of available CPUs: %d\n',intCPUs);
if isempty(gcp('nocreate'))
    hPar = parpool('CPU multicore',20,'SpmdEnabled',false);
end
%% set parameters
sParams = struct;
sParams.strHome = strHome;
sParams.strOutputDir = strOutputDir;
sParams.strStimDriveDir = strStimDriveDir;
sParams.dblDeltaT = 0.5/1000; %seconds
sParams.dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
sStimParams = struct;
sStimParams.dblDeltaT = sParams.dblDeltaT;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0.1;
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 0.5;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [42.5 47.5];%[0:(180/12):179];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
sStimParams.vecSpatialFrequencies = 0.5;%2.^[-4:1];
sStimParams.vecTemporalFrequencies = 0;
sStimParams.intReps = 1; %number of repetitions
sStimParams.dblStimSizeRetDeg = 16;
sStimParams.vecScrPixWidthHeight = [128 128];
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
sStimParams.strStimDriveDir = strStimDriveDir;
sParams.sStimParams = sStimParams;

%% build stimuli
sParams.sStimInputs = loadSimStim(sStimParams);

%% load connectivity
if ~exist('strConnFile','var')
	strConnFile = 'sConn_Col48N2160S637056_2017-03-31.mat';
end
[sConnParams,sConnectivity] = loadConnectivity_xArea(strConnDir,strConnFile);

%assign
sParams.sConnParams = sConnParams;
sParams.sConnectivity = sConnectivity;
sParams.strConnFile = strConnFile;

%% msg
printf('Preparations complete; will be using connectivity matrix <%s> [%s]\n',strConnFile,getTime);

%% run simulation
sSimRun = struct;
parfor intWorker=1:intCPUs
	printf('Starting distributed run on worker %d/%d [%s]\n',intWorker,intCPUs,getTime);
	
	strLicenseErrorID = 'MATLAB:license:checkouterror';
	%run sim
	%boolRunning = true;
	%while boolRunning
	%	try
			sData = getSimDistribLIF_xArea(sParams,dblRunningTime-toc(hTic),intWorker);
	%		boolRunning = false;
	%	catch ME
	%		warning(ME);
	%		boolRunning = true;
	%		pause(60);
	%	end
	%end
	
	%trial + time vars
	sSimRun(intWorker).vecOverallT = sData.vecOverallT;
	sSimRun(intWorker).vecTrialStimType = sData.vecTrialStimType;
	sSimRun(intWorker).vecTrialStimRep = sData.vecTrialStimRep;
	sSimRun(intWorker).vecTrialStartSecs = sData.vecTrialStartSecs;
	sSimRun(intWorker).vecStimStartSecs = sData.vecStimStartSecs;
	sSimRun(intWorker).vecStimStopSecs = sData.vecStimStopSecs;
	sSimRun(intWorker).vecTrialEndSecs = sData.vecTrialEndSecs;
	
	%spiking
	sSimRun(intWorker).cellSpikeTimesLGN_ON = sData.cellSpikeTimesLGN_ON;
	sSimRun(intWorker).cellSpikeTimesLGN_OFF = sData.cellSpikeTimesLGN_OFF;
	sSimRun(intWorker).cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	
end

%save data
strRand = num2str(randi([0 9],10,1))';
strDataFile = [strOutputDir 'Simulation_'  getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('_%f_%s_%s.mat',now,getDate,strRand)];
printf(' .. Saving data to <%s>... [%s]\n',strDataFile,getTime);
save(strDataFile,'sSimRun','sParams');
printf(' .. Data saved! [%s]\n',getTime);

%catch ME
%	printf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);
%	fprintf('Error; ID: %s; msg: %s [%s]\n',ME.identifier,ME.message,getTime);
%end
