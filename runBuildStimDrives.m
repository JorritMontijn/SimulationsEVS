%% msg
clearvars;
fprintf('Starting offline construction of stimulus profile... [%s]\n\n',getTime);

%% prepare stimulus list
sStimParams = struct;
sStimParams.dblDeltaT = 0.5/1000;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0;
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 0.5;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [0:(180/18):179];%[44.5 45.5];%[0:(180/18):179];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
sStimParams.vecSpatialFrequencies = 0.25;%2.^[-4:1];
sStimParams.vecTemporalFrequencies = 2;
sStimParams.vecContrasts = 100;
sStimParams.vecLuminance = 100;
sStimParams.intReps = 1; %number of repetitions
sStimParams.dblStimSizeRetDeg = 5;
sStimParams.vecScrPixWidthHeight = [128 128];
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
sStimParams.strStimDriveDir = 'D:\Simulations\StimDrives\';

%% build stimuli
sStimInputs = loadSimStim(sStimParams);
sStimInputs.vecStimTypeAttention = zeros(size(sStimInputs.vecTrialStimType));
sStimParams.strStimType = [sStimParams.strStimType 'OriDrift'];

%% edit
intDoEdit = 0
if intDoEdit == 1
	%build stimulus drives
	matBlank = sStimInputs.matBlankLGN_ON;
	dblContrast = sStimInputs.cellContrast{1};
	sStimInputs.matBlankLGN_ON = 0*sStimInputs.matBlankLGN_ON;
	sStimInputs.matBlankLGN_OFF = 0*sStimInputs.matBlankLGN_OFF;
	%vecTrialInputStrength = 2.^(-3:0.01:3);
	vecTrialInputStrength = 2.^(-1:0.0025:2);
	vecTrialSpikeProbLGN = vecTrialInputStrength*mean(matBlank(:));
	intTrials = numel(vecTrialInputStrength);
	
	%pre-allocate
	sStimInputs.cellR_ON = [];
	sStimInputs.cellR_OFF = [];
	sStimInputs.vecTrialOris = zeros(1,intTrials);
	sStimInputs.vecTrialOriIdx = ones(1,intTrials);
	sStimInputs.vecTrialSFs = zeros(1,intTrials);
	sStimInputs.vecTrialSFIdx = ones(1,intTrials);
	sStimInputs.vecTrialTFs = zeros(1,intTrials);
	sStimInputs.vecTrialTFIdx = ones(1,intTrials);
	sStimInputs.vecTrialContrasts = zeros(1,intTrials);
	sStimInputs.vecTrialContrastIdx = ones(1,intTrials);
	
	sStimInputs.vecTrialStimType = 1:intTrials;
	sStimInputs.vecTrialStimRep = ones(1,intTrials);
	sStimInputs.vecTrialStartSecs = sStimParams.dblStartFirstTrialSecs:sStimInputs.dblTrialDur:(sStimInputs.dblTrialDur*intTrials+sStimParams.dblStartFirstTrialSecs-sStimInputs.dblTrialDur/2);
	sStimInputs.vecStimStartSecs = sStimInputs.vecTrialStartSecs+sStimParams.dblPreStimBlankDur;
	sStimInputs.vecStimStopSecs = sStimInputs.vecStimStartSecs+sStimParams.dblStimDur;
	sStimInputs.vecTrialEndSecs = sStimInputs.vecStimStopSecs+sStimParams.dblPostStimBlankDur;
	
	sStimInputs.vecStimTypeOris = 0;
	sStimInputs.vecStimTypeSFs = 0;
	sStimInputs.vecStimTypeTFs = 0;
	sStimInputs.vecStimTypeContrasts = 0;
	
	%assign LGN matrices
	for intTrial=1:intTrials
		dblStrength = vecTrialInputStrength(intTrial);
		sStimInputs.cellLGN_ON{intTrial} = dblStrength*matBlank;
		sStimInputs.cellLGN_OFF{intTrial} = dblStrength*matBlank;
		sStimInputs.cellContrast{intTrial} = dblContrast*dblStrength;
	end
	sStimInputs.vecTrialInputStrength = vecTrialInputStrength;
	sStimInputs.vecTrialSpikeProbLGN = vecTrialSpikeProbLGN;
	sStimParams.strStimType = 'ContrastDependency';
end

%% save
strStimFile = sprintf('sStim_%s_x%dR%d_%s.mat',sStimParams.strStimType,numel(sStimInputs.vecTrialStimType),max(sStimInputs.vecTrialStimRep),getDate);
strStimDir = 'D:\Simulations\Stimulation\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strStimFile,strStimDir,getTime);
save([strStimDir strStimFile],'sStimParams','sStimInputs');
