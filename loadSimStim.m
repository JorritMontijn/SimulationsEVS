function sStimInputs = loadSimStim(sStimParams)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%input
	if ~exist('sStimParams','var');sStimParams=struct;end
	if ~isfield(sStimParams,'dblStartFirstTrialSecs'),sStimParams.dblStartFirstTrialSecs = 1;end
	if ~isfield(sStimParams,'dblPreStimBlankDur'),sStimParams.dblPreStimBlankDur = 0.25;end%0.25
	if ~isfield(sStimParams,'dblStimDur'),sStimParams.dblStimDur = 0.5;end%0.5
	if ~isfield(sStimParams,'dblPostStimBlankDur'),sStimParams.dblPostStimBlankDur = 0.25;end%0.25
	if ~isfield(sStimParams,'vecOrientations'),sStimParams.vecOrientations = 44.5 + [0 1];end
	if ~isfield(sStimParams,'vecSpatialFrequencies'),sStimParams.vecSpatialFrequencies = 0.2;end
	if ~isfield(sStimParams,'vecContrasts'),sStimParams.vecContrasts = 100;end
	if ~isfield(sStimParams,'vecLuminance'),sStimParams.vecLuminance = 100;end
	if ~isfield(sStimParams,'intReps'),sStimParams.intReps = 1;end %number of repetitions
	if ~isfield(sStimParams,'dblDeltaT'),sStimParams.dblDeltaT = 1;end %number of repetitions
	if ~isfield(sStimParams,'strStimType'),sStimParams.strStimType = 'SquareGrating';end
	strStimType = sStimParams.strStimType;
	cellParamTypes = {'SquareGrating','SineGrating','Line','NatMov'};
	intStimType = find(ismember(cellParamTypes,strStimType));
	if isempty(intStimType),error([mfilename ':StimTypeError'],sprintf('Stimulus type "%s" is not recognized [%s]',strStimType,getTime));end
	
	
	%check stim type
	dblDeltaT = sStimParams.dblDeltaT;
	dblStartFirstTrialSecs = sStimParams.dblStartFirstTrialSecs;
	dblPreStimBlankDur = sStimParams.dblPreStimBlankDur;%0.25
	dblStimDur = sStimParams.dblStimDur;%0.5
	dblPostStimBlankDur = sStimParams.dblPostStimBlankDur;%0.25
	vecOrientations = sStimParams.vecOrientations;
	vecContrasts = sStimParams.vecContrasts;
	vecLuminance = sStimParams.vecLuminance;
	vecTemporalFrequencies = sStimParams.vecTemporalFrequencies;
	if intStimType == 3
		sStimParams.vecSpatialFrequencies = 1/sStimParams.vecSpatialFrequencies;
	end
	vecSpatialFrequencies = sStimParams.vecSpatialFrequencies;
	intReps = sStimParams.intReps; %number of repetitions
	dblBlankDur = dblPreStimBlankDur + dblPostStimBlankDur;
	
	%build
	dblTrialDur = (dblBlankDur + dblStimDur);
	
	intOris = numel(vecOrientations);
	vecOriIdx = 1:intOris;
	
	intSFs = numel(vecSpatialFrequencies);
	vecSFIdx = 1:intSFs;
	
	intTFs = numel(vecTemporalFrequencies);
	vecTFIdx = 1:intTFs;
	
	intCs = numel(vecContrasts);
	vecCIdx = 1:intCs;
	
	intLums = numel(vecLuminance);
	vecLumIdx = 1:intLums;
	
	%% BUILD STIM COMBINATIONS
	cellParamTypes = {vecOriIdx,vecSFIdx,vecTFIdx,vecCIdx,vecLumIdx};
	matStimTypeCombos = buildStimCombos(cellParamTypes);
	intStimTypes = size(matStimTypeCombos,2);
	vecStimIdx = 1:intStimTypes;
	
	%build trial presentations
	intTrials = intStimTypes*intReps;
	
	if intReps==1 %do not shuffle
		vecTrialStimType = vecStimIdx;
		vecTrialStimRep = ones(1,intStimTypes);
	else
		vecTrialStimType = [];
		vecTrialStimRep = [];
		for intRep=1:intReps
			vecShuffle = randperm(intStimTypes);
			vecTrialStimType = [vecTrialStimType vecStimIdx(vecShuffle)];
			vecTrialStimRep = [vecTrialStimRep intRep*ones(1,intStimTypes)];
		end
	end
	
	vecTrialOriIdx = matStimTypeCombos(1,vecTrialStimType);
	vecTrialSFIdx = matStimTypeCombos(2,vecTrialStimType);
	vecTrialTFIdx = matStimTypeCombos(3,vecTrialStimType);
	vecTrialContrastIdx = matStimTypeCombos(4,vecTrialStimType);
	vecTrialLuminanceIdx = matStimTypeCombos(5,vecTrialStimType);
	
	vecTrialOris = vecOrientations(vecTrialOriIdx);
	vecTrialSFs = vecSpatialFrequencies(vecTrialSFIdx);
	vecTrialTFs = vecTemporalFrequencies(vecTrialTFIdx);
	vecTrialContrasts = vecContrasts(vecTrialContrastIdx);
	vecTrialLuminance = vecLuminance(vecTrialLuminanceIdx);
	
	vecTrialStartSecs = dblStartFirstTrialSecs:dblTrialDur:(dblTrialDur*intTrials+dblStartFirstTrialSecs-dblTrialDur/2);
	vecStimStartSecs = vecTrialStartSecs+dblPreStimBlankDur;
	vecStimStopSecs = vecStimStartSecs+dblStimDur;
	vecTrialEndSecs = vecStimStopSecs+dblPostStimBlankDur;
	fprintf(' .. Created stimulus presentation list: %.0f ms stimulus duration, %.0f ms blanks; %d orientations, %d repeats. Total of %d trials\n',dblStimDur*1000,(dblPreStimBlankDur+dblPostStimBlankDur)*1000,numel(vecOrientations),intReps,intTrials);
	
	%make cell arrays
	cellR_ON = cell(1,intStimTypes);
	cellR_OFF = cell(1,intStimTypes);
	cellLGN_ON = cell(1,intStimTypes);
	cellLGN_OFF = cell(1,intStimTypes);
	cellContrast = cell(1,intStimTypes);
	cellLuminance = cell(1,intStimTypes);
	
	%get first type
	% build/load stimuli and LGN responses
	sStimParams.intFrameRate = 100;
	sStimParams.dblAngleInDeg = vecOrientations(matStimTypeCombos(1,1));
	sStimParams.dblSF = vecSpatialFrequencies(matStimTypeCombos(2,1));
	sStimParams.dblTF = vecTemporalFrequencies(matStimTypeCombos(3,1));
	sStimParams.dblContrast = vecContrasts(matStimTypeCombos(4,1));
	sStimParams.dblLuminance = vecLuminance(matStimTypeCombos(5,1));
	%sStimParams.dblStimSizeRetDeg = 16;
	sStimParams.vecScrPixWidthHeight = [128 128];
	sStimParams.vecScrDegWidthHeight = [25.6 25.6];
	sStimParams.dblDeltaT = dblDeltaT;
	sStimParams.dblStimDur = dblStimDur;
	sStimParams.dblBlankDur = dblBlankDur;
	
	%get stim drive
	sStimDrive = getBottomUpInputs(sStimParams);

	%assign
	cellR_ON{1} = sStimDrive.matR_ON;
	cellR_OFF{1} = sStimDrive.matR_OFF;
	cellLGN_ON{1} = sStimDrive.matLGN_ON;
	cellLGN_OFF{1} = sStimDrive.matLGN_OFF;
	cellContrast{1} = sStimDrive.vecContrast;
	cellLuminance{1} = sStimDrive.vecLuminance;
	varDeltaSyn = sStimDrive.varDeltaSyn;
	sStimParams.varDeltaSyn = varDeltaSyn;
	matBlankLGN_ON = sStimDrive.matLGN_ON(:,:,end);
	matBlankLGN_OFF = sStimDrive.matLGN_OFF(:,:,end);
	dblVisSpacing = sStimDrive.dblVisSpacing;
	
	%run other stimuli
	for intStimType=2:intStimTypes
		% build/load stimuli and LGN responses
		sStimParams.dblAngleInDeg = vecOrientations(matStimTypeCombos(1,intStimType));
		sStimParams.dblSF = vecSpatialFrequencies(matStimTypeCombos(2,intStimType));
		sStimParams.dblTF = vecTemporalFrequencies(matStimTypeCombos(3,intStimType));
		sStimParams.dblContrast = vecContrasts(matStimTypeCombos(4,intStimType));
		sStimParams.dblLuminance = vecLuminance(matStimTypeCombos(5,intStimType));
		
		sStimDrive = getBottomUpInputs(sStimParams);
		
		cellR_ON{intStimType} = sStimDrive.matR_ON;
		cellR_OFF{intStimType} = sStimDrive.matR_OFF;
		cellLGN_ON{intStimType} = sStimDrive.matLGN_ON;
		cellLGN_OFF{intStimType} = sStimDrive.matLGN_OFF;
		cellContrast{intStimType} = sStimDrive.vecContrast;
		cellLuminance{intStimType} = sStimDrive.vecLuminance;
	end
	
	%output
	sStimInputs.dblTrialDur = dblTrialDur;
	sStimInputs.dblVisSpacing = dblVisSpacing;
	sStimInputs.varDeltaSyn = varDeltaSyn;
	sStimInputs.matBlankLGN_ON = matBlankLGN_ON;
	sStimInputs.matBlankLGN_OFF = matBlankLGN_OFF;
	sStimInputs.cellR_ON = cellR_ON;
	sStimInputs.cellR_OFF = cellR_OFF;
	sStimInputs.cellLGN_ON = cellLGN_ON;
	sStimInputs.cellLGN_OFF = cellLGN_OFF;
	sStimInputs.cellContrast = cellContrast;
	sStimInputs.cellLuminance = cellLuminance;
	sStimInputs.vecTrialOris = vecTrialOris;
	sStimInputs.vecTrialOriIdx = vecTrialOriIdx;
	sStimInputs.vecTrialSFs = vecTrialSFs;
	sStimInputs.vecTrialSFIdx = vecTrialSFIdx;
	sStimInputs.vecTrialTFs = vecTrialTFs;
	sStimInputs.vecTrialTFIdx = vecTrialTFIdx;
	sStimInputs.vecTrialContrasts = vecTrialContrasts;
	sStimInputs.vecTrialContrastIdx = vecTrialContrastIdx;
	sStimInputs.vecTrialLuminance = vecTrialLuminance;
	sStimInputs.vecTrialLuminanceIdx = vecTrialLuminanceIdx;
	sStimInputs.vecTrialStimType = vecTrialStimType;
	sStimInputs.vecTrialStimRep = vecTrialStimRep;
	sStimInputs.vecTrialStartSecs = vecTrialStartSecs;
	sStimInputs.vecTrialEndSecs = vecTrialEndSecs;
	sStimInputs.vecStimStartSecs = vecStimStartSecs;
	sStimInputs.vecStimStopSecs = vecStimStopSecs;
	sStimInputs.matStimTypeCombos = matStimTypeCombos;
	sStimInputs.vecStimTypeOris = vecOrientations(matStimTypeCombos(1,:));
	sStimInputs.vecStimTypeSFs = vecSpatialFrequencies(matStimTypeCombos(2,:));
	sStimInputs.vecStimTypeTFs = vecTemporalFrequencies(matStimTypeCombos(3,:));
	sStimInputs.vecStimTypeContrasts = vecContrasts(matStimTypeCombos(4,:));
	sStimInputs.vecStimTypeLuminance = vecLuminance(matStimTypeCombos(5,:));
end

