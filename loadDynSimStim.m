function [sStimParams,sStimInputs] = loadDynSimStim(sStimParams)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	%%
	%input
	if ~exist('sStimParams','var');sStimParams=struct;end
	if ~isfield(sStimParams,'vecScrPixWidthHeight'),sStimParams.vecScrPixWidthHeight = [128 128];end
	if ~isfield(sStimParams,'vecScrDegWidthHeight'),sStimParams.vecScrDegWidthHeight = [25.6 25.6];end
	if ~isfield(sStimParams,'dblStimSizeRetDeg'),sStimParams.dblStimSizeRetDeg = 16;end
	if ~isfield(sStimParams,'dblStartFirstTrialSecs'),sStimParams.dblStartFirstTrialSecs = 1;end
	if ~isfield(sStimParams,'dblPreStimBlankDur'),sStimParams.dblPreStimBlankDur = 0.25;end%0.25
	if ~isfield(sStimParams,'dblStimDur'),sStimParams.dblStimDur = 0.5;end%0.5
	if ~isfield(sStimParams,'dblPostStimBlankDur'),sStimParams.dblPostStimBlankDur = 0.25;end%0.25
	if ~isfield(sStimParams,'vecOrientations'),sStimParams.vecOrientations = 44.5 + [0 1];end
	if ~isfield(sStimParams,'vecSpatialFrequencies'),sStimParams.vecSpatialFrequencies = 0.2;end
	if ~isfield(sStimParams,'vecContrasts'),sStimParams.vecContrasts = 100;end
	if ~isfield(sStimParams,'vecLuminance'),sStimParams.vecLuminance = 100;end
	if ~isfield(sStimParams,'vecPhases'),sStimParams.vecPhases = 0;end
    if ~isfield(sStimParams,'vecGainMean'),sStimParams.vecGainMean = 1;end
    if ~isfield(sStimParams,'intReps'),sStimParams.intReps = 1;end %number of repetitions
	if ~isfield(sStimParams,'dblDeltaT'),sStimParams.dblDeltaT = 1;end %number of repetitions
	if ~isfield(sStimParams,'strStimType'),sStimParams.strStimType = 'SquareGrating';end
	if ~isfield(sStimParams,'boolUseAllCombs'),sStimParams.boolUseAllCombs = true;end
	strStimType = sStimParams.strStimType;
	cellParamTypes = {'SquareGrating','SineGrating','Line','NatMov'};
	intStimType = find(ismember(cellParamTypes,strStimType));
	if isempty(intStimType),error([mfilename ':StimTypeError'],sprintf('Stimulus type "%s" is not recognized [%s]',strStimType,getTime));end
	
	%check stim type
	boolUseAllCombs = sStimParams.boolUseAllCombs;
	dblDeltaT = sStimParams.dblDeltaT;
	dblStartFirstTrialSecs = sStimParams.dblStartFirstTrialSecs;
	dblPreStimBlankDur = sStimParams.dblPreStimBlankDur;%0.25
	dblStimDur = sStimParams.dblStimDur;%0.5
	dblPostStimBlankDur = sStimParams.dblPostStimBlankDur;%0.25
	vecOrientations = sStimParams.vecOrientations;
	vecContrasts = sStimParams.vecContrasts;
	vecLuminances = sStimParams.vecLuminances;
	vecPhases = sStimParams.vecPhases;
	vecGainMean = sStimParams.vecGainMean;
    vecTemporalFrequencies = sStimParams.vecTemporalFrequencies;
	if intStimType == 3
		sStimParams.vecSpatialFrequencies = 1/sStimParams.vecSpatialFrequencies;
	end
	vecSpatialFrequencies = sStimParams.vecSpatialFrequencies;
	intReps = sStimParams.intReps; %number of repetitions
	dblBlankDur = dblPreStimBlankDur + dblPostStimBlankDur;
	
	%build
	dblTrialDur = (dblBlankDur + dblStimDur);
	
	
	%% check noises
	if ~isfield(sStimParams,'vecOrientationNoise'),...
			sStimParams.vecOrientationNoise = zeros(size(sStimParams.vecOrientations));end
	if ~isfield(sStimParams,'vecSpatialFrequencyNoise'),...
			sStimParams.vecSpatialFrequencyNoise = zeros(size(sStimParams.vecSpatialFrequencies));end
	if ~isfield(sStimParams,'vecTemporalFrequencyNoise'),...
			sStimParams.vecTemporalFrequencyNoise = zeros(size(sStimParams.vecTemporalFrequencies));end
	if ~isfield(sStimParams,'vecContrastNoise'),...
			sStimParams.vecContrastNoise = zeros(size(sStimParams.vecContrasts));end
	if ~isfield(sStimParams,'vecLuminanceNoise'),...
			sStimParams.vecLuminanceNoise = zeros(size(sStimParams.vecLuminance));end
	if ~isfield(sStimParams,'vecPhaseNoise'),...
			sStimParams.vecPhaseNoise = zeros(size(sStimParams.vecPhases));end
	if ~isfield(sStimParams,'vecGainNoise'),...
			sStimParams.vecGainNoise = zeros(size(sStimParams.vecGainMean));end
	
	%assign
	vecOrientationNoise = sStimParams.vecOrientationNoise;
	vecSpatialFrequencyNoise = sStimParams.vecSpatialFrequencyNoise;
	vecTemporalFrequencyNoise = sStimParams.vecTemporalFrequencyNoise;
	vecContrastNoise = sStimParams.vecContrastNoise;
	vecLuminanceNoise = sStimParams.vecLuminanceNoise;
	vecPhaseNoise = sStimParams.vecPhaseNoise;
	vecGainNoise = sStimParams.vecGainNoise;
     
	%% BUILD STIM COMBINATIONS
	if boolUseAllCombs
		intOris = numel(vecOrientations);
		vecOriIdx = 1:intOris;
		
		intSFs = numel(vecSpatialFrequencies);
		vecSFIdx = 1:intSFs;
		
		intTFs = numel(vecTemporalFrequencies);
		vecTFIdx = 1:intTFs;
		
		intCs = numel(vecContrasts);
		vecCIdx = 1:intCs;
		
		intLums = numel(vecLuminances);
		vecLumIdx = 1:intLums;
		
		intPhs = numel(vecPhases);
		vecPhIdx = 1:intPhs;
		
        intGs = numel(vecGainMean);
        vecGIdx = 1:intGs;
		
        
		%make combos
		cellParamTypes = {vecOriIdx,vecSFIdx,vecTFIdx,vecCIdx,vecLumIdx,vecPhIdx,vecGIdx};
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
		
		vecStimTypeOriIdx = matStimTypeCombos(1,vecTrialStimType);
		vecStimTypeSFIdx = matStimTypeCombos(2,vecTrialStimType);
		vecStimTypeTFIdx = matStimTypeCombos(3,vecTrialStimType);
		vecStimTypeContrastIdx = matStimTypeCombos(4,vecTrialStimType);
		vecStimTypeLuminanceIdx = matStimTypeCombos(5,vecTrialStimType);
		vecStimTypePhaseIdx = matStimTypeCombos(6,vecTrialStimType);
		vecStimTypeGainIdx = matStimTypeCombos(7,vecTrialStimType);
		
        
        vecStimTypeOris = vecOrientations(vecStimTypeOriIdx);
		vecStimTypeOriNoise = vecOrientationNoise(vecStimTypeOriIdx);
		
		vecStimTypeSFs = vecSpatialFrequencies(vecStimTypeSFIdx);
		vecStimTypeSFNoise = vecSpatialFrequencyNoise(vecStimTypeSFIdx);
		
		vecStimTypeTFs = vecTemporalFrequencies(vecStimTypeTFIdx);
		vecStimTypeTFNoise = vecTemporalFrequencyNoise(vecStimTypeTFIdx);
		
		vecStimTypeContrasts = vecContrasts(vecStimTypeContrastIdx);
		vecStimTypeContrastNoise = vecContrastNoise(vecStimTypeContrastIdx);
		
		vecStimTypeLuminances = vecLuminances(vecStimTypeLuminanceIdx);
		vecStimTypeLuminanceNoise = vecLuminanceNoise(vecStimTypeLuminanceIdx);
		
		vecStimTypePhase = vecPhases(vecStimTypePhaseIdx);
		vecStimTypePhaseNoise = vecPhaseNoise(vecStimTypePhaseIdx);
		
        vecStimTypeGain = vecGainMean(vecStimTypeGainIdx);
		vecStimTypeGainNoise = vecGainNoise(vecStimTypeGainIdx);
		
	else
		%stim idx
		vecStimTypeOriIdx = label2idx(vecOrientations);
		vecStimTypeSFIdx = label2idx(vecSpatialFrequencies);
		vecStimTypeTFIdx = label2idx(vecTemporalFrequencies);
		vecStimTypeContrastIdx =label2idx(vecContrasts);
		vecStimTypeLuminanceIdx = label2idx(vecLuminances);
		vecStimTypePhaseIdx = label2idx(vecPhases);
		vecStimTypeGainIdx = label2idx(vecGainMean);
		
        
		%combos
		matStimTypeCombos = cat(1,vecStimTypeOriIdx,vecStimTypeSFIdx,vecStimTypeTFIdx,vecStimTypeContrastIdx,vecStimTypeLuminanceIdx,vecStimTypePhaseIdx,vecStimTypeGainIdx);
		
		%stim vals
		vecStimTypeOris = vecOrientations;
		vecStimTypeOriNoise = vecOrientationNoise;
		
		vecStimTypeSFs = vecSpatialFrequencies;
		vecStimTypeSFNoise = vecSpatialFrequencyNoise;
		
		vecStimTypeTFs = vecTemporalFrequencies;
		vecStimTypeTFNoise = vecTemporalFrequencyNoise;
		
		vecStimTypeContrasts = vecContrasts;
		vecStimTypeContrastNoise = vecContrastNoise;
		
		vecStimTypeLuminances = vecLuminances;
		vecStimTypeLuminanceNoise = vecLuminanceNoise;
		
		vecStimTypePhase = vecPhases;
		vecStimTypePhaseNoise = vecPhaseNoise;
		
        vecStimTypeGain = vecGainMean;
		vecStimTypeGainNoise = vecGainNoise;
		
		%make combos
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
		
		
	end
	
	%timing
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
	sStimParams.dblLuminance = vecLuminances(matStimTypeCombos(5,1));
	sStimParams.dblPhase = vecPhases(matStimTypeCombos(6,1));
	sStimParams.dblGain = vecGainMean(matStimTypeCombos(7,1));
    
    %sStimParams.dblStimSizeRetDeg = 16;
	%sStimParams.vecScrPixWidthHeight = [128 128];
	%sStimParams.vecScrDegWidthHeight = [25.6 25.6];
	sStimParams.dblDeltaT = dblDeltaT;
	sStimParams.dblStimDur = dblStimDur;
	sStimParams.dblBlankDur = dblBlankDur;
	
	%get stim drive
	sStimDrive = getDynBottomUpInputs(sStimParams);
	
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
	
	%trial timing
	sStimInputs.vecTrialStartSecs = vecTrialStartSecs;
	sStimInputs.vecTrialEndSecs = vecTrialEndSecs;
	sStimInputs.vecStimStartSecs = vecStimStartSecs;
	sStimInputs.vecStimStopSecs = vecStimStopSecs;
	
	%trial stim type
	sStimInputs.vecTrialStimType = vecTrialStimType;
	sStimInputs.vecTrialStimRep = vecTrialStimRep;
	
	%stim type props
	sStimInputs.matStimTypeCombos = matStimTypeCombos;
	sStimInputs.vecStimTypeOris = vecStimTypeOris;
	sStimInputs.vecStimTypeOriNoise = vecStimTypeOriNoise;
	sStimInputs.vecStimTypeOriIdx = vecStimTypeOriIdx;
	sStimInputs.vecStimTypeSFs = vecStimTypeSFs;
	sStimInputs.vecStimTypeSFNoise = vecStimTypeSFNoise;
	sStimInputs.vecStimTypeSFIdx = vecStimTypeSFIdx;
	sStimInputs.vecStimTypeTFs = vecStimTypeTFs;
	sStimInputs.vecStimTypeTFNoise = vecStimTypeTFNoise;
	sStimInputs.vecStimTypeTFIdx = vecStimTypeTFIdx;
	sStimInputs.vecStimTypeContrasts = vecStimTypeContrasts;
	sStimInputs.vecStimTypeContrastNoise = vecStimTypeContrastNoise;
	sStimInputs.vecStimTypeContrastIdx = vecStimTypeContrastIdx;
	sStimInputs.vecStimTypeLuminances = vecStimTypeLuminances;
	sStimInputs.vecStimTypeLuminanceNoise = vecStimTypeLuminanceNoise;
	sStimInputs.vecStimTypeLuminanceIdx = vecStimTypeLuminanceIdx;
	sStimInputs.vecStimTypePhase = vecStimTypePhase;
	sStimInputs.vecStimTypePhaseNoise = vecStimTypePhaseNoise;
	sStimInputs.vecStimTypePhaseIdx = vecStimTypePhaseIdx;
    sStimInputs.vecStimTypeGain = vecStimTypeGain;
	sStimInputs.vecStimTypeGainNoise = vecStimTypeGainNoise;
	sStimInputs.vecStimTypeGainIdx = vecStimTypeGainIdx;
    	
end

