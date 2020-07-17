%runBuildConnectivity	Creates the input that will feed into the network
%						through a simplified retina/LGN filter layer. 
%
%This function can create visual stimuli in retinal degrees, just as you
%would for an actual visual neuroscience experiment. You can set the
%stimulation parameters by editing the fields of the sStimParams structure
%in this file. See the function loadDynSimStim() for more information.
%
%NOTE: I'm planning to update the stimulation scripts to use the same
%functions as I use for creating an actual visual stimulus
%
%Version History:
%2019-02-12 Updated name and help description [by Jorrit Montijn]

%% msg
clearvars;
fprintf('Starting offline construction of stimulus profile... [%s]\n\n',getTime);

%% prepare stimulus list
boolSmallRet = true;
if boolSmallRet,dblDivFactor = 2;else dblDivFactor = 1;end
sStimParams = struct;
sStimParams.dblDeltaT = 0.5/1000;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0;
sStimParams.dblPreStimBlankDur = 0.5;%0.25
sStimParams.dblStimDur = 500/1000;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = 40;%[0:20:359];%[44 46 44 44 46]; %[44 46];%[0:(180/8):179];%[44 46 44 44 46];
sStimParams.vecOrientationNoise = 0*ones(size(sStimParams.vecOrientations));%[5 0 0 0];%[0 5 0 0 5];
sStimParams.vecSpatialFrequencies = 0.25;%[0.25 0.25 0.2 0.25 0.2];%[0.25 0.25 0.2 0.25 0.2];%2.^[-4:1];
sStimParams.vecSpatialFrequencyNoise = 0;%0.05*ones(size(sStimParams.vecSpatialFrequencies));%[0.05 0 0 0];%[0 0 0.1 0 0.1];
sStimParams.vecTemporalFrequencies = 2;%[2 2 2 2.5 2.5];%[2 2 2 2.5 2.5];
sStimParams.vecTemporalFrequencyNoise = 0;%0.5*ones(size(sStimParams.vecSpatialFrequencies));%[0.5 0 0 0];%[0 0 0 1 1];
sStimParams.vecContrasts = [0 2.^(-1:0.5:6) 100];%ones(size(sStimParams.vecOrientations))*100;
sStimParams.vecContrastNoise = 0*ones(size(sStimParams.vecContrasts));
sStimParams.vecLuminances = 100;%ones(size(sStimParams.vecOrientations))*100;
sStimParams.vecLuminanceNoise = 0;%zeros(size(sStimParams.vecOrientations))*100;
sStimParams.vecPhases = 0;%zeros(size(sStimParams.vecOrientations))*100;
sStimParams.vecPhaseNoise = 10^6;%zeros(size(sStimParams.vecOrientations))*100;
sStimParams.vecGainMean = 1;%zeros(size(sStimParams.vecOrientations))*100;
sStimParams.vecGainNoise = 0;%zeros(size(sStimParams.vecOrientations))*100;
sStimParams.intReps = 1; %number of repetitions
sStimParams.dblStimSizeRetDeg = 5;
sStimParams.vecScrPixWidthHeight = [128 128]/dblDivFactor;
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
sStimParams.strStimDriveDir = 'D:\Simulations\StimDrives\';
sStimParams.boolUseAllCombs = true;%false

%% build stimuli
[sStimParams,sStimInputs] = loadDynSimStim(sStimParams);
sStimInputs.vecStimTypeAttention = zeros(size(sStimInputs.vecTrialStimType));
sStimParams.strStimTag= sprintf('ExpRet%dCont%sDrift%d',...
	sStimParams.vecScrPixWidthHeight(1)/2,...
	num2str(numel(sStimParams.vecContrasts)),...'All',...
	numel(sStimParams.vecOrientations));

%% build stimuli v2
%{
%% specific to simulations
%add model parameters
sStimModel.vecScrPixWidthHeight = [128 128]/dblDivFactor; %OLD
sStimParams.dblDeltaT = 0.5/1000;

%% stimulus params
%visual space parameters
sStimParams = struct;
sStimParams.strStimType = 'SquareGrating';
sStimParams.dblSubjectPosX_cm = 0; % cm; relative to center of screen
sStimParams.dblSubjectPosY_cm = 0; % cm; relative to center of screen
sStimParams.dblScreenDistance_cm = 100; % cm; measured
sStimParams.vecUseMask = 1; %[1] if mask to emulate retinal-space, [0] use screen-space

%receptive field size&location parameters
sStimParams.vecStimPosX_deg = 0; % deg; relative to subject
sStimParams.vecStimPosY_deg = 0; % deg; relative to subject
sStimParams.vecStimulusSize_deg = 5;%circular window in degrees [35]
sStimParams.vecSoftEdge_deg = 0.2; %width of cosine ramp  in degrees, [0] is hard edge

%screen variables
sStimParams.intUseScreen = 1; %which screen to use
sStimParams.intCornerTrigger = 0; % integer switch; 0=none,1=upper left, 2=upper right, 3=lower left, 4=lower right
sStimParams.dblCornerSize = 0; % fraction of screen width

sStimParams.dblScreenWidth_cm = 45.5; % cm; measured [51]
sStimParams.dblScreenHeight_cm = 45.5; % cm; measured [29]
sStimParams.dblScreenWidth_deg = 2 * atand(sStimParams.dblScreenWidth_cm / (2 * sStimParams.dblScreenDistance_cm));
sStimParams.dblScreenHeight_deg = 2 * atand(sStimParams.dblScreenHeight_cm / (2 * sStimParams.dblScreenDistance_cm));

%stimulus control variables
sStimParams.intUseParPool = 0; %number of workers in parallel pool; [2]
sStimParams.intUseGPU = 1;
sStimParams.intAntiAlias = 0;
sStimParams.str90Deg = '0 degrees is leftward motion; 90 degrees is upward motion';
sStimParams.vecBackgrounds = 0.5; %background intensity (dbl, [0 1])

sStimParams.intBackground = round(mean(sStimParams.vecBackgrounds)*255);
sStimParams.vecContrasts = 100; %contrast; [0-100]
sStimParams.vecOrientations = 0:20:359;%[357 3 24 45 66 87 93 114 135 156 177 183 204 225 246 267 273 294 315 336]; %orientation (0 is drifting rightward)
sStimParams.vecSpatialFrequencies = 0.08; %Spat Frequency in cyc/deg 0.08
sStimParams.vecTemporalFrequencies = 1; %Temporal frequency in cycles per second (0 = static gratings only)

%build single-repetition list
[sStimParams,sStimObject,sStimTypeList] = getDriftingGratingCombos(sStimParams);
%}
%% save
strStimFile = sprintf('sStim2_%s%s_x%dR%d_%s.mat',sStimParams.strStimType,sStimParams.strStimTag,numel(sStimInputs.vecTrialStimType),max(sStimInputs.vecTrialStimRep),getDate);
strStimDir = 'C:\Code\SimulationsEVS\Stimulation\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strStimFile,strStimDir,getTime);
save([strStimDir strStimFile],'sStimParams','sStimInputs');
