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
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 500/1000;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [42.5 47.5];%[44 46 44 44 46]; %[44 46];%[0:(180/8):179];%[44 46 44 44 46];
sStimParams.vecOrientationNoise = 5*ones(size(sStimParams.vecOrientations));%[5 0 0 0];%[0 5 0 0 5];
sStimParams.vecSpatialFrequencies = 0.25;%[0.25 0.25 0.2 0.25 0.2];%[0.25 0.25 0.2 0.25 0.2];%2.^[-4:1];
sStimParams.vecSpatialFrequencyNoise = 0;%0.05*ones(size(sStimParams.vecSpatialFrequencies));%[0.05 0 0 0];%[0 0 0.1 0 0.1];
sStimParams.vecTemporalFrequencies = 2;%[2 2 2 2.5 2.5];%[2 2 2 2.5 2.5];
sStimParams.vecTemporalFrequencyNoise = 0;%0.5*ones(size(sStimParams.vecSpatialFrequencies));%[0.5 0 0 0];%[0 0 0 1 1];
sStimParams.vecContrasts = 100;%ones(size(sStimParams.vecOrientations))*100;
sStimParams.vecContrastNoise = 0;%zeros(size(sStimParams.vecOrientations))*100;
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
sStimParams.strStimTag= sprintf('ExpRet%dNoise%sOri%dDrift%d',...
	sStimParams.vecScrPixWidthHeight(1)/2,...
	num2str(sStimParams.vecOrientationNoise(1)),...'All',...
	round(range(sStimParams.vecOrientations(1:2))),numel(sStimParams.vecOrientations));

%% save
strStimFile = sprintf('sStim2_%s%s_x%dR%d_%s.mat',sStimParams.strStimType,sStimParams.strStimTag,numel(sStimInputs.vecTrialStimType),max(sStimInputs.vecTrialStimRep),getDate);
strStimDir = 'D:\Simulations\Stimulation\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strStimFile,strStimDir,getTime);
save([strStimDir strStimFile],'sStimParams','sStimInputs');
