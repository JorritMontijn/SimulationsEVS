function [matR_ON,matR_OFF,dblVisSpacing,vecLuminance] = getDynRetinalDriveDirect(sParams,sP)
	%getRetinalDrive Produces retinal ON and OFF field responses to an
	%input image. Syntax:
	%   [matR_ON,matR_OFF] = getRetinalDrive(matImage,sParams)
	%
	%matR_ON and matR_OFF are 3D matrices [x,y,t] containing the retinal
	%drive given the input matImage also of size [x2,y2,t], given
	%parameters defined by their fieldnames in sParams. If a particular
	%field is absent from the input, the default parameters are chosen as
	%follows:
	%
	%dblContrast = 100;				%stimulus contrast [1-100%]
	%dblStepT = 0.5;				%Time-step in ms
	%dblVisSpacingImage = 0.2;		%Visual spacing between pixels in degrees
	%vecArraySize = [21 21];		%Size of OFF/ON field array (nr of cells)
	%dblVisSpacingCells = 0.2;		%Visual spacing between cells in degrees
	%dblDelayCS = 3;				%delay in ms between center and surround
	%dblSigmaCenter = 0.176;		%center size in degrees
	%dblSigmaSurround = 0.53;		%surround size in degrees
	%dblK_center = 17;				%weighting variable
	%dblK_surround = 16;			%weighting variable
	%dblTauCenter = 10;				%in ms, controls exp. drop-off
	%dblTauSurround = 20;			%same as above, for center-field
	%dblR_baseline = 15;			%baseline spiking rate in Hz
	
	%% run
	%get default values
	if ~exist('sP','var'),sP=struct;end
	if ~isfield(sP,'dblLuminance'),dblLuminance = 100;else dblLuminance=sP.dblLuminance;end
	if ~isfield(sP,'dblDelayCS'),dblDelayCS = 3;else dblDelayCS=sP.dblDelayCS;end
	if ~isfield(sP,'dblSigmaCenter'),dblSigmaCenter = 0.176;else dblSigmaCenter=sP.dblSigmaCenter;end
	if ~isfield(sP,'dblSigmaSurround'),dblSigmaSurround = 0.53;else dblSigmaSurround=sP.dblSigmaSurround;end
	if ~isfield(sP,'dblK_center'),dblK_center = 17;else dblK_center=sP.dblK_center;end
	if ~isfield(sP,'dblK_surround'),dblK_surround = 16;else dblK_surround=sP.dblK_surround;end
	if ~isfield(sP,'dblTauCenter'),dblTauCenter = 10;else dblTauCenter=sP.dblTauCenter;end
	if ~isfield(sP,'dblTauSurround'),dblTauSurround = 20;else dblTauSurround=sP.dblTauSurround;end
	if ~isfield(sP,'dblR_baseline'),dblR_baseline = 15;else dblR_baseline=sP.dblR_baseline;end
	
	%% get values
	%transform format
	dblTF = sParams.TF; %temporal frequency; 5
	dblSF = sParams.SF; %spatial frequency; 2
	dblStimSizeRetDeg = sParams.SSRD; %stim size; 5
	intArraySize = sParams.PixWH; %array size; 256
	dblArrayDeg = sParams.DegWH; %array size in visual degrees; 6.4
	dblAngleInDeg = sParams.Ori; %stim orientation; 42.5
	dblPhase = sParams.Phase; %starting phase; 0
    dblDeltaT = sParams.dT; %time step; 0.5/1000
	dblStimDur = sParams.SD; %stim dur; 0.2
	dblBlankDur = sParams.BD; %blank dur, half before, half after; 0.2
	dblVisSpacing = dblArrayDeg/intArraySize;
	vecArraySize = [intArraySize intArraySize];
	
	%% prep
	%get image size
	dblAngle = -deg2rad(dblAngleInDeg);
	intImX = intArraySize;
	intImY = intArraySize;
	intImT = round((dblStimDur + dblBlankDur)/dblDeltaT);
	intBaseT = round((dblBlankDur/2)/dblDeltaT);
	
	%get spatial receptive fields
	vecSpaceX = dblVisSpacing*((-(intImX - 1)/2):(intImX - 1)/2);
	vecSpaceY = dblVisSpacing*((-(intImY - 1)/2):(intImY - 1)/2);
	[matMeshX,matMeshY] = meshgrid(vecSpaceX,vecSpaceY);
	
	%define variable
	psi = dblPhase;
	k_x = 2*pi*dblSF*cos(dblAngle); % But assuming x is in degree and you want a spatial frequency of q cycles/degree, then |k| should be set to 2*pi*q. Note that |k| = sqrt(k_x^2 + k_y^2).
	k_y = 2*pi*dblSF*sin(dblAngle); % But assuming x is in degree and you want a spatial frequency of q cycles/degree, then |k| should be set to 2*pi*q. Note that |k| = sqrt(k_x^2 + k_y^2).
	k_abs = sqrt(k_x^2 + k_y^2);
	x = [matMeshX(:) matMeshY(:)];
	vecT = dblDeltaT:dblDeltaT:dblStimDur; %time after onset
	sigma_c = dblSigmaCenter; %spatial size of the filter (center)
	sigma_s = dblSigmaSurround; %spatial size of the filter (surround)
	omega = dblTF*2*pi;
	tau_c = dblTauCenter; %dblTauCenter or dblTauSurround
	tau_s = dblTauSurround; %dblTauCenter or dblTauSurround
	
	%% run
	matR_center = nan([intImX intImY numel(vecT)]);
	matR_surround = nan([intImX intImY numel(vecT)]);
	%calculate
	for intT=1:numel(vecT)
		t = vecT(intT);
		
		%%
		for intCS=1:2
			if intCS == 1 %center
				sigma = sigma_c; %spatial size of the filter (surround)
				tau = tau_c;
			else %surround
				sigma = sigma_s; %spatial size of the filter (surround)
				tau = tau_s;
			end
			nu = 1/tau;
			nu = nu/dblDeltaT;
			
			preTerm = exp(-((k_abs^2)*(sigma^2))/2)/(1+(omega^2)*(tau^2));
			k = [k_x k_y];
			vecKX = k*x';
			
			I = preTerm*(cos(vecKX - omega*t+psi) - omega*tau*sin(vecKX-omega*t+psi)...
				-exp(-nu*t)*(cos(vecKX+psi) - omega*tau*sin(vecKX+psi)));
			
			I = reshape(I,[intImX intImY]);
			if intCS == 1%center
				matR_center(:,:,intT) = (dblK_center*I);
			else%surround
				matR_surround(:,:,intT) = (dblK_surround*I);
			end
		end
	end
	%% window
	%calculate window variables
	dblWindowDiameter = dblStimSizeRetDeg / dblVisSpacing; %size of non-ramped stimulus
	dblRampWidth = 0.5/dblVisSpacing; %size of cosine ramp where the stimulus goes from maximal intensities to background grey
	
	%build the circular ramp using another function and also create an inverse mask
	matRampGrid = buildCircularCosineRampPix(intImX,dblWindowDiameter,dblRampWidth);
	%matRampGridInverse = abs(matRampGrid - 1);
	
	matR_center = bsxfun(@times,matR_center,matRampGrid); %multiply the ramped mask with the stimulus
	matR_surround = bsxfun(@times,matR_surround,matRampGrid); %multiply the ramped mask with the stimulus
	
	%% final steps
	%select only required pixels
	vecSelectX = (find(vecSpaceX >= dblVisSpacing*((-(vecArraySize(1) - 1)/2)),1) : find(vecSpaceX >= dblVisSpacing*(((vecArraySize(1) - 1)/2)),1));
	vecSelectY = (find(vecSpaceX >= dblVisSpacing*((-(vecArraySize(2) - 1)/2)),1) : find(vecSpaceX >= dblVisSpacing*(((vecArraySize(2) - 1)/2)),1));
	matR_center = matR_center(vecSelectX,vecSelectY,:); %start contrast-dependent response around 1% instead of 10%, 2019-7-17
	matR_surround = matR_surround(vecSelectX,vecSelectY,:); %start contrast-dependent response around 1% instead of 10%, 2019-7-17
	matR_center = cat(3,zeros([intImX intImY intBaseT]),matR_center,zeros([intImX intImY intBaseT])); %start contrast-dependent response around 1% instead of 10%, 2019-7-17
	matR_surround = cat(3,zeros([intImX intImY intBaseT]),matR_surround,zeros([intImX intImY intBaseT])); %start contrast-dependent response around 1% instead of 10%, 2019-7-17
	
	%shift surround by dblDelayCS
	intShiftTau = round((dblDelayCS/1000)/dblDeltaT);
	matR_surround = cat(3,zeros([vecArraySize intShiftTau]),matR_surround(:,:,1:(end-intShiftTau)));
	
	%luminance offset
	dblBeta = 3;
	vecLuminance = dblBeta * max(0,log10(dblLuminance));
	
	%calculate overall response
	matR_ON = max(0,dblR_baseline + matR_center - matR_surround);
	matR_OFF = max(0,dblR_baseline - matR_center + matR_surround);
	
	matR_ON = bsxfun(@mtimes,matR_ON,vecLuminance);
	matR_OFF = bsxfun(@mtimes,matR_OFF,vecLuminance);
	

%end
