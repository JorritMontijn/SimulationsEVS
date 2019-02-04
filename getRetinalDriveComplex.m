function [matR_ON,matR_OFF,dblVisSpacingImage,vecLuminance] = getRetinalDriveComplex(matImage,sParams)
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
	
	%get default values
	if ~exist('sParams','var'),sParams=struct;end
	if ~isfield(sParams,'dblContrast'),dblContrast = 100;else dblContrast=sParams.dblContrast;end
	if ~isfield(sParams,'dblLuminance'),dblLuminance = 100;else dblLuminance=sParams.dblLuminance;end
	if ~isfield(sParams,'dblStepT'),dblStepT = 0.5;else dblStepT=sParams.dblStepT;end
	if ~isfield(sParams,'dblVisSpacingImage'),dblVisSpacingImage = 0.2;else dblVisSpacingImage=sParams.dblVisSpacingImage;end
	if ~isfield(sParams,'vecArraySize'),vecArraySize = [64 64];else vecArraySize=sParams.vecArraySize;end
	if ~isfield(sParams,'dblVisSpacingCells'),dblVisSpacingCells = 0.2;else dblVisSpacingCells=sParams.dblVisSpacingCells;end
	if ~isfield(sParams,'dblDelayCS'),dblDelayCS = 3;else dblDelayCS=sParams.dblDelayCS;end
	if ~isfield(sParams,'dblSigmaCenter'),dblSigmaCenter = 0.176;else dblSigmaCenter=sParams.dblSigmaCenter;end
	if ~isfield(sParams,'dblSigmaSurround'),dblSigmaSurround = 0.53;else dblSigmaSurround=sParams.dblSigmaSurround;end
	if ~isfield(sParams,'dblK_center'),dblK_center = 17;else dblK_center=sParams.dblK_center;end
	if ~isfield(sParams,'dblK_surround'),dblK_surround = 16;else dblK_surround=sParams.dblK_surround;end
	if ~isfield(sParams,'dblTauCenter'),dblTauCenter = 10;else dblTauCenter=sParams.dblTauCenter;end
	if ~isfield(sParams,'dblTauSurround'),dblTauSurround = 20;else dblTauSurround=sParams.dblTauSurround;end
	if ~isfield(sParams,'dblR_baseline'),dblR_baseline = 15;else dblR_baseline=sParams.dblR_baseline;end
	
	%resize input image
	if dblVisSpacingImage ~= dblVisSpacingCells
		matRetinalInput = imresize(matImage,dblVisSpacingImage/dblVisSpacingCells);
		dblVisSpacingImage = dblVisSpacingImage / (dblVisSpacingImage/dblVisSpacingCells);
	else
		matRetinalInput = matImage;
	end
	
	%get image size
	[intImX,intImY,intImT] = size(matRetinalInput);
	
	%get spatial receptive fields
	vecSpaceX = dblVisSpacingImage*((-(intImX - 1)/2):(intImX - 1)/2);
	vecSpaceY = dblVisSpacingImage*((-(intImY - 1)/2):(intImY - 1)/2);
	[matMeshX,matMeshY] = meshgrid(vecSpaceX,vecSpaceY);
	
	%get center RF
	matF_center = (dblK_center / (2*pi*dblSigmaCenter^2) ) * exp( -(matMeshX.^2 + matMeshY.^2) / (2*dblSigmaCenter^2));
	intRangeC = ceil((dblSigmaCenter*3)/dblVisSpacingImage);
	vecSelectC = (-(intRangeC+1):intRangeC) + ceil((intImX - 1)/2) + 1;
	matF_center = matF_center(vecSelectC,vecSelectC);
	
	%get surround RF
	matF_surround = (dblK_surround / (2*pi*dblSigmaSurround^2) ) * exp( -(matMeshX.^2 + matMeshY.^2) / (2*dblSigmaSurround^2));
	intRangeS = ceil((dblSigmaSurround*3)/dblVisSpacingImage);
	vecSelectS = (-(intRangeS+1):intRangeS) + ceil((intImX - 1)/2) + 1;
	matF_surround = matF_surround(vecSelectS,vecSelectS);
	
	%get time vector
	clear vecT;
	vecT(1,1,:) = intImT*dblStepT:-1:1;
	vecG_center = (1/dblTauCenter) * exp(-vecT/dblTauCenter);
	vecG_center(vecG_center < max(vecG_center)/1000) = []; %remove extremely low values
	vecG_surround= (1/dblTauSurround) * exp(-vecT/dblTauSurround);
	vecG_surround(vecG_surround < max(vecG_surround)/1000) = []; %remove extremely low values
	
	%make spatiotemporal filter
	matFilterCenter = bsxfun(@times,matF_center,vecG_center); %(:,:,end) is highest value (i.e., t=now)
	matFilterSurround = bsxfun(@times,matF_surround,vecG_surround); %(:,:,end) is highest value (i.e., t=now)
	
	%convolve
	matRetinalInput = matRetinalInput - 0.5;
	matR_center = nan(size(matRetinalInput));
	matR_surround = nan(size(matRetinalInput));
	intMaxT = size(matRetinalInput,3);
	parfor intT=1:intMaxT
		matR_center(:,:,intT) = convn(matRetinalInput(:,:,intT),matFilterCenter,'same');
		matR_surround(:,:,intT) = convn(matRetinalInput(:,:,intT),matFilterSurround,'same');
		if mod(intT,100)==0,fprintf('Convolving; now at %d/%d [%s]\n',intT,intMaxT,getTime);drawnow;end
	end
	matR_center = matR_center + 0.5;
	matR_surround = matR_surround + 0.5;
	%matRetinalInput = matRetinalInput + 0.5;
	
	%select only required pixels
	vecSelectX = (find(vecSpaceX >= dblVisSpacingCells*((-(vecArraySize(1) - 1)/2)),1) : find(vecSpaceX >= dblVisSpacingCells*(((vecArraySize(1) - 1)/2)),1))-1;
	vecSelectY = (find(vecSpaceX >= dblVisSpacingCells*((-(vecArraySize(2) - 1)/2)),1) : find(vecSpaceX >= dblVisSpacingCells*(((vecArraySize(2) - 1)/2)),1))-1;
	matR_center = matR_center(vecSelectX,vecSelectY,:);
	matR_surround = matR_surround(vecSelectX,vecSelectY,:);
	
	%shift surround by dblDelayCS
	intShiftTau = round(dblDelayCS/dblStepT);
	matR_surround = cat(3,zeros([vecArraySize intShiftTau]),matR_surround(:,:,1:(end-intShiftTau)));
	
	%calculate overall response
	matR_ON = max(0,dblR_baseline + matR_center - matR_surround);
	matR_OFF = max(0,dblR_baseline - matR_center + matR_surround);
	dblBeta = 3;
	
	%contrast dependency
	%{
	vecContrast = zeros(1,1,intMaxT);
	parfor intT=1:intMaxT
		mTemp = matR_ON(:,:,intT);
		vecContrast(intT) = range(mTemp(:));
	end
	vecContrast=(vecContrast./max(vecContrast(:)))*dblContrast;
	vecContrast=dblBeta * max(0,log10(vecContrast));
	%}
	vecLuminance = dblBeta * max(0,log10(dblLuminance));
	matR_ON = bsxfun(@mtimes,matR_ON,vecLuminance);
	matR_OFF = bsxfun(@mtimes,matR_OFF,vecLuminance);
	
	%plot
	%{
	dblMax = max(max(matR_ON(:)),max(matR_OFF(:)));
	matPlot_ON = matR_ON/dblMax;
	matPlot_OFF = matR_OFF/dblMax;
	
	close;
	for intT=1:intImT
		subplot(2,2,1)
		imshow(matRetinalInput(:,:,intT));
		title(sprintf('t=%d',intT))
		subplot(2,2,3)
		imshow(matPlot_ON(:,:,intT));
		subplot(2,2,4)
		imshow(matPlot_OFF(:,:,intT));
		pause
	end
	%}
%end
