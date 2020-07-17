function [matLGN_ON,matLGN_OFF,varDeltaSyn] = getDynResponsesLGN(matR_ON,matR_OFF,sParams)
	%getResponsesLGN Produces LGN ON and OFF field responses to a retinal drive
	%   [matLGN_ON,matLGN_OFF,varDeltaSyn] = getResponsesLGN(matR_ON,matR_OFF,sParams)
	%
	%matR_ON and matR_OFF are 3D matrices [x,y,t] containing the retinal
	%drive from getRetinalDrive() given the parameters defined by their
	%fieldnames in sParams. If a particular field is absent from the input,
	%the default parameters are chosen as follows:
	%
	%varDeltaSyn = gaussrnd(mu=3,sd=1);	%random synaptic delays, of size [x,y,2] for ON and OFF fields
	%dblLGN_baseline = 15;				%baseline spiking rate in Hz
	%dblStepT = 0.5;					%Time-step in ms
	%dblContrast = 100;					%stimulus contrast
	%dblK = 1/1000;						spiking probability
	
	%get size
	vecSize = size(matR_ON);
	
	%get default values
	if ~exist('sParams','var'),sParams=struct;end
	if ~isfield(sParams,'varDeltaSyn')
		vecDeltaSynP = [3 1];
	else
		if numel(sParams.varDeltaSyn) == 2
			vecDeltaSynP = sParams.varDeltaSyn;
		elseif size(sParams.varDeltaSyn,1) == vecSize(1) && size(sParams.varDeltaSyn,2) == vecSize(2) && size(sParams.varDeltaSyn,3) == 2
			varDeltaSyn=sParams.varDeltaSyn;
		else
			error([mfilename ':InputNotRecognized'],['Input form for sParams.varDeltaSyn not recognized: [' sprintf('%d ',size(sParams.varDeltaSyn)) ']']);
		end
	end
	if ~exist('varDeltaSyn','var')
		%create matrix with random delays
		varDeltaSyn = normrnd(vecDeltaSynP(1),vecDeltaSynP(2),[vecSize(1:2) 2]);
		
		%remove negative delays
		varDeltaSyn(varDeltaSyn<0) = -varDeltaSyn(varDeltaSyn<0);
	end
	if ~isfield(sParams,'dblLGN_baseline'),dblLGN_baseline = 1;else dblLGN_baseline=sParams.dblLGN_baseline;end
	if ~isfield(sParams,'dblStepT'),dblStepT = 0.5;else dblStepT=sParams.dblStepT;end
	if ~isfield(sParams,'dblLuminance'),dblLuminance = 100;else dblLuminance=sParams.dblLuminance;end
	if ~isfield(sParams,'dblContrast'),dblContrast = 100;else dblContrast=sParams.dblContrast;end
	if ~isfield(sParams,'dblK'),dblK = dblStepT/1000;else dblK=sParams.dblK;end
	
	%create LGN responses
	matFillZeros = round(varDeltaSyn/dblStepT);
	matFillZeros_ON = matFillZeros(:,:,1);
	matFillZeros_OFF = matFillZeros(:,:,2);
	
	%pre-allocate
	matLGN_ON = nan(size(matR_ON));
	matLGN_OFF = nan(size(matR_OFF));
	
	%assign delays and responses
	for i = 1:vecSize(1)
		for j = 1:vecSize(2)
			matLGN_ON(i,j,:) = cat(3,zeros(1,1,matFillZeros_ON(i,j)),matR_ON(i,j,1:(end-matFillZeros_ON(i,j))));
			matLGN_OFF(i,j,:) = cat(3,zeros(1,1,matFillZeros_OFF(i,j)),matR_OFF(i,j,1:(end-matFillZeros_OFF(i,j))));
		end
	end
	
	%remove initial segment
	matLGN_ON(:,:,1:50) = repmat(matLGN_ON(:,:,51),[1 1 50]);
	matLGN_OFF(:,:,1:50) = repmat(matLGN_OFF(:,:,51),[1 1 50]);
	
	%normalize for luminance
	matLGN_ON = matLGN_ON - min(matLGN_ON(:));
	matLGN_OFF = matLGN_OFF - min(matLGN_OFF(:));
	dblOldMax = max(max(matLGN_ON(:)),max(matLGN_OFF(:)));
	dblNewMax = dblLGN_baseline + 40*max(0,log10(dblLuminance));
	matLGN_ON = (matLGN_ON/dblOldMax)*dblNewMax;
	matLGN_OFF = (matLGN_OFF/dblOldMax)*dblNewMax;
	
	%contrast dependency
	dblC01 = (dblContrast/100);
	matLGN_ON = (1-dblC01)*(dblNewMax/2) + dblC01*matLGN_ON;
	matLGN_OFF = (1-dblC01)*(dblNewMax/2) + dblC01*matLGN_OFF;
	
	%{
	%contrast dependency
	vecContrast = zeros(1,1,vecSize(3));
	parfor intT=1:vecSize(3)
		mTemp = matLGN_ON(:,:,intT);
		vecContrast(intT) = range(mTemp(:));
	end
	matLGN_ON(:,:,vecContrast==0)=dblLGN_baseline;
	matLGN_OFF(:,:,vecContrast==0)=dblLGN_baseline;
	%}
	
	%transform to spike rate
	matLGN_ON = dblK*matLGN_ON;
	matLGN_OFF = dblK*matLGN_OFF;
end
