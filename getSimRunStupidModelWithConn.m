%{
function [dblCurT,... %only output
		vecThisV,boolStimPresent,intPrevTrial,intTrialT,intIter,... %input & output; line 1
		cellSpikeTimesLGN_ON,cellSpikeTimesLGN_OFF,cellSpikeTimesCortex,... %input & output; line 2
		vecSpikeCounterLGN_ON,vecSpikeCounterLGN_OFF,vecSpikeCounterCortex,... %input & output; line 3
		vecSpikeCounterPreAllocatedLGN_ON,vecSpikeCounterPreAllocatedLGN_OFF,vecSpikeCounterPreAllocatedCortex,intPreAllocationSize]... %input & output; line 4
		= getSimulationRun(vecOverallT,dblDeltaT,matCortConn,dblSynSpikeMem,vecCortSynType,intCortexCells,vecCortDelay,vecCortConductance,...%only input
		vecCellThresh,vecTauPeakByType,vecCellV_E,vecCellV_I,vecCellV_AHP,vecCellV_Leak,vecCellCm,vecCellG_Leak,vecCellG_AHP,...
		vecSynConductanceON_to_Cort,vecSynConductanceOFF_to_Cort,vecSynWeightON_to_Cort,vecSynWeightOFF_to_Cort,vecSynDelayON_to_Cort,vecSynDelayOFF_to_Cort,...%only input
		matSynConnON_to_Cort,matSynConnOFF_to_Cort,matBlankLGN_ON,matBlankLGN_OFF,cellLGN_ON,cellLGN_OFF,... %only input
		vecTrialOris,vecTrialOriIdx,vecStimStartSecs,vecTrialEndSecs,...  %only input
		vecThisV,boolStimPresent,intPrevTrial,intTrialT,intIter,...  %input & output; line 1
		cellSpikeTimesLGN_ON,cellSpikeTimesLGN_OFF,cellSpikeTimesCortex,... %input & output; line 2
		vecSpikeCounterLGN_ON,vecSpikeCounterLGN_OFF,vecSpikeCounterCortex,... %input & output; line 3
		vecSpikeCounterPreAllocatedLGN_ON,vecSpikeCounterPreAllocatedLGN_OFF,vecSpikeCounterPreAllocatedCortex,intPreAllocationSize) %input & output; line 4
%}

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%function sData = getSimRunStupidModel(sData)
%check if we're on cluster
%	global boolClust;
%	if isempty(boolClust),boolClust=false;end
clear all
hTic = tic;
boolSaveVm = false;

%% load other matJ
runPlotMatJ();
matJ = matJ';
dblV_thresh_preset = 1;
dblNoiseFactor = 1;
dblFracE = nan;

%% edit matJ
intType = 21;
if intType == -1 %original
	intNeurons = 250;
	
	dblFracE = 0.8;
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
	
	dblConnEE = 6/(intPyramids-1); %to E from E
	dblConnIE = 5.9/(intPyramids); %to I from E
	dblConnEI = -9.5/(intInterneurons); %to E from I
	dblConnII = -9.4/(intInterneurons-1); %to I from I
	
	
	matJ = zeros(intNeurons,intNeurons);
	matJ(vecCellTypes==1,vecCellTypes==1) = dblConnEE; %TO pyramid FROM pyramid
	matJ(vecCellTypes==2,vecCellTypes==1) = dblConnIE; %TO interneuron FROM pyramid
	matJ(vecCellTypes==1,vecCellTypes==2) = dblConnEI; %TO pyramid FROM interneuron
	matJ(vecCellTypes==2,vecCellTypes==2) = dblConnII; %TO interneuron FROM interneuron
	matJ(logical(eye(intNeurons))) = 0;
elseif intType == 0 %no connections
	intNeurons = 5;
	dblFracE = 0.8;
	
	
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
	
	dblConnEE = 6/(intPyramids-1); %to E from E
	dblConnIE = 5.9/(intPyramids); %to I from E
	dblConnEI = -9.5/(intInterneurons); %to E from I
	dblConnII = -9.4/(intInterneurons-1); %to I from I
	
	
	matJ = zeros(intNeurons,intNeurons);
	matJ(vecCellTypes==1,vecCellTypes==1) = dblConnEE; %TO pyramid FROM pyramid
	matJ(vecCellTypes==2,vecCellTypes==1) = dblConnIE; %TO interneuron FROM pyramid
	matJ(vecCellTypes==1,vecCellTypes==2) = dblConnEI; %TO pyramid FROM interneuron
	matJ(vecCellTypes==2,vecCellTypes==2) = dblConnII; %TO interneuron FROM interneuron
	matJ(logical(eye(intNeurons))) = 0;
	
	matJ(:) = 0;
elseif intType == 1 %step-function connectivity; reset=5
	intPyrWidth = 160;
	intInhWidth = 120;
	matJ=getBandJ(intPyrWidth,intInhWidth);
	dblV_thresh_preset = 5;
	
elseif intType == 2 || intType == 21 %hybrid; reset=5
	dblV_thresh_preset = 5;
	matIndEE = logical(double(vecCellTypes==1)' * double(vecCellTypes==1)); %to E from E
	matIndIE = logical(double(vecCellTypes==2)' * double(vecCellTypes==1)); %to I from E
	matIndEI = logical(double(vecCellTypes==1)' * double(vecCellTypes==2)); %to E from I
	matIndII = logical(double(vecCellTypes==2)' * double(vecCellTypes==2)); %to I from I
	matIndIsConn = matJ~=0;
	
	dblConnEE = 6/sConnParams.matConnCortFromTo(1,1); %to E from E
	dblConnIE = 5.9/sConnParams.matConnCortFromTo(1,2); %to I from E
	dblConnEI = -9.5/sConnParams.matConnCortFromTo(2,1); %to E from I
	dblConnII = -9.4/sConnParams.matConnCortFromTo(2,2); %to I from I
	matJ(matIndEE & matIndIsConn) = dblConnEE;
	matJ(matIndIE & matIndIsConn) = dblConnIE;
	matJ(matIndEI & matIndIsConn) = dblConnEI;
	matJ(matIndII & matIndIsConn) = dblConnII;
	
	%resort
	[dummy,vecReorder] = sort(vecCellTypes);
	matJ = matJ(vecReorder,vecReorder);
	
elseif intType == 3 %no noise
	dblNoiseFactor = 0;
	
	dblFracE = 0.8;
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
	
	dblConnEE = 6/(intPyramids-1); %to E from E
	dblConnIE = 5.9/(intPyramids); %to I from E
	dblConnEI = -9.5/(intInterneurons); %to E from I
	dblConnII = -9.4/(intInterneurons-1); %to I from I
	
	
	matJ = zeros(intNeurons,intNeurons);
	matJ(vecCellTypes==1,vecCellTypes==1) = dblConnEE; %TO pyramid FROM pyramid
	matJ(vecCellTypes==2,vecCellTypes==1) = dblConnIE; %TO interneuron FROM pyramid
	matJ(vecCellTypes==1,vecCellTypes==2) = dblConnEI; %TO pyramid FROM interneuron
	matJ(vecCellTypes==2,vecCellTypes==2) = dblConnII; %TO interneuron FROM interneuron
	matJ(logical(eye(intNeurons))) = 0;
	
elseif intType > 40 && intType < 50
	if intType == 41 %reproduce ruben's graph, pop1
		intNeurons = 75;
	elseif intType == 42 %reproduce ruben's graph; pop2
		intNeurons = 125;
	elseif intType == 43 %reproduce ruben's graph; pop3
		intNeurons = 250;
	elseif intType == 44 %reproduce ruben's graph; pop4
		intNeurons = 500;	
	end
	
	dblFracE = 0.8;
	
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
	
	dblConnEE = (3*4/2)/(intPyramids-1); %to E from E
	dblConnIE = 5.4440/(intPyramids); %to I from E
	dblConnEI = -9.5049/(intInterneurons); %to E from I
	dblConnII = -8.9456/(intInterneurons-1); %to I from I
	
	
	matJ = zeros(intNeurons,intNeurons);
	matJ(vecCellTypes==1,vecCellTypes==1) = dblConnEE; %TO pyramid FROM pyramid
	matJ(vecCellTypes==2,vecCellTypes==1) = dblConnIE; %TO interneuron FROM pyramid
	matJ(vecCellTypes==1,vecCellTypes==2) = dblConnEI; %TO pyramid FROM interneuron
	matJ(vecCellTypes==2,vecCellTypes==2) = dblConnII; %TO interneuron FROM interneuron
	matJ(logical(eye(intNeurons))) = 0;

end
%% set parameters
%dblInputG = 30; %Hz??
dblInputG = 180; %Hz; 180!!
dblV_reset = 0;
dblV_thresh = dblV_thresh_preset;
dblDeltaT = 0.1/1000;
dblRootDeltaT = sqrt(dblDeltaT);
dblSigmaInd = dblNoiseFactor*sqrt(76.5)*dblRootDeltaT;
dblSigmaS = dblNoiseFactor*sqrt(3.5)*dblRootDeltaT;
dblTauE = 2/1000;
dblTauI = 3/1000;
if intType == 0
	dblSynMem = dblDeltaT;
else
	dblSynMem = 0.05;  %Synaptic memory -- drop spikes beyond this horizon
end

dblSimT = 1200;
%dblSimT = 12000;
vecOverallT = dblDeltaT:dblDeltaT:dblSimT;

%% build model %%% NOTE: J IS TRANSPOSE OF OTHER MODEL: TO FROM
vecCellThresh = dblV_thresh*ones(intNeurons,1);
vecCellReset = dblV_reset*ones(intNeurons,1);
vecCellResetBy = vecCellReset - vecCellThresh;

if intType == 21
	vecRads = (-pi + pi/intNeurons):(2*pi)/(intNeurons):pi;
	vecInput = dblInputG*dblDeltaT*vonMises(vecRads,[0 1])';
else
	vecInput = (ones(intNeurons,1)*dblInputG)*dblDeltaT;
end

%PSP function
dblSumE = sum(max(0,dblSynMem - (0:dblDeltaT:dblSynMem)).*(1/dblTauE).*exp(-(dblSynMem - (0:dblDeltaT:dblSynMem))/dblTauE),2);
dblSumI = sum(max(0,dblSynMem - (0:dblDeltaT:dblSynMem)).*(1/dblTauI).*exp(-(dblSynMem - (0:dblDeltaT:dblSynMem))/dblTauI),2);
vecTau = [dblTauE dblTauI];
vecSum = [dblSumE dblSumI];
fPSP = @(dblT,vecSpikeTimes,intSynType) sum(max(0,dblT - vecSpikeTimes).*(1/vecTau(intSynType)).*exp(-(dblT - vecSpikeTimes)/vecTau(intSynType)),2)/vecSum(intSynType);

%pre-allocate
vecThisV = zeros(intNeurons,1);
cellSpikeTimes = cell(intNeurons,1);
vecSpikeCounter = zeros(intNeurons,1);
intPreAllocationSize = 100;
vecSpikeCounterPreAllocated = vecSpikeCounter;

%% create recent spike array Cortex
%build cylindrical recent spike storage
intMaxRecentSpikes = 100;
matRecentSpikes = nan(intNeurons,intMaxRecentSpikes);
vecRSCounter = zeros(1,intNeurons);
vecRSPos = ones(1,intNeurons);
dblCurT = 0;

%get recent spikes
for intN = 1:intNeurons
	vecRecentSpikes = cellSpikeTimes{intN}(((cellSpikeTimes{intN} - dblCurT) + dblSynMem) > 0);
	if ~isempty(vecRecentSpikes)
		intNumRecentSpikes = numel(vecRecentSpikes);
		matRecentSpikes(intN,vecRSPos(intN):(vecRSPos(intN)+intNumRecentSpikes-1)) = vecRecentSpikes;
		vecRSPos(intN) = vecRSPos(intN) + intNumRecentSpikes;
		vecRSCounter(intN) = vecRSCounter(intN) + intNumRecentSpikes;
	end
end
%vecTimeSinceSpike = (vecOverallT(1) + 10)*ones(size(vecThisV));

%% pre-allocate Vm matrix
if boolSaveVm
	matVm = nan(intNeurons,numel(vecOverallT));
end

% run
intIter = 0;
hTic = tic;
%fprintf('\n   >>> Starting simulation run [%s]\n',getTime);
dblLastMsg = toc(hTic);
for dblCurT = vecOverallT%(end):dblDeltaT:(vecOverallT(end)+dblDeltaT*100)
	intIter = intIter + 1;
	
	%% get inputs
	%intInputIdx = vecInputIdx(intIter); %[1 x T] with M index values
	%vecInputG = matInput(:,intInputIdx); %[N x M]; neurons by input indices
	
	
	%%
	% remove old spikes
	for intN = find(vecRSCounter>0)
		intPos = getCycPos(vecRSPos(intN),vecRSCounter(intN),intMaxRecentSpikes);
		if dblCurT - matRecentSpikes(intN,intPos) > dblSynMem
			matRecentSpikes(intN,intPos) = nan;
			vecRSCounter(intN) = vecRSCounter(intN) - 1;
		end
	end
	
	%%
	%get PSPs for all spiking neurons
	vecPSPs = zeros(intNeurons,1);
	for intSourceNeuron=find(vecRSCounter>0)
		%get recent spikes
		vecPos = getCycPosVec(vecRSPos(intSourceNeuron),vecRSCounter(intSourceNeuron),intMaxRecentSpikes);
		vecSpikeTimes = matRecentSpikes(intSourceNeuron,vecPos);
		vecSpikeTimes(isnan(vecSpikeTimes))=0;
		intSynType = vecCellTypes(intSourceNeuron);
		
		dblSynG = fPSP(dblCurT,vecSpikeTimes,intSynType);
		if isnan(dblSynG),break;end
		
		%calculate PSPs
		vecPSPs = vecPSPs + dblSynG * matJ(:,intSourceNeuron);
	end
	
	% calculate stochastic components
	vecXiInd = max(min(3,normrnd(0,1,[intNeurons,1])),-3);
	dblXiShared = max(min(3,normrnd(0,1,1)),-3);
	
	% integrate all inputs
	vecPrevV = vecThisV;
	vecThisV = vecPrevV + ...
		vecPSPs + ...
		vecInput + ...
		dblSigmaInd*vecXiInd + ...
		dblSigmaS*dblXiShared;
	
	%check for errors
	if any(isnan(vecThisV) | isinf(vecThisV));
		intErrN = find(isnan(vecThisV) | isinf(vecThisV),1);
		
		fprintf('\n.. ERROR DETECTED! Please check; neuron %d; Currents; PSP: %f, Input: %f, Xi_I: %f, Xi_S: %f [%s]\n',intErrN,vecPSPs(intErrN),vecInput(intErrN),dblSigmaInd*vecXiInd(intErrN),dblSigmaS*dblXiShared,getTime);
		warning([mfilename ':VmOutOfBounds'],'Membrane voltage out of bounds!');
		break;
	end
	
	% spiking
	%vecTimeSinceSpike = vecTimeSinceSpike + dblDeltaT;
	indSpiking = vecThisV>=vecCellThresh;
	if any(indSpiking)
		for intN = find(indSpiking)'
			%check number of spikes and pre-allocated size
			vecSpikeCounter(intN) = vecSpikeCounter(intN) + 1;
			if vecSpikeCounterPreAllocated(intN) < vecSpikeCounter(intN)
				vecSpikeCounterPreAllocated(intN) = vecSpikeCounterPreAllocated(intN) + intPreAllocationSize;
				cellSpikeTimes{intN} = [cellSpikeTimes{intN} nan(1,intPreAllocationSize)];
			end
			
			% add spike time
			cellSpikeTimes{intN}(vecSpikeCounter(intN)) = dblCurT;
			vecThisV(intN) = vecThisV(intN) + vecCellResetBy(intN); %soft reset
			%vecThisV(intN) = vecCellReset(intN); %hard reset
			
			%add recent spikes
			matRecentSpikes(intN,vecRSPos(intN)) = dblCurT;
			vecRSPos(intN) = mod(vecRSPos(intN) + 1,intMaxRecentSpikes);
			vecRSCounter(intN) = vecRSCounter(intN) + 1;
		end
	end
	vecRSPos(vecRSPos==0)=intMaxRecentSpikes; %mod(y,y) is 0, but we want it to be y
	
	
	%save Vm
	if boolSaveVm
		matVm(:,intIter) = vecThisV;
	end
	
	%% message
	if (toc(hTic) - dblLastMsg) > 5
		dblLastMsg = toc(hTic);
		dblFracDone = intIter/numel(vecOverallT);
		dblTimeElapsed = toc(hTic);
		dblTimeRemaining = (dblTimeElapsed/dblFracDone)-dblTimeElapsed;
		strTimeDone = datestr(datenum(datetime('now'))+dblTimeRemaining/(24*60*60));
		
		vecTime = fix(clock);
		strTime = sprintf('%02d:%02d:%02d',vecTime(4:6));
		
		fprintf('Now at t=%.1fs / %.1fs [%s]; estimated completion: [%s]\n',dblCurT,vecOverallT(end),strTime,strTimeDone);
	end
end

%end msg
hToc = toc(hTic);
fprintf('Done! Total run took %f seconds [%s]\n',hToc,getTime);


%% get spiking data
%cellSpikeTimesCortex
dblSimT = vecOverallT(end);
%dblSimT = dblCurT;%vecOverallT(end);
dblStep = 1;
vecBinsTime = 0:dblStep:dblSimT;
%vecBinsTime = 0:(2/dblDeltaT):(dblSimT/dblDeltaT);
%vecBinsTime = 0:10:dblSimT;
intTrials = numel(vecBinsTime)-1;
intNeurons = numel(cellSpikeTimes);
%if ~exist('matModelResp','var')
hTic = tic;
matModelResp = zeros(intNeurons,intTrials);
vecSpikeCounter = cellfun(@numel,cellSpikeTimes);
for intNeuron=1:intNeurons
	vecCounts = histcounts(cellSpikeTimes{intNeuron},vecBinsTime)/mean(diff(vecBinsTime));
	matModelResp(intNeuron,:) = vecCounts;
	
	if toc(hTic) > 5 || intNeuron==1
		hTic = tic;
		fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
	end
end
matCovariance = cov(matModelResp');
matCovarianceReduced = matCovariance(vecSpikeCounter > 0, vecSpikeCounter > 0);
matCorr = corr(matModelResp');
matCorr(diag(diag(true(size(matCorr)))))=0;

%% plot
figure
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

subplot(2,2,1)
imagesc(matJ);colorbar
title(sprintf('Total running time: %.1fs',dblSimT));
freezeColors;

subplot(2,2,2)
%matCovarianceReduced = matCovariance(vecSpikeCounter > 0, vecSpikeCounter > 0);
plot(xmean(matModelResp,2));
ylim([0 max(get(gca,'ylim'))]);
title(sprintf('Type: %d; input current G: %.3f',intType,dblInputG));
freezeColors;

subplot(2,2,3)
%matCovarianceReduced = matCovariance(vecSpikeCounter > 0, vecSpikeCounter > 0);
imagesc(matCovariance);colorbar
title(sprintf('Covariance; Time block size: %.3fs',dblStep));
freezeColors;

subplot(2,2,4)
imagesc(matCorr);
set(gca,'CLim',[-1 1]*max(abs(get(gca,'Clim'))));
colormap(redblue);
colorbar
title(sprintf('Correlation; Time block size: %.3fs',dblStep));

strConn = strConnFile(end-13:end-4);
if intType == 0,strConn = 'NoConn';end
export_fig(['D:\Data\Results\NonLeaky\RandomRunType ' num2str(intType) '_' num2str(dblSimT) 'Secs_Conn' strConn '_Date' getDate() '.pdf'])

%% gather parameters
sParameters = struct;
%sParameters.matConductancesFromTo = matConductancesFromTo;
sParameters.intNeurons = intNeurons;
sParameters.dblInputG = dblInputG; %Hz??
sParameters.dblFracE = dblFracE;
sParameters.dblV_reset = dblV_reset;
sParameters.dblV_thresh = dblV_thresh;
sParameters.dblSigmaInd = dblSigmaInd;
sParameters.dblSigmaS = dblSigmaS;
sParameters.dblTauE = dblTauE;
sParameters.dblTauI = dblTauI;
sParameters.dblSynMem = dblSynMem;
sParameters.dblDeltaT = dblDeltaT;
sParameters.dblSimT = dblSimT;
sParameters.dblInputG = dblInputG;
sParameters.intType = intType;


%% save
strFile = sprintf('PIF_Type%d_Input%d_N%d_T%d_Date%s.mat',intType,round(dblInputG),intNeurons,round(dblSimT),getDate);

fprintf('Processing completed; saving data <%s> [%s]\n',strFile,getTime);
%end
save(['D:\Data\Results\NonLeaky\' strFile],'sParameters','cellSpikeTimes','matModelResp','vecCellTypes','matJ','matCovariance', 'vecBinsTime','vecOverallT');
%}

