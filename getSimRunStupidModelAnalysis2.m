%for intRunFile=1:6
%	for boolDoCorr = [true false]
%% set file & parameters
close all;
clearvars -except intRunFile boolDoCorr;
intRunFile = 11;
boolDoCorr = false;

%% load
if intRunFile == 1
	strFileSim = 'stupidmodelX.mat';
elseif intRunFile == 2 %original
	strFileSim = 'stupidmodelWithOtherConnType-1_2017-09-27.mat';
elseif intRunFile == 3 %no connections
	strFileSim = 'stupidmodelWithOtherConnType0_2017-09-26.mat';
elseif intRunFile == 4 %step function circulant
	strFileSim = 'stupidmodelWithOtherConnType1_2017-09-29.mat';
elseif intRunFile == 5 %stochastic circulant G=180
	strFileSim = 'stupidmodelWithOtherConnType2_2017-09-29.mat';
elseif intRunFile == 51 %stochastic circulant G=30
	strFileSim = 'stupidmodelWithOtherConnType2_2017-09-27.mat';
elseif intRunFile == 6 %no noise
	%strFileSim = 'stupidmodelWithOtherConnType3_2017-09-26.mat';
	strFileSim = 'stupidmodelWithOtherConnType3_2017-10-05.mat';
elseif intRunFile == 7 %no connections
	strFileSim = 'PIF_Type0_Input30_N5_T1200_Date2017-10-04.mat';
elseif intRunFile == 8 %no connections, soft reset
	strFileSim = 'PIF_Type0_Input30_N5_T1200_Date2017-10-06.mat';
	
elseif intRunFile == 11  %no connections, soft reset, dT=0.001ms
	strFileSim = 'PIF_Type21_Input30_N1440_T120_Date2017-10-05.mat';
end

strDir = 'D:\Data\Results\NonLeaky\';
%strDir = 'A:\Dropbox\_PopCoding_Conn_Corr_R2Pred\sources\';
intType = str2double(getFlankedBy(strFileSim,'Type','_'));
load([strDir strFileSim]);
if exist('sParameters','var')
	dblInputG = sParameters.dblInputG;
	dblV_thresh = sParameters.dblV_thresh;
else
	dblInputG = nan;
	dblV_thresh = 1;
	intType=-1;
	
	dblFracE = 0.8;
	intNeurons = size(matJ,1);
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
end

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
matCorrFacDiag = diag(1./sqrt(diag(matCovariance)));
matCorr = matCorrFacDiag * matCovariance * matCorrFacDiag;

%% get initial spiking data
%cellSpikeTimesCortex
dblInitialT = 50;
dblStepInitial = 0.1;
vecBinsInitialTime = 0:dblStepInitial:dblInitialT;
intInitialT = numel(vecBinsInitialTime)-1;
hTic = tic;
matInitialResp = zeros(intNeurons,intInitialT);
vecSpikeCounter = cellfun(@numel,cellSpikeTimes);
for intNeuron=1:intNeurons
	vecCounts = histcounts(cellSpikeTimes{intNeuron},vecBinsInitialTime)/mean(diff(vecBinsInitialTime));
	matInitialResp(intNeuron,:) = vecCounts;
	
	if toc(hTic) > 5 || intNeuron==1
		hTic = tic;
		fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
	end
end
matCovInitial = cov(matInitialResp');
matCorrInitialFacDiag = diag(1./sqrt(diag(matCovInitial)));
matCorrInitial = matCorrInitialFacDiag * matCovInitial * matCorrInitialFacDiag;

%% plot
figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

intMaxT = min([10000 size(matModelResp,2)]);
subplot(2,1,1)
imagesc([0 intMaxT*dblStep],[1 intNeurons],matModelResp(:,1:intMaxT))
xlabel('Time (s)');
ylabel('Neuron')
title(sprintf('Firing rate (Hz); Type=%d; Input G=%.1f; Neurons=%d; Sim Time=%.1fs; ',intType,dblInputG,intNeurons,dblSimT));
colorbar

subplot(2,2,3)
imagesc(matCovariance)
xlabel('Neuron');
ylabel('Neuron')
title('Covariance');
colorbar

subplot(2,2,3)
imagesc(matCovariance)
xlabel('Neuron');
ylabel('Neuron')
title('Covariance');
colorbar

subplot(2,2,4)
imagesc([0 dblInitialT],[1 intNeurons],matInitialResp)
xlabel('Time (s)');
ylabel('Neuron')
title(sprintf('Initial %.1fs response (Hz)',dblInitialT));
colorbar

%% save figure
strFigFile = sprintf('PIF_Anal2_Type%d_Input%d_N%d_T%d_Date%s',intType,round(dblInputG),intNeurons,round(dblSimT),getDate);
export_fig([strDir strFigFile '.pdf'])
