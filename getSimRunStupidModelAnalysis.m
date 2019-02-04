%for intRunFile=1:6
%	for boolDoCorr = [true false]
%% set file & parameters
close all;
clearvars -except intRunFile boolDoCorr;
intRunFile = 3;
boolDoCorr = false;

%% load
if intRunFile == 1
	strFileSim = 'stupidmodelX.mat';
elseif intRunFile == 11 %no connections, hard reset, dT=0.1ms
	strFileSim = 'PIF_Type0_Input30_N5_T1200_Date2017-10-04.mat';
elseif intRunFile == 12 %no connections, soft reset, dT=0.1ms
	strFileSim = 'PIF_Type0_Input30_N5_T12000_Date2017-10-09.mat';
elseif intRunFile == 13 %no connections, soft reset, dT=0.01ms
	strFileSim = 'PIF_Type0_Input30_N5_T12000_Date2017-10-07.mat';
elseif intRunFile == 14  %no connections, soft reset, dT=0.001ms
	strFileSim = 'PIF_Type0_Input30_N5_T1200_Date2017-10-07.mat';
	
elseif intRunFile == 2 %original
	strFileSim = 'stupidmodelWithOtherConnType-1_2017-09-27.mat';
elseif intRunFile == 21 %soft reset, 30
	strFileSim = 'PIF_Type-1_Input30_N250_T1200_Date2017-10-10.mat';
elseif intRunFile == 22 %hard reset, 30
	strFileSim = 'PIF_HardReset_Type-1_Input30_N250_T1200_Date2017-10-10.mat';
elseif intRunFile == 23 %soft reset, 180
	strFileSim = 'PIF_Type-1_Input180_N250_T1200_Date2017-10-10.mat';
elseif intRunFile == 24 %hard reset, 180
	strFileSim = 'PIF_HardReset_Type-1_Input180_N250_T1200_Date2017-10-11.mat';

elseif intRunFile == 3 %stochastic circulant, bump input
	strFileSim = 'PIF_Type21_Input180_N1440_T1200_Date2017-10-12.mat';
	
elseif intRunFile == 4 %step function circulant
	strFileSim = 'stupidmodelWithOtherConnType1_2017-09-29.mat';
elseif intRunFile == 5 %stochastic circulant
	strFileSim = 'stupidmodelWithOtherConnType2_2017-09-29.mat';
elseif intRunFile == 6 %no noise
	%strFileSim = 'stupidmodelWithOtherConnType3_2017-09-26.mat';
	strFileSim = 'stupidmodelWithOtherConnType3_2017-10-05.mat';
	

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
	
end
if ~exist('vecCellTypes','var')
	dblFracE = 0.8;
	intNeurons = size(matJ,1);
	intPyramids = floor(dblFracE*intNeurons);
	intInterneurons = intNeurons-intPyramids;
	vecCellTypes = ones(1,intNeurons);
	vecCellTypes(intPyramids+1:end) = 2;
end

strExtra = '';
if strcmp(getFlankedBy(strFileSim,'PIF_','_Type'),'HardReset')
	strExtra = 'HardReset';
end

%% get spiking data
%cellSpikeTimesCortex
dblSimT = vecOverallT(end);
dblDeltaT = sParameters.dblDeltaT;
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

for intNeuron=1:intNeurons
	cellSpikeTimes{intNeuron} = cellSpikeTimes{intNeuron}(~isnan(cellSpikeTimes{intNeuron}));
	vecCounts = histcounts(cellSpikeTimes{intNeuron},vecBinsTime)/mean(diff(vecBinsTime));
	matModelResp(intNeuron,:) = vecCounts;
	
	if toc(hTic) > 5 || intNeuron==1
		hTic = tic;
		fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
	end
end
matModelResp(:,1) = [];
intTrials = size(matModelResp,2);
matCovariance = cov(matModelResp');
matCorrFacDiag = diag(1./sqrt(diag(matCovariance)));
matCorr = matCorrFacDiag * matCovariance * matCorrFacDiag;

%% get anal mat
dblVarIndep = 76.5;
dblVarShared = 3.5;
if intRunFile == 6
	dblVarIndep = 0;
	dblVarShared = 0;
end

matJ = matJ + eye(size(matJ))*-dblV_thresh; %add reset to diagonal
matCovAnal = getCovAnalFromJ(matJ, dblVarIndep, dblVarShared)/dblStep;

matCorrAnalFacDiag = diag(1./sqrt(diag(matCovAnal)));
matCorrAnal = matCorrAnalFacDiag * matCovAnal * matCorrAnalFacDiag;

%{
%% plot original matrices
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

%export_fig(['D:\Data\Results\NonLeaky\RandomRunType ' num2str(intType) '_' num2str(dblSimT) 'Secs_Conn' strConn '_Date' getDate() '.pdf'])
%}

%% run asymptote analysis
% switch matrices for correlation; cov = corr
if boolDoCorr
	strMetric = 'Corr';
	matCovAnal = matCorrAnal;
	matCovariance = matCorr;
else
	strMetric = 'Cov';
end

%build selection matrices
matSelectPyrPyr = logical(double(vecCellTypes'==1) * double(vecCellTypes==1));
matSelectPyrInt = logical(double(vecCellTypes'==1) * double(vecCellTypes==2) + double(vecCellTypes'==2) * double(vecCellTypes==1));
matSelectIntInt = logical(double(vecCellTypes'==2) * double(vecCellTypes==2));
matSelectLowerTri = tril(true(numel(vecCellTypes)),-1);
matSelectDiag = diag(diag(true(numel(vecCellTypes))));

%keep only spiking cells
vecSpikeCounter = cellfun(@numel,cellSpikeTimes)';
indKeepNeurons = vecSpikeCounter > dblSimT; %keep cells spiking >1 Hz
matSelectSpNr = logical(double(indKeepNeurons') * double(indKeepNeurons));

%get values
vecPyrPyr = matCovariance(matSelectPyrPyr & matSelectLowerTri & matSelectSpNr);
vecPyrInt = matCovariance(matSelectPyrInt & matSelectLowerTri & matSelectSpNr);
vecIntInt = matCovariance(matSelectIntInt & matSelectLowerTri & matSelectSpNr);
vecPyrDiag = matCovariance(matSelectPyrPyr & matSelectDiag & matSelectSpNr);
vecIntDiag = matCovariance(matSelectIntInt & matSelectDiag & matSelectSpNr);

%% calculate subsamples
matModelRespTransp = matModelResp';
vecSubsample = unique(max(2,round((intTrials/10):(intTrials/10):intTrials)));
intNumRandSubSample = 100;
matSaveVals = nan(numel(vecSubsample),intNumRandSubSample,5);
for intSubIdx=1:numel(vecSubsample)
	dblSubVal = vecSubsample(intSubIdx)
	
	for intRandNr = 1:intNumRandSubSample
		vecRandSubSample = randi(intTrials,[1 dblSubVal]);
		
		if boolDoCorr
			matCov = corr(matModelRespTransp(vecRandSubSample,:));
		else
			matCov = cov(matModelRespTransp(vecRandSubSample,:));
		end
		
		vecPyrPyr = matCov(matSelectPyrPyr & matSelectLowerTri & matSelectSpNr);
		vecPyrInt = matCov(matSelectPyrInt & matSelectLowerTri & matSelectSpNr);
		vecIntInt = matCov(matSelectIntInt & matSelectLowerTri & matSelectSpNr);
		vecPyrDiag = matCov(matSelectPyrPyr & matSelectDiag & matSelectSpNr);
		vecIntDiag = matCov(matSelectIntInt & matSelectDiag & matSelectSpNr);
		
		matSaveVals(intSubIdx,intRandNr,1) = nanmean(vecPyrPyr);
		matSaveVals(intSubIdx,intRandNr,2) = nanmean(vecPyrInt);
		matSaveVals(intSubIdx,intRandNr,3) = nanmean(vecIntInt);
		matSaveVals(intSubIdx,intRandNr,4) = nanmean(vecPyrDiag);
		matSaveVals(intSubIdx,intRandNr,5) = nanmean(vecIntDiag);
	end
end
%% get analytic preds
vecAnalCovPyrPyr = matCovAnal(matSelectPyrPyr & matSelectLowerTri & matSelectSpNr);
dblRangePyrPyr = range(vecAnalCovPyrPyr);
dblAnalMeanPyrPyr = mean(vecAnalCovPyrPyr);

vecAnalCovPyrInt = matCovAnal(matSelectPyrInt & matSelectLowerTri & matSelectSpNr);
dblRangePyrInt = range(vecAnalCovPyrInt);
dblAnalMeanPyrInt = mean(vecAnalCovPyrInt);

vecAnalCovIntInt = matCovAnal(matSelectIntInt & matSelectLowerTri & matSelectSpNr);
dblRangeIntInt = range(vecAnalCovIntInt);
dblAnalMeanIntInt = mean(vecAnalCovIntInt);

vecAnalCovPyrDiag = matCovAnal(matSelectPyrPyr & matSelectDiag & matSelectSpNr);
dblRangePyrDiag = range(vecAnalCovPyrDiag);
dblAnalMeanPyrDiag = mean(vecAnalCovPyrDiag);

vecAnalCovIntDiag = matCovAnal(matSelectIntInt & matSelectDiag & matSelectSpNr);
dblRangeIntDiag = range(vecAnalCovIntDiag);
dblAnalMeanIntDiag = mean(vecAnalCovIntDiag);

%% calc means & plot
figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

subplot(2,3,1)
vecMeanCovPyrPyr = xmean(matSaveVals(:,:,1),2);
vecSDCovPyrPyr = 2*xstd(matSaveVals(:,:,1),2);
hold on
plot([1 intTrials],dblAnalMeanPyrPyr*[1 1],'r--')
errorbar(vecSubsample,vecMeanCovPyrPyr,vecSDCovPyrPyr);
hold off
title(sprintf('%s %s, Pyr-Pyr; %.3f (mean) +/- %.3f (2sd)',strExtra,strMetric,vecMeanCovPyrPyr(end),vecSDCovPyrPyr(end)))
xlabel('Time bins used (subsamples)');
ylabel(sprintf('%s, mean+/-2sd over bootstraps (%d)',strMetric,intNumRandSubSample));
xlim([0 intTrials]);
ylim([min([0 min(get(gca,'ylim'))]) max(get(gca,'ylim'))]);

subplot(2,3,2)
vecMeanCovPyrInt = xmean(matSaveVals(:,:,2),2);
vecSDCovPyrInt = 2*xstd(matSaveVals(:,:,2),2);
hold on
plot([1 intTrials],dblAnalMeanPyrInt*[1 1],'r--')
errorbar(vecSubsample,vecMeanCovPyrInt,vecSDCovPyrInt);
hold off
title(sprintf('%s, Pyr-Int; %.3f (mean) +/- %.3f (2sd)',strMetric,vecMeanCovPyrInt(end),vecSDCovPyrInt(end)))
xlabel('Time bins used (subsamples)');
ylabel(sprintf('%s, mean+/-2sd over bootstraps (%d)',strMetric,intNumRandSubSample));
xlim([0 intTrials]);
ylim([min([0 min(get(gca,'ylim'))]) max(get(gca,'ylim'))]);

subplot(2,3,3)
vecMeanCovIntInt = xmean(matSaveVals(:,:,3),2);
vecSDCovIntInt = 2*xstd(matSaveVals(:,:,3),2);
hold on
plot([1 intTrials],dblAnalMeanIntInt*[1 1],'r--')
errorbar(vecSubsample,vecMeanCovIntInt,vecSDCovIntInt);
hold off
title(sprintf('%s, Int-Int; %.3f (mean) +/- %.3f (2sd)',strMetric,vecMeanCovIntInt(end),vecSDCovIntInt(end)))
xlabel('Time bins used (subsamples)');
ylabel(sprintf('%s, mean+/-2sd over bootstraps (%d)',strMetric,intNumRandSubSample));
xlim([0 intTrials]);
ylim([min([0 min(get(gca,'ylim'))]) max(get(gca,'ylim'))]);

subplot(2,3,4)
vecMeanCovPyrDiag = xmean(matSaveVals(:,:,4),2);
vecSDCovPyrDiag = 2*xstd(matSaveVals(:,:,4),2);
hold on
plot([1 intTrials],dblAnalMeanPyrDiag*[1 1],'r--')
errorbar(vecSubsample,vecMeanCovPyrDiag,vecSDCovPyrDiag);
hold off
title(sprintf('%s, Pyr-Diag; %.3f (mean) +/- %.3f (2sd)',strMetric,vecMeanCovPyrDiag(end),vecSDCovPyrDiag(end)))
xlabel('Time bins used (subsamples)');
ylabel(sprintf('%s, mean+/-2sd over bootstraps (%d)',strMetric,intNumRandSubSample));
xlim([0 intTrials]);
ylim([min([0 min(get(gca,'ylim'))]) max(get(gca,'ylim'))]);

subplot(2,3,5)
vecMeanCovIntDiag = xmean(matSaveVals(:,:,5),2);
vecSDCovIntDiag = 2*xstd(matSaveVals(:,:,5),2);
hold on
plot([1 intTrials],dblAnalMeanIntDiag*[1 1],'r--')
errorbar(vecSubsample,vecMeanCovIntDiag,vecSDCovIntDiag);
hold off
title(sprintf('%s, Int-Diag; %.3f (mean) +/- %.3f (2sd)',strMetric,vecMeanCovIntDiag(end),vecSDCovIntDiag(end)))
xlabel('Time bins used (subsamples)');
ylabel(sprintf('%s, mean+/-2sd over bootstraps (%d)',strMetric,intNumRandSubSample));
xlim([0 intTrials]);
ylim([min([0 min(get(gca,'ylim'))]) max(get(gca,'ylim'))]);

subplot(2,3,6)
matJ_Plot = matJ;
matJ_Plot(diag(diag(true(size(matJ))))) = nan;
imagesc(matJ_Plot);colorbar;
title(sprintf('Bin=%.0fs; Type=%d; Input G=%.0f; Neurons=%d; Sim Time=%.0fs; dT=%.3fms',dblStep,intType,dblInputG,intNeurons,dblSimT,dblDeltaT*1000));
xlabel('Neuron');
ylabel('Neuron');
drawnow;

%% save figure
strFigFile = sprintf('%s%s_Type%d_Input%d_N%d_T%d_Bins%d_Date%s',strMetric,strExtra,intType,round(dblInputG),intNeurons,round(dblSimT),intTrials,getDate);
export_fig([strDir strFigFile '.pdf'])

%%
figure;
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;

vecMean = mean(matModelResp,2);
vecVar = var(matModelResp,[],2);

subplot(2,3,1)
scatter(vecMean,vecVar);
set(gca,'xscale','log','yscale','log')
hold on
dblMinAx = min([min(get(gca,'xlim')) min(get(gca,'ylim'))]);
dblMaxAx = max([max(get(gca,'xlim')) max(get(gca,'ylim'))]);
plot([dblMinAx dblMaxAx],[dblMinAx dblMaxAx],'k--')
xlabel('Mean rate')
ylabel('Var rate')
title(sprintf('Bin=%.0fs; Type=%d; Input G=%.0f; Neurons=%d; Sim Time=%.0fs; dT=%.3fms',dblStep,intType,dblInputG,intNeurons,dblSimT,dblDeltaT*1000));

subplot(2,3,2)
vecFano = vecVar ./ vecMean;
histx(vecFano)
ylabel('Number of neurons')
xlabel('Fano factors')
drawnow;
intMaxT = min([10000 size(matModelResp,2)]);

subplot(2,3,3)
imagesc(matCovariance)
xlabel('Neuron');
ylabel('Neuron')
title('Covariance');
colorbar

subplot(2,1,2)
imagesc([0 intMaxT*dblStep],[1 intNeurons],matModelResp(:,1:intMaxT))
xlabel('Time (s)');
ylabel('Neuron')
title(sprintf('Firing rate (Hz)'));
colorbar

strFigFile = sprintf('Fano%s_Type%d_Input%d_N%d_T%d_Bins%d_Date%s',strExtra,intType,round(dblInputG),intNeurons,round(dblSimT),intTrials,getDate);
export_fig([strDir strFigFile '.pdf'])

%	end
%end


