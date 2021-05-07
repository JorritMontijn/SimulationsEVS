%% load data
load(['F:\Data\Results\SimResults\xMaster_Ori2Noise0SquareGrating_000_1537938543_2021-03-23']);

%% plot response of example cell
sSimRun.cellSpikeTimesCortex;
sSimRun.vecStimStartSecs


%%
dblStepSize = (sSimRun.vecOverallT(end)-sSimRun.vecOverallT(1))/(numel(sSimRun.vecOverallT)-1);
cellSpikeTimes = sSimRun.cellSpikeTimesCortex;

vecNrSpikes = cellfun(@numel,sSimRun.cellSpikeTimesCortex);
vecStimOn = sSimRun.vecStimStartSecs;
vecStimOff = sSimRun.vecStimStopSecs;
vecOrientation = sSimRun.vecTrialOris;
dblTrialDur = mean(vecStimOff-vecStimOn);

%%
matSpikeRate = getSpikeCounts(cellSpikeTimes,vecStimOn,vecStimOff)./(vecStimOff-vecStimOn);
vecMeanR = mean(matSpikeRate,2);
[matRespNSR,vecStimTypes,vecUnique] = getStimulusResponses(matSpikeRate,vecOrientation);
matMeanNS = mean(matRespNSR(:,:,2),3);
dblThreshold = 8;%9.5
vecMaxResp = max(matMeanNS,[],2); %in any orientation
vecEI = sParams.sConnectivity.vecCellTypes;

figure
plot(matMeanNS)

%% prep
vecPrefOri = sParams.sConnectivity.vecPrefOri*2;
vecUniques = unique(roundi(diff(sort(vecPrefOri)),6));
dblStep = vecUniques(2);
vecBinEdges = (-pi-(dblStep/2)):dblStep:(pi+(dblStep/2));
vecBinCenters = vecBinEdges(2:end)-(dblStep/2);
matSynFromTo = sParams.sConnectivity.matSynFromTo;
%% plot
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	for intSourceEI=1:2
		if intSourceEI==1
			strSourceEI = 'E';
		else
			strSourceEI = 'I';
		end
		
		vecSourceN = find(vecEI(:)==intSourceEI);
		vecTargetN = find(vecEI(:)==intTargetEI);
		
		vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
		vecSynW = sParams.sConnectivity.vecSynWeight(vecPlotSyn);
		vecSynOriDiff = nan(size(vecSynW));
		for intSynIdx=1:numel(vecSynW)
			vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
			vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
		end
		
		intSubplot = intSubplot + 1;
		subplot(2,2,intSubplot);
		scatter(rad2deg(vecSynOriDiff)/2,vecSynW,20,'.');
		title(sprintf('%s to %s',strSourceEI,strTargetEI));
		xlabel('Orientation difference (degs)');
		ylabel('Synaptic weight');
		ylim([0 20]);
		xlim([-90 90]);
		fixfig;grid off;
	end
end

maxfig(gcf,0.85);
drawnow;
%strFigFile1 = ['SynWeightByPrefOriDiff_' strFigBase];
%export_fig(fullfile(strSourcePath,[strFigFile1 '.tif']));
%export_fig(fullfile(strSourcePath,[strFigFile1 '.jpg']));
%export_fig(fullfile(strSourcePath,[strFigFile1 '.pdf']));

%% plot 2
figure
intSubplot = 0;
for intTargetEI=1:2
	if intTargetEI==1
		strTargetEI = 'E';
	else
		strTargetEI = 'I';
	end
	
	
	for intSourceEI=1:2
		if intSourceEI==1
			strSourceEI = 'E';
		else
			strSourceEI = 'I';
		end
		
		vecSourceN = find(vecEI(:)==intSourceEI);
		vecTargetN = find(vecEI(:)==intTargetEI);
		
		vecPlotSyn = find(ismember(matSynFromTo(:,1),vecSourceN) & ismember(matSynFromTo(:,2),vecTargetN));
		vecSynW = sParams.sConnectivity.vecSynWeight(vecPlotSyn);
		vecSynOriDiff = nan(size(vecSynW));
		for intSynIdx=1:numel(vecSynW)
			vecST = matSynFromTo(vecPlotSyn(intSynIdx),:);
			vecSynOriDiff(intSynIdx) = circ_dist(vecPrefOri(vecST(1)),vecPrefOri(vecST(2)));
		end
		
		vecCounts = histcounts(vecSynOriDiff,vecBinEdges);
		intSubplot = intSubplot + 1;
		
		subplot(2,2,intSubplot);
		
		plot(rad2deg(vecBinCenters)/2,vecCounts);
		title(sprintf('%s to %s',strSourceEI,strTargetEI));
		xlabel('Orientation difference (degs)');
		ylabel('# of synapses');
		%ylim([0 20]);
		xlim([-90 90]);
		fixfig;grid off;
	end
end
maxfig(gcf,0.85);
drawnow;
strFigFile2 = ['SynNumberByPrefOriDiff_' strFigBase];
export_fig(fullfile(strSourcePath,[strFigFile2 '.tif']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.jpg']));
export_fig(fullfile(strSourcePath,[strFigFile2 '.pdf']));
