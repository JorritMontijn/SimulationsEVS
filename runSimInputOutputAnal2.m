%function runSimInputOutputTest
clear all;

%% load connectivity
%if ~exist('sData','var')
%	clearvars;
load('D:\Data\Processed\V1_LIFmodel\F-I_curves.mat');
strConnDir = 'D:\Simulations\Connectivity\';

cellConnFile{1} = 'sConn_Col48N2160S637056_2017-06-14.mat';
cellConnFile{2} = 'sConn_ExcOnlyCol48N2160S377856_2017-06-14.mat';
cellConnFile{3} = 'sConn_InhOnlyCol48N2160S291456_2017-06-14.mat';
cellConnFile{4} = 'sConn_NoRecurCol48N2160S32256_2017-06-14.mat';
hFig = figure;
for intUseConn =1:4
	%intUseConn = 3;
	strConnFile = cellConnFile{intUseConn};
	
	[sConnParams,sData] = loadConnectivity_xArea(strConnDir,strConnFile);
	%end
	
	%% build input
	%set params
	dblDeltaT = 0.0005;
	vecInputG = 400;
	dblDur = 1; %seconds
	dblSynSpikeMem = 0.2;
	
	%get input
	intNeurons = sData.intCortexCells;
	intCellsV1 = sData.intCellsV1;
	%vecThisV = sData.vecCellV_Leak;
	vecThisV = gaussrnd(-56,1,[intNeurons,1]);
	vecGauss = normpdf(1:intCellsV1,intCellsV1/2,intCellsV1/4)';
	vecInput = vecInputG*(vecGauss/max(vecGauss));
	vecInput((intCellsV1+1):intNeurons,:) = 0;
	vecInputIdx = ones(round(dblDur/dblDeltaT),1) * (1:numel(vecInputG));
	vecInputIdx = vecInputIdx(:)';
	vecOverallT = (1:numel(vecInputIdx))*dblDeltaT;
	
	%put in sData
	sData.vecThisV = vecThisV;
	sData.dblSynSpikeMem = dblSynSpikeMem;
	sData.dblDeltaT = dblDeltaT;
	sData.vecOverallT = vecOverallT;
	sData.vecInputIdx = vecInputIdx; %[1 x T] with M index values
	sData.matInput = vecInput; %[N x M]; neurons by input indices
	
	sData.cellSpikeTimesCortex = cell(intNeurons,1);
	sData.vecSpikeCounterCortex = zeros(intNeurons,1);
	sData.intPreAllocationSize = 100;
	%sData.vecCellRefracT=zeros(size(sData.vecCellRefracT));
	
	%% run
	sData = getSimRunNoStim(sData);
	
	%% get spiking data
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	intTrials = numel(vecInputG);
	intNeurons = numel(cellSpikeTimesCortex);
	%vecBinsTime = sort([0 find(diff(vecInputIdx)) numel(vecInputIdx)])*dblDeltaT;
	vecBinsTime = 0:0.001:max(vecOverallT);
	intBins = numel(vecBinsTime)-1;
	%if ~exist('matModelResp','var')
	hTic = tic;
	matModelResp = zeros(intNeurons,intBins);
	
	for intNeuron=1:intNeurons
		vecCounts = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTime);
		matModelResp(intNeuron,:) = vecCounts./diff(vecBinsTime);
		
		if toc(hTic) > 5 || intNeuron==1
			hTic = tic;
			fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
	end
	%end
	vecInhExc = ((sData.vecCellTypes==1)-(sData.vecCellTypes==2))';
	matModelResp = bsxfun(@times,matModelResp,vecInhExc);
	figure;%subplot(1,4,intUseConn);
	imagesc(matModelResp,[-1 1]);colormap(redblue)
	ylim([0 intNeurons]);
	
	xlabel('Time (ms)');
	ylabel('Neuron ID')
	strTag = getFlankedBy(strConnFile,'sConn_','Col');
	if isempty(strTag),strTag = 'FullConn';end
	title(strTag);
	drawnow;
	fixfig
	%% save
	%full screen
	drawnow;
	jFig = get(handle(gcf), 'JavaFrame');
	jFig.setMaximized(true);
	figure(gcf);
	drawnow;
	strFigDir = 'D:\Data\Results\V1_LIFmodel\';
	export_fig([strFigDir 'PopSpiking_' strTag '_' getDate '.tif']);
	export_fig([strFigDir 'PopSpiking_' strTag '_' getDate '.pdf']);
	
	%% analyze
	%vecY1 = fLog(cellFitP{1},vecX);
	
	figure(hFig)
	subplot(2,2,intUseConn)
	vecExc = find(vecInhExc(1:intCellsV1)==1);
	vecInh = find(vecInhExc(1:intCellsV1)==-1);
	plot(1:intCellsV1,0.1*vecInput(1:intCellsV1),'k--')
	hold on
	plot(vecExc,xmean(matModelResp(vecExc,:),2),'r')
	plot(vecInh,-xmean(matModelResp(vecInh,:),2),'b')
	hold off
	title(getFlankedBy(strConnFile,'sConn_','Col'));
	xlabel('Neuron ID')
	ylabel('Mean spiking rate (Hz)')
	fixfig;
	drawnow;
end
%% save
%full screen
drawnow;
jFig = get(handle(gcf), 'JavaFrame');
jFig.setMaximized(true);
figure(gcf);
drawnow;
strFigDir = 'D:\Data\Results\V1_LIFmodel\';
export_fig([strFigDir 'PopSpiking_RespProfile_' getDate '.tif']);
export_fig([strFigDir 'PopSpiking_RespProfile_' getDate '.pdf']);
