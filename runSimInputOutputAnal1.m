%function runSimInputOutputTest
clear all;
%% load connectivity
if ~exist('sData','var')
	clearvars;
	strConnDir = 'D:\Simulations\Connectivity\';
	%strConnFile = 'sConn_Col48N2S0_2017-05-18.mat';
	strConnFile = 'sConn_Col48N2S0_2017-05-23.mat';
	
	[sConnParams,sData] = loadConnectivity_xArea(strConnDir,strConnFile);
end

%% build input
%set params
dblDeltaT = 0.0005;
vecInputG = 0:10:1600;
dblDur = 10; %seconds
dblSynSpikeMem = 0.2;

%get input
intNeurons = sData.intCortexCells;
vecThisV = sData.vecCellV_Leak;
matInput = repmat(vecInputG,[intNeurons 1]);
vecInputIdx = ones(round(dblDur/dblDeltaT),1) * (1:numel(vecInputG));
vecInputIdx = vecInputIdx(:)';
vecOverallT = (1:numel(vecInputIdx))*dblDeltaT;

%put in sData
sData.vecThisV = vecThisV;
sData.dblSynSpikeMem = dblSynSpikeMem;
sData.dblDeltaT = dblDeltaT;
sData.vecOverallT = vecOverallT;
sData.vecInputIdx = vecInputIdx; %[1 x T] with M index values
sData.matInput = matInput; %[N x M]; neurons by input indices
sData.boolSaveVm = true;

sData.cellSpikeTimesCortex = cell(intNeurons,1);
sData.vecSpikeCounterCortex = zeros(intNeurons,1);
sData.intPreAllocationSize = 1000;

%% edit
%sData.vecCellG_Leak = [10;25]; %[35;25]
%sData.vecCellG_AHP = [100;50]; %[50;25]
%sData.vecCellCm = [1;0.5]; %[0.5;0.2]
%sData.vecTauPeakByType =  [0.001 0.0001];%[0.001 0.002]
%sData.vecCellTauPeak = [];
sData.vecCellRefracT = [0.008;0.004];

%% run
sData = getSimRunNoStim(sData);

%% get data
matVm = sData.matVm;
cellSpikeTimesCortex = sData.cellSpikeTimesCortex;

%% get spiking data
%cellSpikeTimesCortex
intTrials = numel(vecInputG);
intNeurons = numel(cellSpikeTimesCortex);
vecBinsTime = sort([0 find(diff(vecInputIdx)) numel(vecInputIdx)])*dblDeltaT;
%if ~exist('matModelResp','var')
	hTic = tic;
	matModelResp = zeros(intNeurons,intTrials);
	
	for intNeuron=1:intNeurons
		vecCounts = histcounts(cellSpikeTimesCortex{intNeuron},vecBinsTime);
		matModelResp(intNeuron,:) = vecCounts./diff(vecBinsTime);
		
		if toc(hTic) > 5 || intNeuron==1
			hTic = tic;
			fprintf('Now at neuron %d/%d [%s]\n',intNeuron,intNeurons,getTime);
		end
	end
%end
return
%% plot
figure
subplot(2,2,1)
plot(vecOverallT,matVm')
legend({'Pyramid','Interneuron'});
title('No refractory periods')

subplot(2,2,2)
plot(dblDeltaT*vecInputG,matModelResp');
legend({'Pyramid','Interneuron'},'Location','Best');
xlabel('Input current (nA)')
ylabel('Firing rate (Hz)')
title('Absolute refractory periods; 8ms (Pyr), 4ms (Int)')
ylim([0 300])
fixfig

return
%% save
strFigDir = 'D:\Data\Results\V1_LIFmodel\';
export_fig([strFigDir 'F-I_curves_' getDate '.tif']);
export_fig([strFigDir 'F-I_curves_' getDate '.pdf']);

%% fit logistic growth
fLog = @(vecP,vecX) max(0,real(vecP(2) * log((max(0,vecX+vecP(1)))*vecP(3))));
vecFitX = dblDeltaT*vecInputG;
vecY = matModelResp(1,:);

lb = [-inf 0 0];
ub = [inf inf inf];
[vecFitP1,resnorm,residual,exitflag] = curvefitfun(fLog,[1 1 1],vecFitX,vecY,lb,ub);
vecFitY = fLog(vecFitP1,vecFitX);

figure
hold on;
scatter(vecFitX,vecY,'rx');
plot(vecFitX,vecFitY,'r')

vecY = matModelResp(2,:);

lb = [-inf 0 0];
ub = [inf inf inf];
[vecFitP2,resnorm,residual,exitflag] = curvefitfun(fLog,[1 1 1],vecFitX,vecY,lb,ub);
vecFitY = fLog(vecFitP2,vecFitX);

scatter(vecFitX,vecY,'bx');
plot(vecFitX,vecFitY,'b')

%% save
cellFitP{1} = vecFitP1;
cellFitP{2} = vecFitP2;

save('D:\Data\Processed\V1_LIFmodel\F-I_curves.mat','fLog','cellFitP')
