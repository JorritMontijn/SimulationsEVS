%function runSimInputOutputTest
clear all;

%% load full simulation
strSimulation = 'xAreaDistributed_493414_2017-05-04'; %contrast
runModelHeader;

%% get mean spiking response across pref oris
intUseStim = 21;%max(vecTrialStimType);
indUseTrials = vecTrialStimType==intUseStim;
dblStimDur = mean(vecStimStopSecs-vecStimStartSecs);
matActHz = double(matModelResp(:,indUseTrials))/dblStimDur;
vecPyrs = vecCellTypes==1;
vecCellsV1 = vecCellArea==1;
vecMeanV1 = xmean(matActHz(vecCellsV1,:),2);
vecSDV1P = xstd(matActHz(vecPyrs&vecCellsV1,:),2);
dblStimOri = vecTrialOris(intUseStim);
vecPrefOriV1P = mod(rad2ang(vecPrefOri(vecPyrs&vecCellsV1)) - vecTrialOris(intUseStim),180);



%% find parameters of input that recover similar profile
vecRespSim = xmean(matActHz(vecCellsV1,:),2);
intRot = find(vecPrefOri(vecPyrs&vecCellsV1)>=ang2rad(dblStimOri),1);
vecRespSimCentered = circshift(vecRespSim,-intRot-1);
vecCellTypesCentered = circshift(vecCellTypes',-intRot-1);

%% load connectivity
%if ~exist('sData','var')
%	clearvars;
%load('D:\Data\Processed\V1_LIFmodel\F-I_curves.mat');
strConnDir = 'D:\Simulations\Connectivity\';


strConnFile = 'sConn_Col48N1440S403200_2017-05-23.mat';

global sConnIO;
[sConnParams,sConnIO] = loadConnectivity_xArea(strConnDir,strConnFile);

%end
vecPrefOri = sConnIO.vecPrefOri;

%% run error minimization between full run and direct current
vecParams = [0.3 0.3 0.1 1];

%vecCellTypesCentered = circshift(vecCellTypes(vecCellsV1),-intRot);
sOpt = optimset('DiffMinChange',0.05);
lb = [0 0 0 0];
ub = [inf inf inf inf];
[vecFitP1,resnorm,residual,exitflag] = lsqcurvefit(@getModelIO,vecParams,vecPrefOri,vecRespSimCentered,lb,ub,sOpt);
%G-input: 279.100354; SD: 483.396418; Offset: 6.757550 [13:42:58]
%vecP = [279.100354/1000 483.396418/1000 6.757550/100];
vecY = getModelIO(vecFitP1,vecPrefOri);
%%
vecParams = [0.15 0.4 1.8 0.8]
vecY2 = getModelIO(vecParams,vecPrefOri);

%%
cla;
plot(vecRespSimCentered)
hold on
plot(vecY2)

%% get spiking data
