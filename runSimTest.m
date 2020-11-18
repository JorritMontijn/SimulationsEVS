%% load data
strConnOld = 'Conn256N1200_2020-10-07.mat'; %old
strConnNew = 'Conn256N1200_2020-10-26.mat'; %new
strStimOld = 'Ret256Noise0.0Ori5_x2R1_2020-07-17.mat'; %old
strStimNew = 'Ret256Noise0.0Ori5_x2R1_2020-07-17.mat'; %new
load(['F:\Code\Simulations\SimulationsEVS\Connectivity\' strConnNew ]);
load(['F:\Code\Simulations\SimulationsEVS\Stimulation\' strStimNew]);

%% plot response of example cell
intNeuron=25;
dblPrefSF = sConnectivity.vecPrefSF(intNeuron);
indSyn = sConnectivity.matSynConnOFF_to_Cort(:,2)==intNeuron;
vecSourceOff = sConnectivity.matSynConnOFF_to_Cort(indSyn,1);
vecW = sConnectivity.vecSynWeightOFF_to_Cort(indSyn);
vecC = sConnectivity.vecSynConductanceOFF_to_Cort(indSyn);
matGrid = zeros(sConnParams.vecSizeInput);
[vecRow,vecCol]=ind2sub(sConnParams.vecSizeInput,vecSourceOff);
matGrid = getFillGrid(matGrid,vecRow,vecCol,vecC);
figure
imagesc(matGrid)
colorbar
title(sprintf('Sum: %e',sum(matGrid(:))))
return
%% test

tic
sStimDriveOld = getDynBottomUpInputs(sStimParams,1);
toc

%%
tic
sStimDrivePeter = getDynBottomUpInputs(sStimParams,0);
toc
%%
figure
%for t=190:size(sStimDriveOld.matLGN_ON,3)
vecT = (sStimParams.dblDeltaT:sStimParams.dblDeltaT:(sStimParams.dblStimDur*2))-(sStimParams.dblStimDur/2);
vecIntT = [210 400 500 600 700];
for intPlot=1:5
	t = vecIntT(intPlot);
	
	subplot(3,5,intPlot)
imagesc(sStimDriveOld.matImage(:,:,t));colorbar
title(sprintf('Image; old;t=%.3fs',vecT(t)))

subplot(3,5,intPlot+5)
imagesc(sStimDriveOld.matLGN_ON(:,:,t));colorbar
title(sprintf('ON; old;t=%.3fs',vecT(t)))

subplot(3,5,intPlot+10)
imagesc(sStimDrivePeter.matLGN_ON(:,:,t));colorbar
title(sprintf('ON; Peter;t=%.3fs',vecT(t)))
drawnow
end

%%
figure
subplot(2,1,1)
plot(squeeze(sStimDriveOld.matLGN_ON(round(end/2),round(end/2),:)))

subplot(2,1,2)
plot(squeeze(sStimDrivePeter.matLGN_ON(round(end/2),round(end/2),:)))
sprintf('done')
