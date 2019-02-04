
%% run Vm sim
runSimBatchCompiled('1,sConn_Col48N2160S637056_2017-09-15.mat,sStim_SquareGratingOri_x2R1_2017-09-15.mat,0,OriNewVm')

%%
clf;
subplot(2,3,1)
imagesc(sData.matVm(:,:,1));
colormap(parula);freezeColors;
title('Vm');colorbar;
colorbar

subplot(2,3,2)
imagesc(sData.matVm(:,:,2));freezeColors;
title('LGN input');colorbar;

subplot(2,3,3)
imagesc(sData.matVm(:,:,3));freezeColors;
title('E input');colorbar;

subplot(2,3,4)
imagesc(sData.matVm(:,:,4));freezeColors;
title('I input');colorbar;


subplot(2,3,5)
matDiff = sData.matVm(:,:,3) - sData.matVm(:,:,4);
dblMax = max(abs(matDiff(:)));
imagesc(matDiff,dblMax*[-1 1]);
colormap(redblue);freezeColors;
title('E-I')
colorbar
%}
%% V1 fields
load('D:\Simulations\Connectivity\sConn_Col48N2160S637056_2017-09-15.mat');

intCellsV1 = sConnectivity.intCellsV1;
vecCellTypes = sConnectivity.vecCellTypes;
cellType = {'Pyramid','Interneuron'};
figure
matFieldsV1 = sConnectivity.matPrefGabors;
vecRand = randperm(intCellsV1,6);
for i=1:6
	subplot(2,3,i)
	intRand = vecRand(i);
	dblMax = max(max(abs(matFieldsV1(:,:,intRand))));
imagesc(matFieldsV1(:,:,intRand),dblMax*[-1 1]);
title(sprintf('V1 neuron %d; %s',intRand,cellType{vecCellTypes(intRand)}));
end
colormap(redblue)

%% V2 fields
intCellsV1 = sConnectivity.intCellsV1;
intCellsV2 = sConnectivity.intCellsV2;
vecCellTypes = sConnectivity.vecCellTypes;
cellType = {'Pyramid','Interneuron'};
figure
matFieldsV2 = sConnectivity.matFieldsV2;
vecRand = randperm(size(matFieldsV2,3),6);
for i=1:6
	subplot(2,3,i)
	intRand = vecRand(i);
	dblMax = max(max(abs(matFieldsV2(:,:,intRand))));
imagesc(matFieldsV2(:,:,intRand),dblMax*[-1 1]);
title(sprintf('V2 neuron %d (#%d); %s',intRand,intRand+intCellsV1,cellType{vecCellTypes(intRand+intCellsV1)}));
end
colormap(redblue)

%% connectivity matrix
intNeurons = sConnectivity.intCortexCells;
vecSynWeight = sConnectivity.vecSynWeight;
vecSynExcInh = sConnectivity.vecSynExcInh;
matSynFromTo = sConnectivity.matSynFromTo;
vecSynConductance = sConnectivity.vecSynConductance;

vecSynWeight = ones(size(vecSynWeight));
vecSign = vecSynExcInh;%vecCortSynType;
vecSign(vecSign==2) = -1;
matConn2D = zeros(intNeurons,intNeurons);
%matConn2D = getFillGrid(matConn2D,matCortConn(:,1),matCortConn(:,2),vecCortConductance.*vecSign);
matConn2D = getFillGrid(matConn2D,matSynFromTo(:,1),matSynFromTo(:,2),vecSynConductance.*vecSign.*vecSynWeight);

vecTypesV1 = vecCellTypes(1:intCellsV1);
[d,vecReorderV1] = sort(vecTypesV1);
vecIdxV2 = (intCellsV1+1):intNeurons;
matConn2D = matConn2D([vecReorderV1 vecIdxV2],[vecReorderV1 vecIdxV2]);


intV1Pyrs = sum(d==1);
matConn2DV2 = matConn2D(1:intV1Pyrs,vecIdxV2);
vecMeanAngleInputs = nan(1,intCellsV2);
for intV2=1:intCellsV2
vecProjFromV1 = find(matConn2DV2(:,intV2)>0);
vecRad = (vecProjFromV1/intV1Pyrs)*2*pi;
vecMeanAngleInputs(intV2) = mod(circ_mean(vecRad),2*pi);
end
[d,vecReorderV2] = sort(vecMeanAngleInputs);
vecReorderV2 = vecReorderV2 + intCellsV1;
vecIdxV1 = 1:intCellsV1;
vecTypesV2 = vecCellTypes(vecReorderV2);
matConn2D = matConn2D([vecIdxV1 vecReorderV2],[vecIdxV1 vecReorderV2]);

[d,vecReReorderV2] = sort(vecTypesV2);
vecReReorderV2 = vecReReorderV2 + intCellsV1;
matConn2D = matConn2D([vecIdxV1 vecReReorderV2],[vecIdxV1 vecReReorderV2]);

imagesc(matConn2D,[-0.6 0.6+eps]);
colormap(redblue);
xlabel('Post-synaptic neuron');
ylabel('Pre-synaptic neuron');
title(['Connectivity matrix']);
cH = colorbar;
cH.Label.String = 'Conductance (G)';
xlim([0 intNeurons]);
ylim([0 intNeurons]);
fixfig;
grid off;
