%function runSimLIF

%% starting message
clearvars;
fprintf(' \nStarting [xArea] Leaky Integrate and Fire simulation of Retina, LGN and V1 [%s]\n\n',getTime);
%try objParPool = parpool('CPU multicore',10);end%delete(gcp('nocreate'))

%% set parameters
global boolClust;
boolClust = false;
dblDeltaT = 0.5/1000; %seconds
dblSynSpikeMem = 0.2; %synaptic spike memory in seconds; older spikes are ignored when calculating PSPs

%% prepare stimulus list
sStimParams = struct;
sStimParams.strStimType = 'SquareGrating'; %{'SquareGrating','SineGrating','Line'}
sStimParams.dblStartFirstTrialSecs = 0.1;
sStimParams.dblPreStimBlankDur = 0.1;%0.25
sStimParams.dblStimDur = 0.3;%0.5
sStimParams.dblPostStimBlankDur = 0.1;%0.25
sStimParams.vecOrientations = [42.5 47.5];%[0:(180/12):179];%[0:(360/16):359]; [42.5 47.5];
sStimParams.vecSpatialFrequencies = 0.5;%2.^[-4:1];
sStimParams.vecTemporalFrequencies = 3;
sStimParams.intReps = 2; %number of repetitions
sStimParams.dblDeltaT = dblDeltaT;
sStimParams.dblStimSizeRetDeg = 16;
sStimParams.vecScrPixWidthHeight = [128 128];
sStimParams.vecScrDegWidthHeight = [25.6 25.6];
	
%get stimuli
sStimInputs = loadSimStim(sStimParams);

figure
colormap(grey)
matPlotOFF = sStimInputs.cellLGN_OFF{end};
matPlotON = sStimInputs.cellLGN_ON{end};
matPlot = matPlotON - matPlotOFF;
dblMin = min(matPlot(:));
dblMax = max(matPlot(:));
for i=1:10:size(matPlot,3)
	imagesc(matPlot(:,:,i),[dblMin dblMax]);
	title(sprintf('%d',i));
	drawnow
end
close;

%extract output
dblStimDur = sStimParams.dblStimDur;
dblTrialDur = sStimInputs.dblTrialDur;
varDeltaSyn = sStimInputs.varDeltaSyn;
dblVisSpacing = sStimInputs.dblVisSpacing;
matBlankLGN_ON = sStimInputs.matBlankLGN_ON;
matBlankLGN_OFF = sStimInputs.matBlankLGN_OFF;
cellLGN_ON = sStimInputs.cellLGN_ON;
cellLGN_OFF = sStimInputs.cellLGN_OFF;
vecTrialOris = sStimInputs.vecTrialOris;
vecTrialOriIdx = sStimInputs.vecTrialOriIdx;
vecTrialSFs = sStimInputs.vecTrialSFs;
vecTrialSFIdx = sStimInputs.vecTrialSFIdx;
vecTrialStimType = sStimInputs.vecTrialStimType;
vecStimStartSecs = sStimInputs.vecStimStartSecs;
vecStimStopSecs = sStimInputs.vecStimStartSecs+sStimParams.dblStimDur;
vecTrialStartSecs = vecStimStartSecs-sStimParams.dblPreStimBlankDur;
vecTrialEndSecs = sStimInputs.vecTrialEndSecs;
vecStimTypeSFs = sStimInputs.vecStimTypeSFs;
vecStimTypeOris = sStimInputs.vecStimTypeOris;

dblSimDur = max(vecTrialEndSecs); %seconds
fprintf(' .. Parameters: time step %.1f ms for %.3f seconds\n',dblDeltaT*1000,dblSimDur);

%% get LGN to V1 parameters
%get spatial receptive fields
intImX = size(varDeltaSyn,1);
intImY = size(varDeltaSyn,2);

vecSpaceX = dblVisSpacing*((-(intImX - 1)/2):(intImX - 1)/2);
vecSpaceY = dblVisSpacing*((-(intImY - 1)/2):(intImY - 1)/2);
[matMeshX,matMeshY] = meshgrid(vecSpaceX,vecSpaceY);

%create preferred stimulus features for all cells in V1
intColumns = 48; %48 / 252
vecDefinitionV1PrefOri = 0:pi/intColumns:(pi-pi/intColumns);
vecDefinitionV1SpatFreq = 2.^[-4:1];
vecDefinitionV1CellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
intCellsV1 = numel(vecDefinitionV1PrefOri) * numel(vecDefinitionV1SpatFreq) * numel(vecDefinitionV1CellTypes);

%distribute RF center locations
intCellsPerSpatialDim = ceil(sqrt(intCellsV1));
vecPrefDistroX = linspace(min(vecSpaceX)/2,max(vecSpaceX)/2,intCellsPerSpatialDim);
vecPrefDistroY = linspace(min(vecSpaceY)/2,max(vecSpaceY)/2,intCellsPerSpatialDim);
[matPrefRF_X,matPrefRF_Y] = meshgrid(vecPrefDistroX,vecPrefDistroY);
vecRandRF = randperm(intCellsPerSpatialDim.^2,intCellsV1);
dblSpatialStep = min(diff(vecPrefDistroX));
vecPrefRF_X_V1 = matPrefRF_X(vecRandRF)+(dblSpatialStep)*(rand(1,intCellsV1)-0.5);
vecPrefRF_Y_V1 = matPrefRF_Y(vecRandRF)+(dblSpatialStep)*(rand(1,intCellsV1)-0.5);
matPrefGabors = nan(numel(vecSpaceX),numel(vecSpaceY),intCellsV1);


%create cell-based parameters
vecPrefOriV1 = nan(1,intCellsV1);
vecPrefSFV1 = nan(1,intCellsV1);
vecCellTypesV1 = nan(1,intCellsV1);
intCounter=0;
intIncrement = numel(vecDefinitionV1CellTypes);
for dblOri=vecDefinitionV1PrefOri
	for dblSF = vecDefinitionV1SpatFreq
		vecPrefOriV1((intCounter+1):(intCounter+intIncrement)) = dblOri;
		vecPrefSFV1((intCounter+1):(intCounter+intIncrement)) = dblSF;
		vecCellTypesV1((intCounter+1):(intCounter+intIncrement)) = vecDefinitionV1CellTypes;
		intCounter = intCounter + intIncrement;
	end
end

%connection definition LGN
vecCellTypeLabels = unique(vecCellTypesV1);
vecConnsPerTypeON = [24 16]; %[pyramid interneuron]
vecConnsPerTypeOFF = [24 16]; %[pyramid interneuron]

dblSigmaX = 0.7; %length of gabor response
dblSigmaY = 0.47; %width of gabor response
intCellsV1 = length(vecCellTypesV1);
vecCellAreaV1 = ones(1,intCellsV1); %1, V1; 2, V2
vecPrefPsiV1 = 2*pi*rand(1,intCellsV1);%phase offset
intConnsLGN_to_CortOFF = sum(vecConnsPerTypeOFF(vecDefinitionV1CellTypes))*(numel(vecDefinitionV1PrefOri) * numel(vecDefinitionV1SpatFreq));
intConnsLGN_to_CortON = sum(vecConnsPerTypeON(vecDefinitionV1CellTypes))*(numel(vecDefinitionV1PrefOri) * numel(vecDefinitionV1SpatFreq));

%create connection parameters from LGN to cortex
vecConductance_FromLGN_ToCort = [5.5 6]*0.3; %to [pyramid interneuron]
vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

matSynConnON_to_Cort = nan(intConnsLGN_to_CortON,2);
vecSynWeightON_to_Cort = nan(intConnsLGN_to_CortON,1);
vecSynConductanceON_to_Cort = nan(intConnsLGN_to_CortON,1);
vecSynDelayON_to_Cort = nan(intConnsLGN_to_CortON,1);

matSynConnOFF_to_Cort = nan(intConnsLGN_to_CortOFF,2);
vecSynWeightOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);
vecSynConductanceOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);
vecSynDelayOFF_to_Cort = nan(intConnsLGN_to_CortOFF,1);

intC_ON = 1;
intC_OFF = 1;
for intNeuron=1:intCellsV1
	%get cell parameters
	dblPsi = vecPrefPsiV1(intNeuron); %cell's phase offset
	dblTheta=vecPrefOriV1(intNeuron); %cell's preferred orientation
	dblPrefSF = vecPrefSFV1(intNeuron);
	dblPrefRF_X = vecPrefRF_X_V1(intNeuron);
	dblPrefRF_Y = vecPrefRF_Y_V1(intNeuron);
	
	intCellType = vecCellTypesV1(intNeuron); %pyramidal cell or interneuron
	intConnsON = vecConnsPerTypeON(intCellType);
	intConnsOFF = vecConnsPerTypeOFF(intCellType);
	
	%get assignment vectors
	vecConnsAssignON = intC_ON:(intC_ON+intConnsON-1);
	vecConnsAssignOFF = intC_OFF:(intC_OFF+intConnsOFF-1);
	intC_ON = intC_ON + intConnsON;
	intC_OFF = intC_OFF + intConnsOFF;
	
	%% create Gabor
	% get rotation matrices
	matX_theta=(matMeshX*cos(dblTheta)+matMeshY*sin(dblTheta))+dblPrefRF_X;
	matY_theta=(-matMeshX*sin(dblTheta)+matMeshY*cos(dblTheta))+dblPrefRF_Y;
	
	%get gabor
	matG = exp(-.5*(matX_theta.^2/dblSigmaX^2+matY_theta.^2/dblSigmaY^2)).*cos(2*pi*dblPrefSF*matX_theta+dblPsi);
	matPrefGabors(:,:,intNeuron) = matG;
	
	%get connection probabilities
	matProbConn_ON = max(0,matG) / sum(sum(max(0,matG)));
	matProbConn_OFF = max(0,-matG) / sum(sum(max(0,-matG)));
	
	%{
	%% plot
	subplot(2,2,1)
	imagesc(matG);
	subplot(2,2,3)
	imagesc(matProbConn_ON);
	subplot(2,2,4)
	imagesc(matProbConn_OFF);
	%}
	
	%% ON fields
	%ON; linearize, get random connections, and reshape to matrix
	vecProbLinON = matProbConn_ON(:);
	vecProbCumSumON = cumsum(vecProbLinON);
	%matConnON = false(size(matProbConn_ON));
	vecConnsON = [];
	while numel(vecConnsON) < intConnsON
		vecConnsON = unique([vecConnsON sum(bsxfun(@gt,rand(1,intConnsON),vecProbCumSumON),1)+1]);
		if numel(vecConnsON) > intConnsON
			vecConnsON = vecConnsON(randperm(numel(vecConnsON),intConnsON));
		end
	end
	%[I,J] = ind2sub([intImX intImY],vecConnsON);
	%matConnON = logical(getFillGrid(zeros(intImX,intImY),I,J,true));
	
	%put in 2D connection matrix
	matConnections2D = [vecConnsON' repmat(intNeuron,[length(vecConnsON) 1])];
	matSynConnON_to_Cort(vecConnsAssignON,:) = matConnections2D;
	
	%put synaptic properties in vector
	vecSynWeightON_to_Cort(vecConnsAssignON) = vecProbLinON(vecConnsON);
	vecSynConductanceON_to_Cort(vecConnsAssignON) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignON));
	vecSynDelayON_to_Cort(vecConnsAssignON) = max(0,normrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignON)));
	
	%% OFF fields
	%OFF; linearize, get random connections
	vecProbLinOFF = matProbConn_OFF(:);
	vecProbCumSumOFF = cumsum(vecProbLinOFF);
	%matConnOFF = false(size(matProbConn_OFF));
	vecConnsOFF = [];
	while numel(vecConnsOFF) < intConnsOFF
		vecConnsOFF = unique([vecConnsOFF sum(bsxfun(@gt,rand(1,intConnsOFF),vecProbCumSumOFF),1)+1]);
		if numel(vecConnsOFF) > intConnsOFF
			vecConnsOFF = vecConnsOFF(randperm(numel(vecConnsOFF),intConnsOFF));
		end
	end
	%[I,J] = ind2sub([intImX intImY],vecConnsOFF);
	%matConnOFF = getFillGrid(zeros(intImX,intImY),I,J,true);
	
	%put in 2D connection matrix
	matConnections2D = [vecConnsOFF' repmat(intNeuron,[length(vecConnsOFF) 1])];
	matSynConnOFF_to_Cort(vecConnsAssignOFF,:) = matConnections2D;
	
	%put synaptic properties in vector
	vecSynWeightOFF_to_Cort(vecConnsAssignOFF) = vecProbLinOFF(vecConnsOFF);
	vecSynConductanceOFF_to_Cort(vecConnsAssignOFF) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignOFF));
	vecSynDelayOFF_to_Cort(vecConnsAssignOFF) = max(0,normrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignOFF)));
	
end

fprintf(' .. Created LGN to V1 connectivity; %d active OFF field synapses, and %d active ON field synapses [%s]\n',intConnsLGN_to_CortOFF,intConnsLGN_to_CortON,getTime);

%% define LIF model parameters
%% neuron parameters
%pyramids
sParamCells{1}.dblThresh = -55; %mV
sParamCells{1}.dblTau_Peak = 1/1000; %ms
sParamCells{1}.dblV_E = 0; %excitatory reversal potential
sParamCells{1}.dblV_I = -70; %inhibitory reversal potential
sParamCells{1}.dblV_AHP = -90;%after-hyperpol reversal potential
sParamCells{1}.dblV_Leak = -65; %resting potential
sParamCells{1}.dblCm = 0.5; %capacitance
sParamCells{1}.dblG_Leak = 25; %nS
sParamCells{1}.dblG_AHP = 40; %nS

%interneurons
sParamCells{2}.dblThresh = -55; %mV
sParamCells{2}.dblTau_Peak = 2/1000; %ms
sParamCells{2}.dblV_E = 0; %excitatory reversal potential
sParamCells{2}.dblV_I = -70; %inhibitory reversal potential
sParamCells{2}.dblV_AHP = -90;%after-hyperpol reversal potential
sParamCells{2}.dblV_Leak = -65; %resting potential
sParamCells{2}.dblCm = 0.2; %capacitance
sParamCells{2}.dblG_Leak = 20; %nS
sParamCells{2}.dblG_AHP = 20; %nS

%define tau by cell type
vecTauPeakByType = [sParamCells{1}.dblTau_Peak sParamCells{2}.dblTau_Peak];

%% connection parameters
%number of connections
dblScalingFactor = 4; %4
matConnCortFromTo(1,:) = [40 40]*dblScalingFactor; %from pyramid to [pyr inter]
matConnCortFromTo(2,:) = [30 30]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
matConductancesFromTo(1,:) = [1.1 1.6]/dblScalingFactor; %from pyramid to [pyr inter]
matConductancesFromTo(2,:) = [1.5 1.0]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
dblDelayMeanCortToCort = 3/1000; %in ms
dblDelaySDCortToCort = 1/1000; %in ms

%msg
fprintf(' .. Created LIF model parameters [%s]\n',getTime);

%% create cortical connectivity and parameter vectors
%connection probability excitatory cells
vecConnProbSD = ang2rad([7.5 80]); %degrees difference in pref ori, for [pyramid interneuron]

%pre-allocate cortical synaptic connection matrix
vecTotalConnectionsPerType = sum(matConnCortFromTo,1);
intTotalV1Connections = sum(vecTotalConnectionsPerType(vecCellTypesV1));
matSynV1FromTo = nan(intTotalV1Connections,2); %[From x To]
vecSynV1Conductance = nan(intTotalV1Connections,1); %synaptic conductance
vecSynV1Weight = ones(intTotalV1Connections,1); %synaptic weight
vecSynV1Delay = nan(intTotalV1Connections,1); %synaptic delay
vecSynV1ExcInh = nan(intTotalV1Connections,1); %synaptic sign; inhibitory or excitatory
vecSynV1Type = nan(intTotalV1Connections,1); %synaptic type; intra- or inter-areal

%pre-allocate cell parameters
vecCellThresh = nan(intCellsV1,1);
vecCellTauPeak = nan(intCellsV1,1);
vecCellV_E = nan(intCellsV1,1);
vecCellV_I = nan(intCellsV1,1);
vecCellV_AHP = nan(intCellsV1,1);
vecCellV_Leak = nan(intCellsV1,1);
vecCellCm =	nan(intCellsV1,1);
vecCellG_Leak = nan(intCellsV1,1);
vecCellG_AHP = nan(intCellsV1,1);

%get random connections
intC = 1;
for intNeuron=1:intCellsV1
	dblPrefOri = vecPrefOriV1(intNeuron);
	vecOriDiff = abs(circ_dist(vecPrefOriV1*2,dblPrefOri*2)/2)';
	
	%get this cell's type
	intTargetCellType = vecCellTypesV1(intNeuron); %pyramidal cell or interneuron
	vecConnsFromTypesToThisType = matConnCortFromTo(:,intTargetCellType); %from [pyr inter]
	
	%assign parameter values
	vecCellThresh(intNeuron) = sParamCells{intTargetCellType}.dblThresh;
	vecCellTauPeak(intNeuron) = sParamCells{intTargetCellType}.dblTau_Peak; %ms
	vecCellV_E(intNeuron) = sParamCells{intTargetCellType}.dblV_E; %excitatory reversal potential
	vecCellV_I(intNeuron) = sParamCells{intTargetCellType}.dblV_I; %inhibitory reversal potential
	vecCellV_AHP(intNeuron) = sParamCells{intTargetCellType}.dblV_AHP;%after-hyperpol reversal potential
	vecCellV_Leak(intNeuron) = sParamCells{intTargetCellType}.dblV_Leak; %resting potential
	vecCellCm(intNeuron) = sParamCells{intTargetCellType}.dblCm; %capacitance
	vecCellG_Leak(intNeuron) = sParamCells{intTargetCellType}.dblG_Leak; %nS
	vecCellG_AHP(intNeuron) = sParamCells{intTargetCellType}.dblG_AHP; %nS
	
	%loop through subtypes to make connections
	for intSubIdx=1:length(vecCellTypeLabels)
		intSourceCellType = vecCellTypeLabels(intSubIdx);
		indThisType = vecCellTypesV1==intSourceCellType;
		intConnsFromThisType = vecConnsFromTypesToThisType(intSourceCellType);
		dblConnProbOriSD = vecConnProbSD(intSourceCellType);
		vecConnProb = normpdf(vecOriDiff,0,dblConnProbOriSD);
		vecConnProb(~indThisType) = 0; %remove other types
		vecConnProb(intNeuron) = 0; %remove the neuron itself
		vecConnProb = vecConnProb/sum(vecConnProb); %normalize so integral is 1
		vecProbCumSum = cumsum(vecConnProb);
		
		%get assignment vectors
		vecConnsAssign = intC:(intC+intConnsFromThisType-1);
		intC = intC + intConnsFromThisType;
		
		%make connections
		vecConns = [];
		while numel(vecConns) < intConnsFromThisType
			vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsFromThisType),vecProbCumSum),1)+1]);
			if numel(vecConns) > intConnsFromThisType
				vecConns = vecConns(randperm(numel(vecConns),intConnsFromThisType));
			end
		end
		
		%assign to connection matrix
		dblG = matConductancesFromTo(intSourceCellType,intTargetCellType);
		matSynV1FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
		vecSynV1Conductance(vecConnsAssign) = dblG;
		vecSynV1Delay(vecConnsAssign) = max(0,normrnd(dblDelayMeanCortToCort,dblDelaySDCortToCort,size(vecConnsAssign)));
		vecSynV1ExcInh(vecConnsAssign) = intSourceCellType;
		vecSynV1Type(vecConnsAssign) = 1; %1; intra-areal, 2 inter-areal
	end
end
fprintf(' .. Created V1 connectivity; %d active corticocortical synapses between %d cells [%s]\n',intTotalV1Connections,intCellsV1,getTime);

%% build V1-V2 connectivity
%set parameters for V2
dblPercProjectingV1_to_V2 = 10; %El-Shamayleh et al (Movshon), 2013;
dblSpatialDropoffV1V2 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
dblSpatialDropoffInterneuronsV2 = 3; %for interneurons
intCellsV2 = round(intCellsV1/2);

%create cell-based parameters
vecDefinitionV2CellTypes = [1 1 1 1 2];
vecCellTypesV2 = repmat(vecDefinitionV2CellTypes,[1 ceil(intCellsV2/numel(vecDefinitionV2CellTypes))]);
vecCellTypesV2(intCellsV2+1:end) = [];
	
%LGN projects to V2 (Bullier & Kennedy, 1983), but how much?

%Schmid et al. (Smirnakis), 2009; 70-95% reduction in V2 activity after lesion of V1

%pulvinar plays major role in projecting to V2, and received most of its input from V1 (Sincich & Horton, 2005)
%V2<->V1; 14,000 feedforward and 11,000 feedback connections per mm^2; if FF connections are 100%, FB are 80%  (Sincich & Horton, 2005; Rockland, 1997)
%connections from V1->V2 are 20-25 fold larger than LGN->V1 (van Essen, 2005)

%V2 cells respond to angles and features; simple model: combining two V1
%neurons with similar receptive field locations (Ito and Komatsu, 2004)


%% calculate similarity between receptive fields of V1 cells
fprintf(' .. Building V2 exemplar fields for %d cells... [%s]\n',intCellsV2,getTime);
matTuningSimilarity=nan(intCellsV1,intCellsV1);
parfor i1=1:intCellsV1
	for i2=1:intCellsV1
		matTuningSimilarity(i1,i2) = corr2(matPrefGabors(:,:,i1),matPrefGabors(:,:,i2)); %#ok<PFBNS>
	end
end
matTuningSimilarity = abs(tril(matTuningSimilarity,-1));
matTuningSimilarity(vecCellTypesV1==2,:) = 0; %remove interneurons
matTuningSimilarity(:,vecCellTypesV1==2) = 0; %remove interneurons
matTuningSimilarity(matTuningSimilarity==0)=nan;
matTotalSim = matTuningSimilarity.^2;% + matTuningSimilarity;

%% create RF exemplars for V2 cells
%matTuningSimilarity = matProbConnPairV1_to_V2;
%select r
vecProbComb=(abs(matTotalSim(:))./nansum(abs(matTotalSim(:))));
[vecProbSorted,vecSort]=sort(vecProbComb);
vecSelectProb = cumsum(vecProbSorted)/nansum(vecProbSorted);
matExemplarFieldsV2 = nan(size(matPrefGabors,1),size(matPrefGabors,2),intCellsV2);
vecPrefRF_X_V2 = nan(1,intCellsV2);
vecPrefRF_Y_V2 = nan(1,intCellsV2);

for intFieldV2=find(vecCellTypesV2==1)
	intComb = vecSort(find(vecSelectProb>rand(1),1,'first'));
	[I,J] = ind2sub([intCellsV1 intCellsV1],intComb);
	matGaborExemplar = matPrefGabors(:,:,I)+matPrefGabors(:,:,J);
	matExemplarFieldsV2(:,:,intFieldV2) = matGaborExemplar;
	vecCoM = calcCenterOfMass(abs(matGaborExemplar));
	vecPrefRF_X_V2(intFieldV2) = vecCoM(2);
	vecPrefRF_Y_V2(intFieldV2) = vecCoM(1);

	
	%imagesc(matGaborExemplar)
	%title(sprintf('Dist: %.3f; sim: %.3f; tot: %.6f',matProbConnPairV1_to_V2(I,J),matTuningSimilarity(I,J),matTotalSim(I,J)))
	%drawnow;
end
vecPrefRF_X_V2 = dblVisSpacing*(vecPrefRF_X_V2 - intImX/2);
vecPrefRF_Y_V2 = dblVisSpacing*(vecPrefRF_Y_V2 - intImY/2);

%% build similarity of V1 fields to V2 exemplars
fprintf(' .. Calculating similarity of V1 Gabors to V2 exemplars... [%s]\n',getTime);
matSimilarityGaborExamplar = nan(intCellsV1,intCellsV2);
parfor intCellV1=1:intCellsV1
	for intCellV2=1:intCellsV2
		matSimilarityGaborExamplar(intCellV1,intCellV2) = corr2(matPrefGabors(:,:,intCellV1),matExemplarFieldsV2(:,:,intCellV2)) %#ok<PFBNS>
	end
end
%msg
fprintf(' .. Created V2 fields, building V1 => V2 connectivity... [%s]\n',getTime);


%% build preferred RF location for interneurons
%distribute RF center locations
vecCellFractionsV2 = [0.8 0.2]; %[pyramid interneuron]
intInterneuronsV2 = round(intCellsV2*vecCellFractionsV2(2));
intInterneuronsPerSpatialDimV2 = ceil(sqrt(intInterneuronsV2));
vecPrefDistroX = linspace(min(vecSpaceX)/2,max(vecSpaceX)/2,intInterneuronsPerSpatialDimV2);
vecPrefDistroY = linspace(min(vecSpaceY)/2,max(vecSpaceY)/2,intInterneuronsPerSpatialDimV2);
[matPrefRF_X,matPrefRF_Y] = meshgrid(vecPrefDistroX,vecPrefDistroY);
vecRandRF = randperm(intInterneuronsPerSpatialDimV2.^2,intInterneuronsV2);
dblSpatialStep = min(diff(vecPrefDistroX));
vecPrefRF_X_IntV2 = matPrefRF_X(vecRandRF)+(dblSpatialStep)*(rand(1,intInterneuronsV2)-0.5);
vecPrefRF_Y_IntV2 = matPrefRF_Y(vecRandRF)+(dblSpatialStep)*(rand(1,intInterneuronsV2)-0.5);
%assign to V2 vector
vecPrefRF_X_V2(vecCellTypesV2==2) = vecPrefRF_X_IntV2;
vecPrefRF_Y_V2(vecCellTypesV2==2) = vecPrefRF_Y_IntV2;


%% create V1=>V2 connections
%still to do: make V1->V2 projections, V2 recurrence, and feedback
%connection parameters
vecConnsPerTypeV1V2 = [48 32];%[48 32]; %pyr/int
dblInterArealFactor = 2;
matConductancesFromToV1V2(1,:) = [1.1 1.6]*dblInterArealFactor; %from pyramid to [pyr inter]
matConductancesFromToV1V2(2,:) = [1.5 1.0]*dblInterArealFactor; %from inter to [pyr inter]

%synaptic delays
dblDelayMeanV1ToV2 = 4/1000; %in ms
dblDelaySDV1ToV2 = 1/1000; %in ms

%pre-allocate synaptic variables
intTotalV1V2Connections = sum(vecConnsPerTypeV1V2(vecCellTypesV2));
matSynV1V2FromTo = nan(intTotalV1V2Connections,2); %[From x To]
vecSynV1V2Conductance = nan(intTotalV1V2Connections,1); %synaptic conductance
vecSynV1V2Weight = nan(intTotalV1V2Connections,1); %synaptic strength
vecSynV1V2Delay = nan(intTotalV1V2Connections,1); %synaptic delay
vecSynV1V2ExcInh = nan(intTotalV1V2Connections,1); %synaptic sign; inhibitory or excitatory
vecSynV1V2Type = nan(intTotalV1V2Connections,1); %synaptic type; intra- or inter-areal
matFieldsV2 = nan(size(matExemplarFieldsV2));

%pre-allocate cell parameters
vecCellV2Thresh = nan(intCellsV2,1);
vecCellV2TauPeak = nan(intCellsV2,1);
vecCellV2V_E = nan(intCellsV2,1);
vecCellV2V_I = nan(intCellsV2,1);
vecCellV2V_AHP = nan(intCellsV2,1);
vecCellV2V_Leak = nan(intCellsV2,1);
vecCellV2Cm =	nan(intCellsV2,1);
vecCellV2G_Leak = nan(intCellsV2,1);
vecCellV2G_AHP = nan(intCellsV2,1);

%get random connections
intC = 1;
for intCellV2=1:intCellsV2
	%get cell parameters
	dblPrefRF_X = vecPrefRF_X_V2(intCellV2);
	dblPrefRF_Y = vecPrefRF_Y_V2(intCellV2);
	intNeuron = intCellV2+intCellsV1;
	
	%type
	intTargetCellType = vecCellTypesV2(intCellV2); %pyramidal cell or interneuron
	cellStrType = {'pyramidal ','inter'};

	%assign parameter values
	vecCellV2Thresh(intCellV2) = sParamCells{intTargetCellType}.dblThresh;
	vecCellV2TauPeak(intCellV2) = sParamCells{intTargetCellType}.dblTau_Peak; %ms
	vecCellV2V_E(intCellV2) = sParamCells{intTargetCellType}.dblV_E; %excitatory reversal potential
	vecCellV2V_I(intCellV2) = sParamCells{intTargetCellType}.dblV_I; %inhibitory reversal potential
	vecCellV2V_AHP(intCellV2) = sParamCells{intTargetCellType}.dblV_AHP;%after-hyperpol reversal potential
	vecCellV2V_Leak(intCellV2) = sParamCells{intTargetCellType}.dblV_Leak; %resting potential
	vecCellV2Cm(intCellV2) = sParamCells{intTargetCellType}.dblCm; %capacitance
	vecCellV2G_Leak(intCellV2) = sParamCells{intTargetCellType}.dblG_Leak; %nS
	vecCellV2G_AHP(intCellV2) = sParamCells{intTargetCellType}.dblG_AHP; %nS
	
	%msg
	%fprintf('Building V1=>V2 connections for V2 %sneuron %d/%d [%s]\n',cellStrType{intTargetCellType},intCellV2,intCellsV2,getTime);
	
	%if interneuron, RF center distance based, otherwise exemplar-field
	%similarity based
	if intTargetCellType==1
		%exemplar-based
		vecTuningSimilarity = matSimilarityGaborExamplar(:,intCellV2);
		vecConnProb = vecTuningSimilarity;
		vecConnProb(vecConnProb<0)=0;
	elseif intTargetCellType==2
		%RF center based
		vecDistanceRF = sqrt((vecPrefRF_Y_V1-dblPrefRF_Y).^2 + (vecPrefRF_X_V1-dblPrefRF_X).^2); %euclidian distance
		vecConnProb = normpdf(vecDistanceRF',0,dblSpatialDropoffV1V2);
	end
	
	%select projection cells
	intSourceCellType = 1; %only pyramids project
	indProjectingCells = vecCellTypesV1==intSourceCellType;
	vecConnProb(~indProjectingCells) = 0; %remove other types
	vecConnProb = vecConnProb/sum(vecConnProb); %normalize so integral is 1
	vecProbCumSum = cumsum(vecConnProb);
	
	%get assignment vectors
	intConnsToThisType = vecConnsPerTypeV1V2(intTargetCellType);
	vecConnsAssign = intC:(intC+intConnsToThisType-1);
	intC = intC + intConnsToThisType;
	
	%make connections
	vecConns = [];
	while numel(vecConns) < intConnsToThisType
		vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsToThisType),vecProbCumSum),1)+1]);
		if numel(vecConns) > intConnsToThisType
			vecConns = vecConns(randperm(numel(vecConns),intConnsToThisType));
		end
	end
	vecTheseWeights = vecConnProb(vecConns);
	vecTheseWeights = (vecTheseWeights./mean(vecTheseWeights));
	
	%assign to connection matrix
	dblG = matConductancesFromToV1V2(intSourceCellType,intTargetCellType);
	matSynV1V2FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
	vecSynV1V2Conductance(vecConnsAssign) = dblG; %synaptic conductance
	vecSynV1V2Weight(vecConnsAssign) = vecTheseWeights; %synaptic conductance
	vecSynV1V2Delay(vecConnsAssign) = max(0,normrnd(dblDelayMeanV1ToV2,dblDelaySDV1ToV2,size(vecConnsAssign))); %synaptic delay
	vecSynV1V2ExcInh(vecConnsAssign) = intSourceCellType; %synaptic sign; inhibitory or excitatory
	vecSynV1V2Type(vecConnsAssign) = 2; %1; intra-areal, 2 inter-areal
	
	%build resultant V2 fields
	matFieldsV2(:,:,intCellV2) = xmean(bsxfun(@mtimes,matPrefGabors(:,:,vecConns),reshape(vecTheseWeights,[1 1 intConnsToThisType])),3);
	%{
	subplot(2,2,1);
	imagesc(matExemplarFieldsV2(:,:,intCellV2));
	subplot(2,2,2);
	imagesc(matFieldsV2(:,:,intCellV2));
	subplot(2,2,3);
	imagesc(xmean(matPrefGabors(:,:,vecConns),3));
	drawnow;pause;
	%}
end
fprintf(' .. Created V1 => V2 connectivity; %d active inter-areal synapses from %d V1 to %d V2 cells [%s]\n',intTotalV1V2Connections,intCellsV1,intCellsV2,getTime);

%% build RF similarities V2
vecPrefRF_X_V2 = nan(1,intCellsV2);
vecPrefRF_Y_V2 = nan(1,intCellsV2);
matSimilarityFieldsV2 = zeros(intCellsV2,intCellsV2);
parfor intCellV2_1=1:intCellsV2
	
	vecCoM = calcCenterOfMass(abs(matFieldsV2(:,:,intCellV2_1)));
	vecPrefRF_X_V2(intCellV2_1) = vecCoM(2);
	vecPrefRF_Y_V2(intCellV2_1) = vecCoM(1);

	for intCellV2_2=1:intCellsV2
		matSimilarityFieldsV2(intCellV2_1,intCellV2_2) = corr2(matFieldsV2(:,:,intCellV2_1),matFieldsV2(:,:,intCellV2_2)); %#ok<PFBNS>
	end
end
vecPrefRF_X_V2 = dblVisSpacing*(vecPrefRF_X_V2 - intImX/2);
vecPrefRF_Y_V2 = dblVisSpacing*(vecPrefRF_Y_V2 - intImY/2);
matSimilarityFieldsV2(matSimilarityFieldsV2<0)=0;
matSimilarityFieldsV2(diag(diag(true(size(matSimilarityFieldsV2))))) = 0;

%msg
fprintf(' .. Created V2 field similarities, building V2 => V2 connectivity... [%s]\n',getTime);


%% V2 => V2
%pre-allocate cortical synaptic connection matrix
intTotalV2Connections = sum(vecTotalConnectionsPerType(vecCellTypesV2));
matSynV2FromTo = nan(intTotalV2Connections,2); %[From x To]
vecSynV2Conductance = nan(intTotalV2Connections,1); %synaptic conductance
vecSynV2Weight = ones(intTotalV2Connections,1); %synaptic weight
vecSynV2Delay = nan(intTotalV2Connections,1); %synaptic delay
vecSynV2ExcInh = nan(intTotalV2Connections,1); %synaptic sign; inhibitory or excitatory
vecSynV2Type = nan(intTotalV2Connections,1); %synaptic type; intra- or inter-areal

%get random connections
intC = 1;
for intCellV2=1:intCellsV2
	%get this cell's type
	intTargetCellType = vecCellTypesV2(intCellV2); %pyramidal cell or interneuron
	vecConnsFromTypesToThisType = matConnCortFromTo(:,intTargetCellType); %from [pyr inter]
	
	%get cell parameters
	dblPrefRF_X = vecPrefRF_X_V2(intCellV2);
	dblPrefRF_Y = vecPrefRF_Y_V2(intCellV2);
	intNeuron = intCellV2+intCellsV1;
	
	%loop through subtypes to make connections
	for intSubIdx=1:length(vecCellTypeLabels)
		intSourceCellType = vecCellTypeLabels(intSubIdx);
		indThisType = vecCellTypesV2==intSourceCellType;
		intConnsFromThisType = vecConnsFromTypesToThisType(intSourceCellType);
		
		%if interneuron, RF center distance based, otherwise exemplar-field
		%similarity based
		if intSourceCellType==1
			%exemplar-based
			vecTuningSimilarity = matSimilarityFieldsV2(:,intCellV2);
			vecConnProb = vecTuningSimilarity;
			vecConnProb(vecConnProb<0)=0;
		elseif intSourceCellType==2
			%RF center based
			vecDistanceRF = sqrt((vecPrefRF_Y_V2-dblPrefRF_Y).^2 + (vecPrefRF_X_V2-dblPrefRF_X).^2); %euclidian distance
			vecConnProb = normpdf(vecDistanceRF',0,dblSpatialDropoffInterneuronsV2);
		end
		
		vecConnProb(~indThisType) = 0; %remove other types
		vecConnProb(intNeuron) = 0; %remove the neuron itself
		vecConnProb = vecConnProb/sum(vecConnProb); %normalize so integral is 1
		vecProbCumSum = cumsum(vecConnProb);
		
		%get assignment vectors
		vecConnsAssign = intC:(intC+intConnsFromThisType-1);
		intC = intC + intConnsFromThisType;
		
		%make connections
		vecConns = [];
		while numel(vecConns) < intConnsFromThisType
			vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsFromThisType),vecProbCumSum),1)+1]);
			if numel(vecConns) > intConnsFromThisType
				vecConns = vecConns(randperm(numel(vecConns),intConnsFromThisType));
			end
		end
		vecConns = vecConns + intCellsV1;
		
		%assign to connection matrix
		dblG = matConductancesFromTo(intSourceCellType,intTargetCellType);
		matSynV2FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
		vecSynV2Conductance(vecConnsAssign) = dblG;
		vecSynV2Weight(vecConnsAssign) = 1;
		vecSynV2Delay(vecConnsAssign) = max(0,normrnd(dblDelayMeanCortToCort,dblDelaySDCortToCort,size(vecConnsAssign)));
		vecSynV2ExcInh(vecConnsAssign) = intSourceCellType;
		vecSynV2Type(vecConnsAssign) = 1; %1; intra-areal, 2 inter-areal
	end
end
fprintf(' .. Created V2 connectivity; %d active corticocortical synapses between %d cells [%s]\n',intTotalV2Connections,intCellsV2,getTime);

%% aggregate variables
%tuning
intCortexCells = intCellsV1+intCellsV2;
vecCellTypes = cat(2,vecCellTypesV1,vecCellTypesV2); %1, pyr, 2, int
vecCellArea = cat(2,vecCellAreaV1,2*ones(1,intCellsV2)); %1, V1; 2, V2
vecPrefPsi = cat(2,vecPrefPsiV1,nan(1,intCellsV2));
vecPrefOri = cat(2,vecPrefOriV1,nan(1,intCellsV2));
vecPrefSF = cat(2,vecPrefSFV1,nan(1,intCellsV2));
vecPrefRF_X = cat(2,vecPrefRF_X_V1,vecPrefRF_X_V2);
vecPrefRF_Y = cat(2,vecPrefRF_Y_V1,vecPrefRF_Y_V2);

%synapses
matSynFromTo = cat(1,matSynV1FromTo,matSynV1V2FromTo,matSynV2FromTo);
vecSynConductance = cat(1,vecSynV1Conductance,vecSynV1V2Conductance,vecSynV2Conductance);
vecSynWeight = cat(1,vecSynV1Weight,vecSynV1V2Weight,vecSynV2Weight);
vecSynDelay = cat(1,vecSynV1Delay,vecSynV1V2Delay,vecSynV2Delay);
vecSynExcInh = cat(1,vecSynV1ExcInh,vecSynV1V2ExcInh,vecSynV2ExcInh);
vecSynType = cat(1,vecSynV1Type,vecSynV1V2Type,vecSynV2Type);

%cell parameters
vecCellThresh = cat(1,vecCellThresh,vecCellV2Thresh);
vecCellTauPeak = cat(1,vecCellTauPeak,vecCellV2TauPeak);
vecCellV_E =cat(1,vecCellV_E,vecCellV2V_E);
vecCellV_I =cat(1,vecCellV_I,vecCellV2V_I);
vecCellV_AHP = cat(1,vecCellV_AHP,vecCellV2V_AHP);
vecCellV_Leak = cat(1,vecCellV_Leak,vecCellV2V_Leak);
vecCellCm =	cat(1,vecCellCm,vecCellV2Cm);
vecCellG_Leak =cat(1,vecCellG_Leak,vecCellV2G_Leak);
vecCellG_AHP = cat(1,vecCellG_AHP,vecCellV2G_AHP);

%% prepare variables for simulation
% run stitched together trials
vecThisV = ones(intCortexCells,1)*-65;
vecPrevV = vecThisV;
boolStimPresent = false;
intPrevTrial = 0;
intTrialT = 0;
vecOverallT = dblDeltaT:dblDeltaT:dblSimDur;
fprintf('\n   >>> Starting simulation run [%s]\n',getTime);
intIter = 0;

intLGN_CellsPerStream = size(cellLGN_ON{1},1)*size(cellLGN_ON{1},2);
cellSpikeTimesLGN_ON = cell(intLGN_CellsPerStream,1);
cellSpikeTimesLGN_OFF = cell(intLGN_CellsPerStream,1);
cellSpikeTimesCortex = cell(intCortexCells,1);

vecSpikeCounterLGN_ON = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterLGN_OFF = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterCortex = zeros(intCortexCells,1);

vecSpikeCounterPreAllocatedLGN_ON = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterPreAllocatedLGN_OFF = zeros(intLGN_CellsPerStream,1);
vecSpikeCounterPreAllocatedCortex = zeros(intCortexCells,1);

intPreAllocationSize = 1000;

%% assign to structure
sData.intCellsV2 = intCellsV2;
sData.intCellsV1 = intCellsV1;
sData.intCortexCells = intCortexCells;
sData.vecOverallT = vecOverallT;
sData.dblDeltaT = dblDeltaT;

sData.matSynFromTo = matSynFromTo;
sData.dblSynSpikeMem = dblSynSpikeMem;
sData.vecSynExcInh = vecSynExcInh;
sData.vecSynDelay = vecSynDelay;
sData.vecSynWeight = vecSynWeight;
sData.vecSynConductance = vecSynConductance;
sData.vecSynType = vecSynType;

sData.vecCellThresh = vecCellThresh;
sData.vecTauPeakByType = vecTauPeakByType;
sData.vecCellV_E = vecCellV_E;
sData.vecCellV_I = vecCellV_I;
sData.vecCellV_AHP = vecCellV_AHP;
sData.vecCellV_Leak = vecCellV_Leak;
sData.vecCellCm = vecCellCm;
sData.vecCellG_Leak = vecCellG_Leak;
sData.vecCellG_AHP = vecCellG_AHP;
sData.vecSynConductanceON_to_Cort = vecSynConductanceON_to_Cort;
sData.vecSynConductanceOFF_to_Cort = vecSynConductanceOFF_to_Cort;
sData.vecSynWeightON_to_Cort = vecSynWeightON_to_Cort;
sData.vecSynWeightOFF_to_Cort = vecSynWeightOFF_to_Cort;
sData.vecSynDelayON_to_Cort = vecSynDelayON_to_Cort;
sData.vecSynDelayOFF_to_Cort = vecSynDelayOFF_to_Cort;
sData.matSynConnON_to_Cort = matSynConnON_to_Cort;
sData.matSynConnOFF_to_Cort = matSynConnOFF_to_Cort;
sData.matBlankLGN_ON = matBlankLGN_ON;
sData.matBlankLGN_OFF = matBlankLGN_OFF;
sData.cellLGN_ON = cellLGN_ON;
sData.cellLGN_OFF = cellLGN_OFF;
sData.vecThisV = vecThisV;
sData.boolStimPresent = double(boolStimPresent);
sData.intPrevTrial = intPrevTrial;
sData.intTrialT = intTrialT;
sData.intIter = intIter;
sData.cellSpikeTimesLGN_ON = cellSpikeTimesLGN_ON;
sData.cellSpikeTimesLGN_OFF = cellSpikeTimesLGN_OFF;
sData.cellSpikeTimesCortex = cellSpikeTimesCortex;
sData.vecSpikeCounterLGN_ON = vecSpikeCounterLGN_ON;
sData.vecSpikeCounterLGN_OFF = vecSpikeCounterLGN_OFF;
sData.vecSpikeCounterCortex = vecSpikeCounterCortex;
sData.intPreAllocationSize = intPreAllocationSize;

%cell preferences
sData.vecCellArea = vecCellArea;
sData.vecCellTypes = vecCellTypes;
sData.vecPrefPsi = vecPrefPsi; %cell's phase offset
sData.vecPrefOri=vecPrefOri; %cell's preferred orientation
sData.vecPrefSF = vecPrefSF;
sData.vecPrefRF_X = vecPrefRF_X;
sData.vecPrefRF_Y = vecPrefRF_Y;
sData.dblSigmaX = dblSigmaX; %length of gabor response
sData.dblSigmaY = dblSigmaY; %width of gabor response
sData.matPrefGabors = matPrefGabors;
sData.matFieldsV2 = matFieldsV2;

%stim params
sData.sStimParams = sStimParams;
sData.dblStimDur = dblStimDur;

sData.dblTrialDur = dblTrialDur;
sData.varDeltaSyn = varDeltaSyn;
sData.dblVisSpacing = dblVisSpacing;

sData.vecTrialOris = vecTrialOris;
sData.vecTrialOriIdx = vecTrialOriIdx;
sData.vecStimTypeOris = vecStimTypeOris;

sData.vecTrialSFs = vecTrialSFs;
sData.vecTrialSFIdx = vecTrialSFIdx;
sData.vecStimTypeSFs = vecStimTypeSFs;

sData.vecTrialStimType = vecTrialStimType;
sData.vecTrialStartSecs = vecTrialStartSecs;
sData.vecStimStartSecs = vecStimStartSecs;
sData.vecStimStopSecs = vecStimStopSecs;
sData.vecTrialEndSecs = vecTrialEndSecs;

%% run simulation
sData = getSimulationRun(sData);
%sData = getSimulationRun_xArea(sData);

%% remove placeholder spike entries
cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
for intN=1:numel(cellSpikeTimesCortex)
	if numel(cellSpikeTimesCortex{intN}) > 1
		cellSpikeTimesCortex{intN} = cellSpikeTimesCortex{intN}(~isnan(cellSpikeTimesCortex{intN}));
	else
		cellSpikeTimesCortex{intN} = [];
	end
end
sData.cellSpikeTimesCortex = cellSpikeTimesCortex;

%% save data
save(['D:\Data\Processed\V1_LIFmodel\Simulation_'  getFlankedBy(mfilename,'_',[]) sStimParams.strStimType sprintf('%s.mat',getDate)],'sData','sStimInputs','sStimParams');

%end
%{
%plot spike times
figure
hold on
for intN=1:intCortexCells
	if numel(cellSpikeTimesCortex{intN}) > 0
		line([cellSpikeTimesCortex{intN};cellSpikeTimesCortex{intN}],repmat([intN-0.5;intN+0.5],[1 numel(cellSpikeTimesCortex{intN})]),'Color','k')
	end
end
%}
