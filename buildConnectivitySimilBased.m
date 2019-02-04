function sConnectivity = buildConnectivitySimilBased(sConnParams)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%% check if we use non-uniform weights for cortical connections
	boolUseWeights = sConnParams.boolUseWeights;
	boolUseRFs = sConnParams.boolUseRFs;

	%% msg
	printf('Starting creation of connectivity structure [%s]\n',getTime);
	
	%% extract data
	intImX = sConnParams.vecSizeInput(1);
	intImY = sConnParams.vecSizeInput(2);
	
	vecDefinitionV1PrefOri = sConnParams.vecDefinitionV1PrefOri;
	vecDefinitionV1SpatFreq = sConnParams.vecDefinitionV1SpatFreq;
	vecDefinitionV1CellTypes = sConnParams.vecDefinitionV1CellTypes;
	dblVisSpacing = sConnParams.dblVisSpacing;
	
	%% get LGN to V1 parameters
	%get spatial receptive fields
	vecSpaceX = dblVisSpacing*((-(intImX - 1)/2):(intImX - 1)/2);
	vecSpaceY = dblVisSpacing*((-(intImY - 1)/2):(intImY - 1)/2);
	[matMeshX,matMeshY] = meshgrid(vecSpaceX,vecSpaceY);
	
	%create preferred stimulus features for all cells in V1
	intCellsV1 = sConnParams.intCellsV1;
	
	%distribute RF center locations
	intCellsPerSpatialDim = ceil(sqrt(intCellsV1));
	vecPrefDistroX = linspace(min(vecSpaceX)/3,max(vecSpaceX)/3,intCellsPerSpatialDim);
	vecPrefDistroY = linspace(min(vecSpaceY)/3,max(vecSpaceY)/3,intCellsPerSpatialDim);
	[matPrefRF_X,matPrefRF_Y] = meshgrid(vecPrefDistroX,vecPrefDistroY);
	vecRandRF = randperm(intCellsPerSpatialDim.^2,intCellsV1);
	dblSpatialStep = min(diff(vecPrefDistroX));
	vecPrefRF_X_V1 = matPrefRF_X(vecRandRF)+(dblSpatialStep)*(rand(1,intCellsV1)-0.5);
	vecPrefRF_Y_V1 = matPrefRF_Y(vecRandRF)+(dblSpatialStep)*(rand(1,intCellsV1)-0.5);
	matPrefGabors = nan(numel(vecSpaceX),numel(vecSpaceY),intCellsV1);
	if ~boolUseRFs
		vecPrefRF_X_V1(:) = 0;
		vecPrefRF_Y_V1(:) = 0;
	end
	
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
	vecConnsPerTypeON = sConnParams.vecConnsPerTypeON; %[pyramid interneuron]
	vecConnsPerTypeOFF = sConnParams.vecConnsPerTypeOFF; %[pyramid interneuron]
	
	dblSigmaX = sConnParams.dblSigmaX; %length of gabor response
	dblSigmaY = sConnParams.dblSigmaY; %width of gabor response
	intCellsV1 = length(vecCellTypesV1);
	vecCellAreaV1 = ones(1,intCellsV1); %1, V1; 2, V2
	vecPrefPsiV1 = 2*pi*rand(1,intCellsV1);%phase offset
	intConnsLGN_to_CortOFF = sum(vecConnsPerTypeOFF(vecDefinitionV1CellTypes))*(numel(vecDefinitionV1PrefOri) * numel(vecDefinitionV1SpatFreq));
	intConnsLGN_to_CortON = sum(vecConnsPerTypeON(vecDefinitionV1CellTypes))*(numel(vecDefinitionV1PrefOri) * numel(vecDefinitionV1SpatFreq));
	
	%create connection parameters from LGN to cortex
	vecConductance_FromLGN_ToCort = sConnParams.vecConductance_FromLGN_ToCort; %to [pyramid interneuron]
	vecMeanSynDelayFromLGN_ToCort = sConnParams.vecMeanSynDelayFromLGN_ToCort; %to [pyramid interneuron]
	vecSDSynDelayFromLGN_ToCort = sConnParams.vecSDSynDelayFromLGN_ToCort; %to [pyramid interneuron]
	
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
		dblThetaRot=vecPrefOriV1(intNeuron); %cell's preferred orientation
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
		dblTheta=dblThetaRot;
		matX_theta=(matMeshX*cos(dblTheta)+matMeshY*sin(dblTheta))+dblPrefRF_X;
		matY_theta=(-matMeshX*sin(dblTheta)+matMeshY*cos(dblTheta))+dblPrefRF_Y;
		
		%get gabor
		matG_rot = exp(-.5*(matX_theta.^2/dblSigmaX^2+matY_theta.^2/dblSigmaY^2)).*cos(2*pi*dblPrefSF*matX_theta+dblPsi);
		%matG_rot = imrotate(matG,rad2ang(dblThetaRot),'bilinear','crop');
		
		
		
		matPrefGabors(:,:,intNeuron) = matG_rot;
		
		%get connection probabilities
		matProbConn_ON = max(0,matG_rot) / sum(sum(max(0,matG_rot)));
		matProbConn_OFF = max(0,-matG_rot) / sum(sum(max(0,-matG_rot)));
		
		%{
	%% plot
	subplot(2,3,1)
	%imagesc(matG);
	subplot(2,3,2)
	imagesc(matG_rot);
	subplot(2,3,4)
	imagesc(matProbConn_ON);
	subplot(2,3,5)
	imagesc(matProbConn_OFF);
	title(sprintf('Neuron %d',intNeuron));
		pause
		%}
		
		%% ON fields
		%ON; linearize, get random connections, and reshape to matrix
		vecProbLinON = matProbConn_ON(:)/sum(matProbConn_ON(:));
		vecProbCumSumON = cumsum(vecProbLinON);
		%matConnON = false(size(matProbConn_ON));
		vecConnsON = [];
		intLoopTracker = 0;
		while numel(vecConnsON) < intConnsON && intLoopTracker < 20
			vecConnsON = unique([vecConnsON sum(bsxfun(@gt,rand(1,intConnsON),vecProbCumSumON),1)+1]);
			intLoopTracker = intLoopTracker + 1;
		end
		if numel(vecConnsON) < intConnsON
			[vecConnProbSorted,vecReorder] = sort(vecProbLinON,'descend');
			vecConnsON = vecReorder(1:intConnsON)';
		end
		if numel(vecConnsON) > intConnsON
			vecConnsON = vecConnsON(randperm(numel(vecConnsON),intConnsON));
		end%[I,J] = ind2sub([intImX intImY],vecConnsON);
		%matConnON = logical(getFillGrid(zeros(intImX,intImY),I,J,true));
		
		%put in 2D connection matrix
		matConnections2D = [vecConnsON' repmat(intNeuron,[length(vecConnsON) 1])];
		matSynConnON_to_Cort(vecConnsAssignON,:) = matConnections2D;
		
		%put synaptic properties in vector
		vecSynWeightON_to_Cort(vecConnsAssignON) = vecProbLinON(vecConnsON)./(mean(vecProbLinON(vecConnsON)));
		vecSynConductanceON_to_Cort(vecConnsAssignON) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignON));
		vecSynDelayON_to_Cort(vecConnsAssignON) = max(0,gaussrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignON)));
		
		%% OFF fields
		%OFF; linearize, get random connections
		vecProbLinOFF = matProbConn_OFF(:)/sum(matProbConn_OFF(:));
		vecProbCumSumOFF = cumsum(vecProbLinOFF);
		%matConnOFF = false(size(matProbConn_OFF));
		vecConnsOFF = [];
		intLoopTracker = 0;
		while numel(vecConnsOFF) < intConnsOFF && intLoopTracker < 20
			vecConnsOFF = unique([vecConnsOFF sum(bsxfun(@gt,rand(1,intConnsOFF),vecProbCumSumOFF),1)+1]);
			intLoopTracker = intLoopTracker + 1;
		end
		if numel(vecConnsOFF) < intConnsOFF
			[vecConnProbSorted,vecReorder] = sort(vecProbLinOFF,'descend');
			vecConnsOFF = vecReorder(1:intConnsOFF)';
		end
		if numel(vecConnsOFF) > intConnsOFF
			vecConnsOFF = vecConnsOFF(randperm(numel(vecConnsOFF),intConnsOFF));
		end
		
		%[I,J] = ind2sub([intImX intImY],vecConnsOFF);
		%matConnOFF = getFillGrid(zeros(intImX,intImY),I,J,true);
		
		%put in 2D connection matrix
		matConnections2D = [vecConnsOFF' repmat(intNeuron,[length(vecConnsOFF) 1])];
		matSynConnOFF_to_Cort(vecConnsAssignOFF,:) = matConnections2D;
		
		%put synaptic properties in vector
		vecSynWeightOFF_to_Cort(vecConnsAssignOFF) = vecProbLinOFF(vecConnsOFF)./(mean(vecProbLinOFF(vecConnsOFF)));
		vecSynConductanceOFF_to_Cort(vecConnsAssignOFF) = vecConductance_FromLGN_ToCort(intCellType)*ones(size(vecConnsAssignOFF));
		vecSynDelayOFF_to_Cort(vecConnsAssignOFF) = max(0,gaussrnd(vecMeanSynDelayFromLGN_ToCort(intCellType),vecSDSynDelayFromLGN_ToCort(intCellType),size(vecConnsAssignOFF)));
		
	end
	printf(' .. Created LGN to V1 connectivity; %d active OFF field synapses, and %d active ON field synapses [%s]\n',intConnsLGN_to_CortOFF,intConnsLGN_to_CortON,getTime);
	
	%% define LIF model parameters
	%% neuron parameters
	%{
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
	%}
	sParamCells{1}.dblThresh = -55; %mV
	sParamCells{1}.dblTau_Peak = 1/1000; %ms
	sParamCells{1}.dblV_E = 0; %excitatory reversal potential
	sParamCells{1}.dblV_I = -70; %inhibitory reversal potential
	sParamCells{1}.dblV_AHP = -90;%after-hyperpol reversal potential
	sParamCells{1}.dblV_Leak = -65; %resting potential
	sParamCells{1}.dblCm = 0.5; %capacitance
	sParamCells{1}.dblG_Leak = 25; %nS
	sParamCells{1}.dblG_AHP = 40; %nS
	sParamCells{1}.dblRefracT = 0.5/1000; %seconds
	

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
	sParamCells{2}.dblRefracT = 0.5/1000; %seconds
	
	vecTauPeakByType = nan(1,numel(sParamCells));
	for intType=1:numel(sParamCells)
		vecTauPeakByType(intType) = sParamCells{intType}.dblTau_Peak;
	end
	
	%% connection parameters
	matConnCortFromTo = sConnParams.matConnCortFromTo; %from pyramid to [pyr inter]
	
	%conductances
	matConductancesFromTo = sConnParams.matConductancesFromTo; %from pyramid to [pyr inter]
	
	%synaptic delays
	dblDelayMeanCortToCort = sConnParams.dblDelayMeanCortToCort; %in ms
	dblDelaySDCortToCort = sConnParams.dblDelaySDCortToCort; %in ms
	
	%msg
	printf(' .. Created LIF model parameters [%s]\n',getTime);
	
	%% calculate similarity between receptive fields of V1 cells
	printf(' .. Building V1 field similarities for %dx%d cells... [%s]\n',intCellsV1,intCellsV1,getTime);
	matTuningSimilarity=nan(intCellsV1,intCellsV1);
	parfor i1=1:intCellsV1
		a = matPrefGabors(:,:,i1);
		a = a - (sum(a(:),'double') / numel(a));
		for i2=1:intCellsV1
			b = matPrefGabors(:,:,i2); %#ok<PFBNS>
			b = b - (sum(b(:),'double') / numel(b));
			r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
			matTuningSimilarity(i1,i2) = r;
		end
	end
	matTuningSimilarity(diag(diag(true(intCellsV1)))) = 0;
	
	%% create cortical connectivity and parameter vectors
	%locality of connection probability between cells for different types
	vecLocalityLambda = sConnParams.vecLocalityLambda;
	printf(' .. Building V1 connectivity; locality parameters are: [%s\b] [%s]\n',sprintf('%.2f ',vecLocalityLambda),getTime);
	
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
	vecCellRefracT = nan(intCellsV1,1);
	
	%get random connections
	intC = 1;
	for intNeuron=1:intCellsV1
		%dblPrefOri = vecPrefOriV1(intNeuron);
		%vecOriDiff = abs(circ_dist(vecPrefOriV1*2,dblPrefOri*2)/2)';
		
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
		vecCellRefracT(intNeuron) = sParamCells{intTargetCellType}.dblRefracT; %secs
		
		
		%loop through subtypes to make connections
		for intSubIdx=1:length(vecCellTypeLabels)
			intSourceCellType = vecCellTypeLabels(intSubIdx);
			indThisType = vecCellTypesV1==intSourceCellType;
			intConnsFromThisType = vecConnsFromTypesToThisType(intSourceCellType);
			dblLocalityLambda = vecLocalityLambda(intSourceCellType);
			%%
			intChooseSimN = round(intConnsFromThisType*(1-abs(dblLocalityLambda)));
			vecSims = matTuningSimilarity(:,intNeuron);
			vecSims(vecSims < 0) = 0; %remove anti-correlated neurons
			vecSims(~indThisType) = 0; %remove other types
			vecSims(intNeuron) = 0; %remove the neuron itself
			vecConnProb = vecSims/sum(vecSims);
			vecProbCumSum = cumsum(vecConnProb);
			
			%get assignment vectors
			vecConnsAssign = intC:(intC+intConnsFromThisType-1);
			intC = intC + intConnsFromThisType;
			
			% get sim-based connections
			%make connections
			vecConns = [];
			intLoopTracker = 0;
			while numel(vecConns) < intChooseSimN && intLoopTracker < 20
				vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intChooseSimN),vecProbCumSum),1)+1]);
				intLoopTracker = intLoopTracker + 1;
			end
			if numel(vecConns) < intChooseSimN
				error
			end
			if numel(vecConns) > intChooseSimN
				vecConns = vecConns(randperm(numel(vecConns),intChooseSimN));
			end
			
			
			% get other connections
			%get required number of connections
			intChooseOtherN = intConnsFromThisType - numel(vecConns);
			if dblLocalityLambda > 0 % local
				%set chosen connections to 0
				vecConnProb(vecConns) = 0;
				[vecConnProbSorted,vecReorder] = sort(vecConnProb,'descend');
				vecConns = [vecConns vecReorder(1:intChooseOtherN)'];
			else %non-local
				indPossibleConnections = indThisType;
				indPossibleConnections(vecConns) = 0;
				indPossibleConnections(intNeuron) = 0; %remove the neuron itself
				vecPossConns = find(indPossibleConnections);
				vecConns = [vecConns vecPossConns(randperm(numel(vecPossConns),intChooseOtherN))];
			end
			%hist(vecConns)
			
			%% assign to connection matrix
			vecTheseWeights = vecConnProb(vecConns);
			vecTheseWeights = (vecTheseWeights./mean(vecTheseWeights));
			dblG = matConductancesFromTo(intSourceCellType,intTargetCellType);
			matSynV1FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
			vecSynV1Conductance(vecConnsAssign) = dblG;
			if boolUseWeights
				vecSynV1Weight(vecConnsAssign) = vecTheseWeights;
			else
				vecSynV1Weight(vecConnsAssign) = 1;
			end
			vecSynV1Delay(vecConnsAssign) = max(0,gaussrnd(dblDelayMeanCortToCort,dblDelaySDCortToCort,size(vecConnsAssign)));
			vecSynV1ExcInh(vecConnsAssign) = intSourceCellType;
			vecSynV1Type(vecConnsAssign) = 1; %1; intra-areal, 2 inter-areal
		end
	end
	printf(' .. Created V1 connectivity; %d active corticocortical synapses between %d cells [%s]\n',intTotalV1Connections,intCellsV1,getTime);
	
	%% build V1-V2 connectivity
	%set parameters for V2
	%dblPercProjectingV1_to_V2 = 10; %El-Shamayleh et al (Movshon), 2013;
	dblSpatialDropoffV1V2 = sConnParams.dblSpatialDropoffV1V2; %gausspdf(vecX,0,0.8); zandvakili&kohn, 2015
	dblSpatialDropoffInterneuronsV2 = sConnParams.dblSpatialDropoffInterneuronsV2; %for interneurons
	intCellsV2 = sConnParams.intCellsV2;
	if intCellsV2>0
		%create cell-based parameters
		vecDefinitionV2CellTypes = sConnParams.vecDefinitionV2CellTypes;
		vecCellTypesV2 = sConnParams.vecCellTypesV2;
		
		%LGN projects to V2 (Bullier & Kennedy, 1983), but how much?
		
		%Schmid et al. (Smirnakis), 2009; 70-95% reduction in V2 activity after lesion of V1
		
		%pulvinar plays major role in projecting to V2, and received most of its input from V1 (Sincich & Horton, 2005)
		%V2<->V1; 14,000 feedforward and 11,000 feedback connections per mm^2; if FF connections are 100%, FB are 80%  (Sincich & Horton, 2005; Rockland, 1997)
		%connections from V1->V2 are 20-25 fold larger than LGN->V1 (van Essen, 2005)
		
		%V2 cells respond to angles and features; simple model: combining two V1
		%neurons with similar receptive field locations (Ito and Komatsu, 2004)
		
		
		%% create RF exemplars for V2 cells
		printf(' .. Building V2 exemplar fields for %d cells... [%s]\n',intCellsV2,getTime);
		%transform r
		matTotalSim = matTuningSimilarity;
		matTotalSim = abs(tril(matTotalSim,-1));
		matTotalSim(vecCellTypesV1==2,:) = 0; %remove interneurons
		matTotalSim(:,vecCellTypesV1==2) = 0; %remove interneurons
		matTotalSim(matTotalSim==0)=nan;
		matTotalSim = matTotalSim.^2;% + matTuningSimilarity;
		
		%select r
		matTotalSim(isnan(matTotalSim))=0;
		vecProbComb=(abs(matTotalSim(:))./sum(abs(matTotalSim(:))));
		[vecProbSorted,vecSort]=sort(vecProbComb);
		vecProbSorted(isnan(vecProbSorted))=0;
		vecSelectProb = cumsum(vecProbSorted)/sum(vecProbSorted);
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
		printf(' .. Calculating similarity of V1 Gabors to V2 exemplars... [%s]\n',getTime);
		matSimilarityGaborExamplar = nan(intCellsV1,intCellsV2);
		parfor intCellV1=1:intCellsV1
			a = matPrefGabors(:,:,intCellV1);
			a = a - (sum(a(:),'double') / numel(a));
			for intCellV2=1:intCellsV2
				b = matExemplarFieldsV2(:,:,intCellV2); %#ok<PFBNS>
				b = b - (sum(b(:),'double') / numel(b));
				r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
				matSimilarityGaborExamplar(intCellV1,intCellV2) = r;
			end
		end
		%msg
		printf(' .. Created V2 fields, building V1 => V2 connectivity... [%s]\n',getTime);
		
		
		%% build preferred RF location for interneurons
		%distribute RF center locations
		vecCellFractionsV2 = sConnParams.vecCellFractionsV2; %[pyramid interneuron]
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
		%connection parameters
		vecConnsPerTypeV1V2 = sConnParams.vecConnsPerTypeV1V2;%[48 32]; %pyr/int
		matConductancesFromToV1V2 = sConnParams.matConductancesFromToV1V2;
		
		%synaptic delays
		dblDelayMeanV1ToV2 = sConnParams.dblDelayMeanV1ToV2; %in ms
		dblDelaySDV1ToV2 = sConnParams.dblDelaySDV1ToV2; %in ms
		
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
		vecCellV2RefracT = nan(intCellsV2,1);
		
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
			vecCellV2RefracT(intCellV2) = sParamCells{intTargetCellType}.dblRefracT; %secs
			
			%msg
			%printf('Building V1=>V2 connections for V2 %sneuron %d/%d [%s]\n',cellStrType{intTargetCellType},intCellV2,intCellsV2,getTime);
			
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
				vecConnProb = gausspdf(vecDistanceRF',0,dblSpatialDropoffV1V2);
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
			intLoopTracker = 0;
			while numel(vecConns) < intConnsToThisType && intLoopTracker < 20
				vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsToThisType),vecProbCumSum),1)+1]);
				intLoopTracker = intLoopTracker + 1;
			end
			if numel(vecConns) < intConnsToThisType
				[vecConnProbSorted,vecReorder] = sort(vecConnProb,'descend');
				vecConns = vecReorder(1:intConnsToThisType)';
			end
			if numel(vecConns) > intConnsToThisType
				vecConns = vecConns(randperm(numel(vecConns),intConnsToThisType));
			end
			
			%assign to connection matrix
			vecTheseWeights = vecConnProb(vecConns);
			vecTheseWeights = (vecTheseWeights./mean(vecTheseWeights));
			dblG = matConductancesFromToV1V2(intSourceCellType,intTargetCellType);
			matSynV1V2FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
			vecSynV1V2Conductance(vecConnsAssign) = dblG; %synaptic conductance
			if boolUseWeights
				vecSynV1V2Weight(vecConnsAssign) = vecTheseWeights; %synaptic weights
			else
				vecSynV1V2Weight(vecConnsAssign) = 1;
			end
			vecSynV1V2Delay(vecConnsAssign) = max(0,gaussrnd(dblDelayMeanV1ToV2,dblDelaySDV1ToV2,size(vecConnsAssign))); %synaptic delay
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
		printf(' .. Created V1 => V2 connectivity; %d active inter-areal synapses from %d V1 to %d V2 cells [%s]\n',intTotalV1V2Connections,intCellsV1,intCellsV2,getTime);
		
		%% build RF similarities V2
		vecPrefRF_X_V2 = nan(1,intCellsV2);
		vecPrefRF_Y_V2 = nan(1,intCellsV2);
		matSimilarityFieldsV2 = zeros(intCellsV2,intCellsV2);
		parfor intCellV2_1=1:intCellsV2
			vecCoM = calcCenterOfMass(abs(matFieldsV2(:,:,intCellV2_1)));
			vecPrefRF_X_V2(intCellV2_1) = vecCoM(2);
			vecPrefRF_Y_V2(intCellV2_1) = vecCoM(1);
			
			a = matFieldsV2(:,:,intCellV2_1);
			a = a - (sum(a(:),'double') / numel(a));
			for intCellV2_2=1:intCellsV2
				b = matFieldsV2(:,:,intCellV2_2); %#ok<PFBNS>
				b = b - (sum(b(:),'double') / numel(b));
				r = sum(sum(a.*b))/sqrt(sum(sum(a.*a))*sum(sum(b.*b)));
				
				matSimilarityFieldsV2(intCellV2_1,intCellV2_2) = r;
			end
		end
		vecPrefRF_X_V2 = dblVisSpacing*(vecPrefRF_X_V2 - intImX/2);
		vecPrefRF_Y_V2 = dblVisSpacing*(vecPrefRF_Y_V2 - intImY/2);
		matSimilarityFieldsV2(matSimilarityFieldsV2<0)=0;
		matSimilarityFieldsV2(diag(diag(true(size(matSimilarityFieldsV2))))) = 0;
		
		%msg
		printf(' .. Created V2 field similarities, building V2 => V2 connectivity... [%s]\n',getTime);
		
		
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
					vecConnProb = gausspdf(vecDistanceRF',0,dblSpatialDropoffInterneuronsV2);
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
				intLoopTracker = 0;
				while numel(vecConns) < intConnsFromThisType && intLoopTracker < 20
					vecConns = unique([vecConns sum(bsxfun(@gt,rand(1,intConnsFromThisType),vecProbCumSum),1)+1]);
					intLoopTracker = intLoopTracker + 1;
				end
				if numel(vecConns) < intConnsFromThisType
					[vecConnProbSorted,vecReorder] = sort(vecConnProb,'descend');
					vecConns = vecReorder(1:intConnsFromThisType)';
				end
				if numel(vecConns) > intConnsFromThisType
					vecConns = vecConns(randperm(numel(vecConns),intConnsFromThisType));
				end
				
				%assign to connection matrix
				vecTheseWeights = vecConnProb(vecConns);
				vecTheseWeights = (vecTheseWeights./mean(vecTheseWeights));
				vecConns = vecConns + intCellsV1;
				dblG = matConductancesFromTo(intSourceCellType,intTargetCellType);
				matSynV2FromTo(vecConnsAssign,:) = [vecConns' repmat(intNeuron,[length(vecConns) 1])]; %[From x To]
				vecSynV2Conductance(vecConnsAssign) = dblG;
				if boolUseWeights
					vecSynV2Weight(vecConnsAssign) = vecTheseWeights; %synaptic weights
				else
					vecSynV2Weight(vecConnsAssign) = 1;
				end
				vecSynV2Delay(vecConnsAssign) = max(0,gaussrnd(dblDelayMeanCortToCort,dblDelaySDCortToCort,size(vecConnsAssign)));
				vecSynV2ExcInh(vecConnsAssign) = intSourceCellType;
				vecSynV2Type(vecConnsAssign) = 1; %1; intra-areal, 2 inter-areal
			end
		end
		printf(' .. Created V2 connectivity; %d active corticocortical synapses between %d cells [%s]\n',intTotalV2Connections,intCellsV2,getTime);
		
		%% V2 => V1
		if isfield(sConnParams,'matConnsPerTypeV2V1')
			printf(' .. Building V2=>V1 feedback connections... [%s]\n',getTime);
			
			%connection parameters
			dblSpatialDropoffV2V1 = sConnParams.dblSpatialDropoffV2V1; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
			matConnsPerTypeV2V1 = sConnParams.matConnsPerTypeV2V1;%[48 32]; %from pyr to pyr/int
			matConductancesFromToV2V1 = sConnParams.matConductancesFromToV2V1;
			
			%synaptic delays
			dblDelayMeanV2ToV1 = sConnParams.dblDelayMeanV2ToV1; %in ms
			dblDelaySDV2ToV1 = sConnParams.dblDelaySDV2ToV1; %in ms
			
			
			%pre-allocate cortical synaptic connection matrix
			vecTotalConnectionsPerTypeV2V1 = sum(matConnsPerTypeV2V1,1);
			intTotalV2V1Connections = sum(vecTotalConnectionsPerTypeV2V1(vecCellTypesV1));
			matSynV2V1FromTo = nan(intTotalV2V1Connections,2); %[From x To]
			vecSynV2V1Conductance = nan(intTotalV2V1Connections,1); %synaptic conductance
			vecSynV2V1Weight = ones(intTotalV2V1Connections,1); %synaptic weight
			vecSynV2V1Delay = nan(intTotalV2V1Connections,1); %synaptic delay
			vecSynV2V1ExcInh = nan(intTotalV2V1Connections,1); %synaptic sign; inhibitory or excitatory
			vecSynV2V1Type = 2*ones(intTotalV2V1Connections,1); %synaptic type; intra- or inter-areal
			
			
			%get random connections
			intC = 1;
			for intCellV1=1:intCellsV1
				%get this cell's type
				intTargetCellType = vecCellTypesV1(intCellV1); %pyramidal cell or interneuron
				vecConnsFromTypesToThisType = matConnsPerTypeV2V1(:,intTargetCellType); %from [pyr inter]
				
				%get cell parameters
				dblTargetPrefRF_X = vecPrefRF_X_V1(intCellV1);
				dblTargetPrefRF_Y = vecPrefRF_Y_V1(intCellV1);
				
				%loop through subtypes to make connections
				for intSubIdx=1:numel(vecConnsFromTypesToThisType) %only pyramids
					intSourceCellType = vecCellTypeLabels(intSubIdx);
					indThisType = vecCellTypesV2==intSourceCellType;
					intConnsFromThisType = vecConnsFromTypesToThisType(intSourceCellType);
					if intConnsFromThisType == 0,continue;end
					
					% RF center distance based
					vecDistanceRF = sqrt((vecPrefRF_Y_V2-dblTargetPrefRF_Y).^2 + (vecPrefRF_X_V2-dblTargetPrefRF_X).^2); %euclidian distance
					vecConnProb = gausspdf(vecDistanceRF',0,dblSpatialDropoffV2V1);
					vecConnProb(~indThisType) = 0; %remove other types
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
					dblG = matConductancesFromToV2V1(intSourceCellType,intTargetCellType);
					matSynV2V1FromTo(vecConnsAssign,:) = [vecConns' repmat(intCellV1,[length(vecConns) 1])]; %[From x To]
					vecSynV2V1Conductance(vecConnsAssign) = dblG;
					vecSynV2V1Weight(vecConnsAssign) = 1;
					vecSynV2V1Delay(vecConnsAssign) = max(0,gaussrnd(dblDelayMeanV2ToV1,dblDelaySDV2ToV1,size(vecConnsAssign)));
					vecSynV2V1ExcInh(vecConnsAssign) = intSourceCellType;
					vecSynV2V1Type(vecConnsAssign) = 1; %1; intra-areal, 2 inter-areal
				end
			end
			
			%msg
			printf(' .. Created feedback connectivity; %d active corticocortical synapses from %d V2 to %d V1 cells [%s]\n',intTotalV2V1Connections,intCellsV2,intCellsV1,getTime);
			
		else
			%build dummies
			matSynV2V1FromTo = []; %[From x To]
			vecSynV2V1Conductance = []; %synaptic conductance
			vecSynV2V1Weight = []; %synaptic weight
			vecSynV2V1Delay = []; %synaptic delay
			vecSynV2V1ExcInh = []; %synaptic sign; inhibitory or excitatory
			vecSynV2V1Type = []; %synaptic type; intra- or inter-areal
		end
	else
		%build dummies
		matSynV2V1FromTo = []; %[From x To]
		vecSynV2V1Conductance = []; %synaptic conductance
		vecSynV2V1Weight = []; %synaptic weight
		vecSynV2V1Delay = []; %synaptic delay
		vecSynV2V1ExcInh = []; %synaptic sign; inhibitory or excitatory
		vecSynV2V1Type = []; %synaptic type; intra- or inter-areal
		
		vecCellTypesV2 = [];
		vecPrefRF_X_V2 = [];
		vecPrefRF_Y_V2 = [];
		
		matSynV1V2FromTo = []; %[From x To]
		vecSynV1V2Conductance = []; %synaptic conductance
		vecSynV1V2Weight = []; %synaptic weight
		vecSynV1V2Delay = []; %synaptic delay
		vecSynV1V2ExcInh = []; %synaptic sign; inhibitory or excitatory
		vecSynV1V2Type = []; %synaptic type; intra- or inter-areal
		
		matSynV2FromTo = []; %[From x To]
		vecSynV2Conductance = []; %synaptic conductance
		vecSynV2Weight = []; %synaptic weight
		vecSynV2Delay = []; %synaptic delay
		vecSynV2ExcInh = []; %synaptic sign; inhibitory or excitatory
		vecSynV2Type = []; %synaptic type; intra- or inter-areal
		
		vecCellV2Thresh = [];
		vecCellV2TauPeak = [];
		vecCellV2V_E = [];
		vecCellV2V_I = [];
		vecCellV2V_AHP = [];
		vecCellV2V_Leak = [];
		vecCellV2Cm =	[];
		vecCellV2G_Leak = [];
		vecCellV2G_AHP = [];
		vecCellV2RefracT = [];
		
		matFieldsV2 = [];
	end
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
	matSynFromTo = cat(1,matSynV1FromTo,matSynV1V2FromTo,matSynV2FromTo,matSynV2V1FromTo);
	vecSynConductance = cat(1,vecSynV1Conductance,vecSynV1V2Conductance,vecSynV2Conductance,vecSynV2V1Conductance);
	vecSynWeight = cat(1,vecSynV1Weight,vecSynV1V2Weight,vecSynV2Weight,vecSynV2V1Weight);
	vecSynDelay = cat(1,vecSynV1Delay,vecSynV1V2Delay,vecSynV2Delay,vecSynV2V1Delay);
	vecSynExcInh = cat(1,vecSynV1ExcInh,vecSynV1V2ExcInh,vecSynV2ExcInh,vecSynV2V1ExcInh);
	vecSynType = cat(1,vecSynV1Type,vecSynV1V2Type,vecSynV2Type,vecSynV2V1Type);
	
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
	%vecCellRefracT = cat(1,vecCellRefracT,vecCellV2RefracT);
	
	%% combine data
	sConnectivity = struct;
	%summary variables
	sConnectivity.intCellsV2 = intCellsV2;
	sConnectivity.intCellsV1 = intCellsV1;
	sConnectivity.intCortexCells = intCortexCells;
	
	%LGN->cortex connectivity
	sConnectivity.vecSynConductanceON_to_Cort = vecSynConductanceON_to_Cort;
	sConnectivity.vecSynConductanceOFF_to_Cort = vecSynConductanceOFF_to_Cort;
	sConnectivity.vecSynWeightON_to_Cort = vecSynWeightON_to_Cort;
	sConnectivity.vecSynWeightOFF_to_Cort = vecSynWeightOFF_to_Cort;
	sConnectivity.vecSynDelayON_to_Cort = vecSynDelayON_to_Cort;
	sConnectivity.vecSynDelayOFF_to_Cort = vecSynDelayOFF_to_Cort;
	sConnectivity.matSynConnON_to_Cort = matSynConnON_to_Cort;
	sConnectivity.matSynConnOFF_to_Cort = matSynConnOFF_to_Cort;
	
	%cortical connectivity
	sConnectivity.matSynFromTo = matSynFromTo;
	sConnectivity.vecSynExcInh = vecSynExcInh;
	sConnectivity.vecSynDelay = vecSynDelay;
	sConnectivity.vecSynWeight = vecSynWeight;
	sConnectivity.vecSynConductance = vecSynConductance;
	sConnectivity.vecSynType = vecSynType;
	
	%cortical cell parameters
	sConnectivity.vecCellThresh = vecCellThresh;
	sConnectivity.vecCellTauPeak = vecCellTauPeak;
	sConnectivity.vecCellV_E = vecCellV_E;
	sConnectivity.vecCellV_I = vecCellV_I;
	sConnectivity.vecCellV_AHP = vecCellV_AHP;
	sConnectivity.vecCellV_Leak = vecCellV_Leak;
	sConnectivity.vecCellCm = vecCellCm;
	sConnectivity.vecCellG_Leak = vecCellG_Leak;
	sConnectivity.vecCellG_AHP = vecCellG_AHP;
	%sConnectivity.vecCellRefracT = vecCellRefracT;
	sConnectivity.sParamCells = sParamCells;
	sConnectivity.vecTauPeakByType = vecTauPeakByType;
	
	%cell preferences
	sConnectivity.vecCellArea = vecCellArea;
	sConnectivity.vecCellTypes = vecCellTypes;
	sConnectivity.vecPrefPsi = vecPrefPsi; %cell's phase offset
	sConnectivity.vecPrefOri=vecPrefOri; %cell's preferred orientation
	sConnectivity.vecPrefSF = vecPrefSF;
	sConnectivity.vecPrefRF_X = vecPrefRF_X;
	sConnectivity.vecPrefRF_Y = vecPrefRF_Y;
	sConnectivity.dblSigmaX = dblSigmaX; %length of gabor response
	sConnectivity.dblSigmaY = dblSigmaY; %width of gabor response
	sConnectivity.matPrefGabors = matPrefGabors;
	sConnectivity.matFieldsV2 = matFieldsV2;
	%end
	
