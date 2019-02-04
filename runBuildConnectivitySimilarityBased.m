%% msg
fprintf('Starting offline construction of connectivity profile... [%s]\n\n',getTime);

%% set connectivity parameters
sConnParams = struct;
sConnParams.intCellsV1 = 3600;
sConnParams.vecSizeInput = [32 32];
sConnParams.dblVisSpacing = 6.4/sConnParams.vecSizeInput(1);

% set connectivity parameters; synaptic weights/RF location variation
sConnParams.boolUseWeights = false;
sConnParams.boolUseRFs = true;
sConnParams.boolUseSFs = true;

%connection definition LGN
sConnParams.vecConnsPerTypeON = 2*[36 24]; %[pyramid interneuron]
sConnParams.vecConnsPerTypeOFF = 2*[36 24]; %[pyramid interneuron]

sConnParams.dblSigmaX = 0.7; %length of gabor response
sConnParams.dblSigmaY = 0.7; %width of gabor response
sConnParams.vecConductance_FromLGN_ToCort = [1.1 1.2]*0.24; %to [pyramid interneuron]
sConnParams.vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
sConnParams.vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

%V1 def
if sConnParams.boolUseSFs
	sConnParams.vecDefinitionV1SpatFreq = 2.^[-3 -2 -1 0];%2.^[-3:1:1];
else
	sConnParams.vecDefinitionV1SpatFreq = 2.^[-2];%2.^[-3:1:1];
end
sConnParams.vecDefinitionV1CellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
sConnParams.intColumns = sConnParams.intCellsV1 / (numel(sConnParams.vecDefinitionV1SpatFreq) * numel(sConnParams.vecDefinitionV1CellTypes)); %48 / 252 / 120 / 600
if ~isint(sConnParams.intColumns),error([mfilename ':ColumnsNotInteger'],'Number of cells requested is not divisable by mumber of tuning property combinations');end
sConnParams.vecDefinitionV1PrefOri = 0:pi/sConnParams.intColumns:(pi-pi/sConnParams.intColumns);

%cortical connectivity
%number of connections
dblScalingFactor = 4; %4
sConnParams.matConnCortFromTo(1,:) = [40 40]*dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConnCortFromTo(2,:) = [30 30]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
sConnParams.matConductancesFromTo(1,:) = 1*[0.96 1.28]/dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromTo(2,:) = 1*[1.50 1.00]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanCortToCort = 3/1000; %in ms
sConnParams.dblDelaySDCortToCort = 1/1000; %in ms

% locality of connectivity;
%1=most local, 0=proportional to similarity, -1=equal probability for all
sConnParams.vecLocalityLambda = [0 -0.5]; %[pyramid interneuron]

%V2 params
sConnParams.dblSpatialDropoffV1V2 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
sConnParams.dblSpatialDropoffInterneuronsV2 = 3; %for interneurons
sConnParams.intCellsV2 = 1200;%sConnParams.intCellsV1;%round(sConnParams.intCellsV1/2);

%create cell-based parameters
sConnParams.vecDefinitionV2CellTypes = [1 1 1 1 2];
sConnParams.vecCellTypesV2 = repmat(sConnParams.vecDefinitionV2CellTypes,[1 ceil(sConnParams.intCellsV2/numel(sConnParams.vecDefinitionV2CellTypes))]);
sConnParams.vecCellTypesV2(sConnParams.intCellsV2+1:end) = [];
sConnParams.vecCellFractionsV2 = [0.8 0.2]; %[pyramid interneuron]

%V1=>V2
dblInterArealFactor = 2.5; %used to be 2.6 (2.8)
sConnParams.vecConnsPerTypeV1V2 = [48 32];%[48 32]; %from pyr to pyr/int
sConnParams.matConductancesFromToV1V2(1,:) = [1.1 1.6]*dblInterArealFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromToV1V2(2,:) = [1.5 1.0]*dblInterArealFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanV1ToV2 = 4/1000; %in ms
sConnParams.dblDelaySDV1ToV2 = 1/1000; %in ms

%% build connectivity
sConnectivity = buildConnectivitySimilBased(sConnParams);
return
%% save
if sConnParams.boolUseRFs && sConnParams.boolUseSFs && ~sConnParams.boolUseWeights
	strConnFile = sprintf('sConnSimil2_Ret%dCol%dN%dS%d_%s.mat',...
		sConnParams.vecSizeInput(1),sConnParams.intColumns,sConnectivity.intCortexCells,...
		numel(sConnectivity.vecSynExcInh),getDate);
else
	strConnFile = sprintf('sConnSimil_Ret%dLP%.1fLI%.1fCol%dN%dS%dW%dRF%dSF%d_%s.mat',...
		sConnParams.vecSizeInput(1),sConnParams.vecLocalityLambda(1),sConnParams.vecLocalityLambda(2),...
		sConnParams.intColumns,sConnectivity.intCortexCells,numel(sConnectivity.vecSynExcInh),...
		sConnParams.boolUseWeights,sConnParams.boolUseRFs,sConnParams.boolUseSFs,getDate);
end
strConnDir = 'D:\Simulations\Connectivity\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strConnFile,strConnDir,getTime);
save([strConnDir strConnFile],'sConnParams','sConnectivity');


