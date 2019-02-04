%% msg
fprintf('Starting offline construction of connectivity profile... [%s]\n\n',getTime);

%% set connectivity parameters, or load pre-existing connectivity structure
intColumns = 32; %48 / 252
sConnParams = struct;
sConnParams.dblVisSpacing = 0.2;
sConnParams.vecSizeInput = [64 64];

%connection definition LGN
sConnParams.vecConnsPerTypeON = [24 16]; %[pyramid interneuron]
sConnParams.vecConnsPerTypeOFF = [24 16]; %[pyramid interneuron]

sConnParams.dblSigmaX = 0.7; %length of gabor response
sConnParams.dblSigmaY = 0.47; %width of gabor response
sConnParams.vecConductance_FromLGN_ToCort = [5.5 6]*0.18; %to [pyramid interneuron]
sConnParams.vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
sConnParams.vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

%V1 def
sConnParams.vecDefinitionV1PrefOri = 0:pi/intColumns:(pi-pi/intColumns);
sConnParams.vecDefinitionV1SpatFreq = 2.^[-4:1:1];
sConnParams.vecDefinitionV1CellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
sConnParams.intCellsV1 = numel(sConnParams.vecDefinitionV1PrefOri) * numel(sConnParams.vecDefinitionV1SpatFreq) * numel(sConnParams.vecDefinitionV1CellTypes);
	
%cortical connectivity
%number of connections
dblScalingFactor = 3; %4
sConnParams.matConnCortFromTo(1,:) = [30 30]*dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConnCortFromTo(2,:) = [40 40]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
sConnParams.matConductancesFromTo(1,:) = [1.1 1.6]/dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromTo(2,:) = [2.0 1.0]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanCortToCort = 3/1000; %in ms
sConnParams.dblDelaySDCortToCort = 1/1000; %in ms

%connection probability excitatory cells
sConnParams.vecConnProbSD = ang2rad([7.5 80]); %degrees difference in pref ori, for [pyramid interneuron]

%V2 params
sConnParams.dblSpatialDropoffV1V2 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
sConnParams.dblSpatialDropoffInterneuronsV2 = 3; %for interneurons
sConnParams.intCellsV2 = round(sConnParams.intCellsV1);

%create cell-based parameters
sConnParams.vecDefinitionV2CellTypes = [1 1 1 1 2];
sConnParams.vecCellTypesV2 = repmat(sConnParams.vecDefinitionV2CellTypes,[1 ceil(sConnParams.intCellsV2/numel(sConnParams.vecDefinitionV2CellTypes))]);
sConnParams.vecCellTypesV2(sConnParams.intCellsV2+1:end) = [];
sConnParams.vecCellFractionsV2 = [0.8 0.2]; %[pyramid interneuron]
	
%V1=>V2
dblInterArealFactor = 2;
sConnParams.vecConnsPerTypeV1V2 = [60 50];%[48 32]; %pyr/int
sConnParams.matConductancesFromToV1V2(1,:) = [1.1 1.6]*dblInterArealFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromToV1V2(2,:) = [1.5 1.0]*dblInterArealFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanV1ToV2 = 4/1000; %in ms
sConnParams.dblDelaySDV1ToV2 = 1/1000; %in ms

%V2=>V1
dblInterArealFeedbackFactor = 1;
sConnParams.dblSpatialDropoffV2V1 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
sConnParams.matConnsPerTypeV2V1(1,:) = [24 24];%from pyramid to [pyr inter]
sConnParams.matConnsPerTypeV2V1(2,:) = [0 0];%from interneuron to [pyr inter]
sConnParams.matConductancesFromToV2V1(1,:) = [1.1 1.6]*dblInterArealFeedbackFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromToV2V1(2,:) = [1.5 1.0]*dblInterArealFeedbackFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanV2ToV1 = 4/1000; %in ms
sConnParams.dblDelaySDV2ToV1 = 1/1000; %in ms

%% build connectivity
sConnectivity = buildConnectivity_xArea(sConnParams);

%% save
strConnFile = sprintf('sConn_Col%dN%dS%d_%s.mat',intColumns,sConnectivity.intCortexCells,numel(sConnectivity.vecSynExcInh),getDate);
strConnDir = 'D:\Simulations\Connectivity\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strConnFile,strConnDir,getTime);
save([strConnDir strConnFile],'sConnParams','sConnectivity');


