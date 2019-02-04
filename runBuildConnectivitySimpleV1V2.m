%% msg
fprintf('Starting offline construction of connectivity profile... [%s]\n\n',getTime);

%% set connectivity parameters, or load pre-existing connectivity structure
intColumns = 120; %48 / 252
sConnParams = struct;
sConnParams.dblVisSpacing = 0.1;
sConnParams.vecSizeInput = [64 64];

%connection definition LGN
sConnParams.vecConnsPerTypeON = [24 16]*1.5; %[pyramid interneuron]
sConnParams.vecConnsPerTypeOFF = [24 16]*1.5; %[pyramid interneuron]

sConnParams.dblSigmaX = 0.7; %length of gabor response
sConnParams.dblSigmaY = 0.7; %width of gabor response
sConnParams.vecConductance_FromLGN_ToCort = ([5.5 6]*0.3)/1.5; %to [pyramid interneuron]
sConnParams.vecMeanSynDelayFromLGN_ToCort = [10 5]/1000; %to [pyramid interneuron]
sConnParams.vecSDSynDelayFromLGN_ToCort = [7 3]/1000; %to [pyramid interneuron]

%V1 def
sConnParams.vecDefinitionV1PrefOri = 0:pi/intColumns:(pi-pi/intColumns);
sConnParams.vecDefinitionV1SpatFreq = 2.^[-2 -1];%2.^[-3:1:1];
sConnParams.vecDefinitionV1CellTypes = [1 1 1 1 2]; %[1=pyramid 2=interneuron]
sConnParams.intCellsV1 = numel(sConnParams.vecDefinitionV1PrefOri) * numel(sConnParams.vecDefinitionV1SpatFreq) * numel(sConnParams.vecDefinitionV1CellTypes);
	
%cortical connectivity
%number of connections
dblScalingFactor = 4; %4
sConnParams.matConnCortFromTo(1,:) = [40 40]*dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConnCortFromTo(2,:) = [30 30]*dblScalingFactor; %from interneuron to [pyr inter]

%conductances
sConnParams.matConductancesFromTo(1,:) = 0.8*[1.1 1.6]/dblScalingFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromTo(2,:) = [1.5 1.0]/dblScalingFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanCortToCort = 3/1000; %in ms
sConnParams.dblDelaySDCortToCort = 1/1000; %in ms

%connection probability excitatory cells
sConnParams.vecConnProbSD = ang2rad([7.5 360]); %degrees difference in pref ori, for [pyramid interneuron]

%V2 params
sConnParams.dblSpatialDropoffV1V2 = 0.8; %normpdf(vecX,0,0.8); zandvakili&kohn, 2015
sConnParams.dblSpatialDropoffInterneuronsV2 = 3; %for interneurons
sConnParams.intCellsV2 = sConnParams.intCellsV1;%round(sConnParams.intCellsV1/2);

%create cell-based parameters
sConnParams.vecDefinitionV2CellTypes = [1 1 1 1 2];
sConnParams.vecCellTypesV2 = repmat(sConnParams.vecDefinitionV2CellTypes,[1 ceil(sConnParams.intCellsV2/numel(sConnParams.vecDefinitionV2CellTypes))]);
sConnParams.vecCellTypesV2(sConnParams.intCellsV2+1:end) = [];
sConnParams.vecCellFractionsV2 = [0.8 0.2]; %[pyramid interneuron]
	
%V1=>V2
dblInterArealFactor = 0.8; %used to be 2.6 (2.8)
sConnParams.vecConnsPerTypeV1V2 = [48 32];%[48 32]; %from pyr to pyr/int
sConnParams.matConductancesFromToV1V2(1,:) = [1.1 1.6]*dblInterArealFactor; %from pyramid to [pyr inter]
sConnParams.matConductancesFromToV1V2(2,:) = [1.5 1.0]*dblInterArealFactor; %from inter to [pyr inter]

%synaptic delays
sConnParams.dblDelayMeanV1ToV2 = 4/1000; %in ms
sConnParams.dblDelaySDV1ToV2 = 1/1000; %in ms

%% build connectivity
sConnectivity = buildConnectivity_SimpleV1V2(sConnParams);

%% save
strConnFile = sprintf('sConnSimple_Col%dN%dS%d_%s.mat',intColumns,sConnectivity.intCortexCells,numel(sConnectivity.vecSynExcInh),getDate);
strConnDir = 'D:\Simulations\Connectivity\';

fprintf('Saving file [%s] to [%s]... [%s]\n',strConnFile,strConnDir,getTime);
save([strConnDir strConnFile],'sConnParams','sConnectivity');


