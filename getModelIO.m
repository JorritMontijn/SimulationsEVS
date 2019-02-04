function vecRespHz = getModelIO(vecParams,vecPrefOri)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	%% get connectivity
	global sConnIO;
	
	%% get params
	dblInputG = vecParams(1)*1000;
	dblGaussSD = vecParams(2)*1000;
	dblOffset = vecParams(3)*100;
	dblFracInputToInh = vecParams(4);
	
	%% build input
	%set params
	dblDeltaT = 0.0005;
	dblDur = 0.5; %seconds
	dblSynSpikeMem = 0.2;
	
	%get input
	intNeurons = sConnIO.intCortexCells;
	vecThisV = gaussrnd(-56,1,[intNeurons,1]);
	vecGauss = normpdf(1:intNeurons,intNeurons/2,dblGaussSD)';
	vecInput = dblOffset + dblInputG*(vecGauss/max(vecGauss));
	vecInputIdx = ones(round(dblDur/dblDeltaT),1) * (1:numel(dblInputG));
	vecInputIdx = vecInputIdx(:)';
	vecOverallT = (1:numel(vecInputIdx))*dblDeltaT;
	indInh = (sConnIO.vecCellTypes==2);
	vecInput(indInh) = vecInput(indInh)*dblFracInputToInh;
	
	
	%put in sData
	sData = sConnIO;
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
	
	%% get output
	cellSpikeTimesCortex = sData.cellSpikeTimesCortex;
	vecRespHz = (cellfun(@numel,cellSpikeTimesCortex) - cellfun(@sum,cellfun(@isnan,cellSpikeTimesCortex,'UniformOutput',false)))/dblDur;
	
	%%msg
	fprintf('G-input: %f; SD: %f; Offset: %f; Ratio I/P: %f [%s]\n',dblInputG,dblGaussSD,dblOffset,dblFracInputToInh,getTime);
	
end

