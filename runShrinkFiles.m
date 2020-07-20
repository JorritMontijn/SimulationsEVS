strPath = 'F:\Data\Results\SimResults\';

sDir = dir([strPath '*.mat']);

for intFile=1:numel(sDir)
	strFile = sDir(intFile).name;
	fprintf('Processing %s (%d/%d) [%s]\n',strFile,intFile,numel(sDir),getTime);
	load([strPath strFile]);
	if isfield(sSimRun,'cellSpikeTimesLGN_ON')
	sSimRun = rmfield(sSimRun,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	end
	if isfield(sSimRun,'cellSpikeTimesLGN_ON')
	sSimRun = rmfield(sSimRun,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	end
	if isfield(sParams,'sConnectivity')
	sParams = rmfield(sParams,{'sConnectivity'});
	end
	
	save([strPath strFile],'sSimRun','sParams','-v7');
end
