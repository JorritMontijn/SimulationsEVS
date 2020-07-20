strPath = 'F:\Code\Simulations\SimulationsEVS\SimResults\';

sDir = dir([strPath 'Simulation*.mat']);

for intFile=1:numel(sDir)
	strFile = sDir(intFile).name;
	fprintf('Processing %s (%d/%d) [%s]\n',strFile,intFile,numel(sDir),getTime);
	load([strPath strFile]);
	sSimRun = rmfield(sSimRun,{'cellSpikeTimesLGN_ON','cellSpikeTimesLGN_OFF'});
	save([strPath strFile],'sSimRun','sParams','-v7');
end
