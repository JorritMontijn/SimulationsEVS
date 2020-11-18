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

%%
return
%% type 2
strPath = 'F:\Data\Results\SimResults\';
sDir = dir([strPath '*.mat']);

for intFile=1:numel(sDir)
	strFile = sDir(intFile).name;
	fprintf('Processing %s (%d/%d) [%s]\n',strFile,intFile,numel(sDir),getTime);
	sLoad = load([strPath strFile]);

	%remove ON/OFF params
	sLoad.matSynConnOFF_to_Cort = [];
	sLoad.matSynConnON_to_Cort = [];
	sLoad.vecSynConductanceOFF_to_Cort = [];
	sLoad.vecSynConductanceON_to_Cort = [];
	sLoad.vecSynDelayOFF_to_Cort = [];
	sLoad.vecSynDelayON_to_Cort = [];
	sLoad.vecSynWeightOFF_to_Cort = [];
	sLoad.vecSynWeightON_to_Cort = [];
	
	save([strPath 'RedFile' strFile],'-struct','sLoad','-v7.3');
	
end
	%{
	%}