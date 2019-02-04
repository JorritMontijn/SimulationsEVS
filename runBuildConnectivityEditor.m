clear all;
strFile = 'sConn_Col48N1800S530880_2017-11-13.mat';
[sConnParams,sConnectivity] = loadConnectivity_xArea('D:\Simulations\Connectivity\',strFile);

intRemType = 0; %1, rem exc; 2, rem inh; 3, rem both
if intRemType == 0 %rem none
	strTag = 'FullConn';
else
	if intRemType == 1 %rem exc; inh only
		sConnParams.matConnCortFromTo = [0 0;120 120];
		indRem=(sConnectivity.vecSynType==1) & sConnectivity.vecSynExcInh == 1;
		strTag = 'InhOnly';
	elseif intRemType == 2 %rem inh; exc only
		sConnParams.matConnCortFromTo = [160 160; 0 0];
		indRem=(sConnectivity.vecSynType==1) & sConnectivity.vecSynExcInh == 2;
		strTag = 'ExcOnly';
	elseif intRemType == 3 %rem both; no recur
		sConnParams.matConnCortFromTo = [0 0;0 0];
		indRem=(sConnectivity.vecSynType==1);
		strTag = 'NoRecur';
	end
	
	
	sConnectivity.matSynFromTo(indRem,:)=[];
	sConnectivity.vecSynExcInh(indRem)=[];
	sConnectivity.vecSynDelay(indRem)=[];
	sConnectivity.vecSynWeight(indRem)=[];
	sConnectivity.vecSynConductance(indRem)=[];
	sConnectivity.vecSynType(indRem)=[];
end

strConnFile = ['sConn_' strTag 'Col48N2160S' num2str(numel(sConnectivity.vecSynExcInh)) '_20' getFlankedBy(strFile,'_20','.mat') '.mat'];
strConnDir = 'D:\Simulations\Connectivity\';
fprintf('Saving file [%s] to [%s]... [%s]\n',strConnFile,strConnDir,getTime);
save([strConnDir strConnFile],'sConnParams','sConnectivity');
