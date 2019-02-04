
strDir = 'D:\Data\Results\NonLeaky\';
sDir = dir([strDir '*_Type4?_*']);

intPopSizes = numel(sDir);
vecSizePops = nan(1,intPopSizes);
vecCorrPops = nan(1,intPopSizes);
for intPopSizeIdx=1:intPopSizes
	strFile = sDir(intPopSizeIdx).name;
	sLoad = load([strDir strFile]);
	
	matCorr = corr(sLoad.matModelResp');
	vecCorr = matCorr(tril(true(size(matCorr)),-1));
	
	vecCorrPops(intPopSizeIdx) = mean(vecCorr);
	vecSizePops(intPopSizeIdx) = size(matCorr,1);
end

figure
plot(vecSizePops,vecCorrPops,'o-');
xlim([0 max(get(gca,'xlim'))]);
ylim([0 max(get(gca,'ylim'))]);
