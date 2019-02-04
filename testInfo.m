intDists = 5;
vecDists = 2*(1:intDists);
intIters=10;
matDprimes = nan(intDists,intIters,4);

for intDist=1:intDists
	intDist
	dblDist = vecDists(intDist);
	for intIter=1:intIters
		intNeurons = 100;
		intTrials = 100;
		
		vecMu1 = zeros(1,intNeurons);
		matCov1 = eye(intNeurons);
		matResp1 = mvnrnd(vecMu1,matCov1,intTrials)';
		
		vecMu2 = ones(1,intNeurons);
		vecMu2 = vecMu2/norm(vecMu2);
		vecMu2 = dblDist*vecMu2;
		matCov2 = eye(intNeurons);
		matResp2 = mvnrnd(vecMu2,matCov2,intTrials)';
		
		boolLinear = 0;
		matData = cat(2,matResp1,matResp2);
		vecTrialTypes = [ones(1,intTrials) 2*ones(1,intTrials)];
		
		%calc direct
		[dblPredA,matPredA,dblDprimeSquaredOffDiagonal,matDprimeSquared,matDprimeSquared_diagonal] = ...
			getSeparation(matData,vecTrialTypes,boolLinear);
		
		dblSubFac =(2*intNeurons)/(intTrials);
		dblProdFacRaw = ((2*intTrials-intNeurons-3)/(2*intTrials-2));
		dblInputI = dblDprimeSquaredOffDiagonal*dblProdFacRaw-dblSubFac;
		
		matDprimes(intDist,intIter,1) = sqrt(dblInputI);
		
		%calc log reg
		%set parameters
		sParamsAnal = struct;
		sParamsAnal.vecUseStimTypes = [1 2];
		sParamsAnal.dblLambda = 1;
		sParamsAnal.boolDirectI = true;
		sParamsAnal.intIters = 3;
		sParamsAnal.vecGroupSizes = intNeurons;
		sParamsAnal.boolVerbose = false;
		sParamsAnal.dblDiffTheta = 1;
		sParamsAnal.boolBiasCorrection = 1;
		
		%do analysis
		[matFisherFull,sOut] = doFisherAnal(matData,vecTrialTypes,sParamsAnal);
		
		dblI_CV = mean(sOut.matI_LogReg_bc_CV(:));
		dblI = mean(sOut.matI_LogReg_bc(:));
		dblI_Dir = mean(sOut.matI_Direct_bc(:));
		
		matDprimes(intDist,intIter,2) = sqrt(dblI);
		matDprimes(intDist,intIter,3) = sqrt(dblI_CV);
		matDprimes(intDist,intIter,4) = sqrt(dblI_Dir);
		
	end
end

vecMeanDirect1 = xmean(matDprimes(:,:,1),2);
vecSDDirect1 = xstd(matDprimes(:,:,1),2);

vecMeanLR = xmean(matDprimes(:,:,2),2);
vecSDLR = xstd(matDprimes(:,:,2),2);

vecMeanLRCV = xmean(matDprimes(:,:,3),2);
vecSDLRCV = xstd(matDprimes(:,:,3),2);

vecMeanDirect2 = xmean(matDprimes(:,:,4),2);
vecSDDirect2 = xstd(matDprimes(:,:,4),2);
%%
figure
hold on
errorbar(vecDists,vecMeanDirect1,vecSDDirect1/sqrt(intIters),'b--');
errorbar(vecDists,vecMeanLR,vecSDLR/sqrt(intIters),'r-');
errorbar(vecDists,vecMeanLRCV,vecSDLRCV/sqrt(intIters),'r--');
errorbar(vecDists,vecMeanDirect2,vecSDDirect2/sqrt(intIters),'b-');
legend({'Direct1','LR','LR CV','Direct2'},'Location','Best');
grid on

title(sprintf('Input vs extracted Fisher I; neurons=%d;trials=%d',intNeurons,intTrials));
xlabel('Input distance between S1 and S2 (d'')');
ylabel('Extracted d''');
fixfig

