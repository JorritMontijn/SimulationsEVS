
vecNeurons = [5 432 887];
matResp3 = matResp(vecNeurons,:);
matStimResp3 = doMatRespTransform(matResp3,cellSelect)
matMeanStimResp3 = squeeze(mean(matStimResp3,2));
matSDStimResp3 = squeeze(std(matStimResp3,[],2));

vecSD = mean(matSDStimResp3,2);

[handle,matX_final,matY_final,matZ_final,matC] = plotTube3D(matMeanStimResp3(:,1),matMeanStimResp3(:,2),matMeanStimResp3(:,3),vecSD/10,1:12)


tubeplot(matMeanStimResp3',vecSD/10)
