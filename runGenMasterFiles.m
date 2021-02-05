vecNoise = 0:0.2:3;
strOldPath = cd('F:\Data\Results\SimResults\todo\');

%ori tuning
strArgs = ['time=0,conn=Conn256N1200_2021-01-08.mat,stim=Ret256Noise0.0Ori160_x9R1_2020-11-02.mat,idx=0,tag=OriTuningNoise0'];
runSimulation(strArgs);
%noise
for dblNoise=vecNoise
	strArgs = ['time=0,conn=Conn256N1200_2021-01-08.mat,stim=Ret256Noise' sprintf('%1.1f',dblNoise) 'Ori5_x2R1_2020-07-17.mat,idx=0,tag=Ori5Noise' sprintf('%02d',round(dblNoise*10))];
	runSimulation(strArgs);
end

cd(strOldPath);