vecNoise = 0:0.2:3;
strOldPath = cd('F:\Data\Results\SimResults\todo\');
for dblNoise=vecNoise
	strArgs = ['time=0,conn=Conn256N1200_2020-10-29.mat,stim=Ret256Noise' sprintf('%1.1f',dblNoise) 'Ori5_x2R1_2020-07-17.mat,idx=0,tag=Ori5Noise' sprintf('%02d',dblNoise*10)];
	runSimulation(strArgs);
end

cd(strOldPath);