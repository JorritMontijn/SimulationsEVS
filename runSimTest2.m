clear all;clc;%close all;clc;
strConnNew = 'Conn256N1200_2021-03-24'; %new
dblNoise = 0;
strStimNew = sprintf('Ret256Noise%.1fOri5_x2R1_2020-07-17.mat',dblNoise); %new


strInput = ['time=1,conn=' strConnNew...
	',stim=' strStimNew...
	',idx=0,tag=Ori2Noise' num2str(dblNoise)];
[sData,sSimRun]=runSimulation(strInput);

%{
Elapsed: 5.7s; now at t=0.095s; mean rate (Hz): 3.064 (V1 Pyr); 9.744 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:13:53]
Elapsed: 10.7s; now at t=0.144s; mean rate (Hz): 4.087 (V1 Pyr); 12.456 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:13:58]
Elapsed: 15.7s; now at t=0.197s; mean rate (Hz): 4.394 (V1 Pyr); 13.409 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:03]
Elapsed: 20.8s; now at t=0.247s; mean rate (Hz): 4.706 (V1 Pyr); 14.018 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:08]
Elapsed: 25.8s; now at t=0.301s; mean rate (Hz): 4.840 (V1 Pyr); 14.442 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:13]
Elapsed: 30.8s; now at t=0.355s; mean rate (Hz): 5.651 (V1 Pyr); 16.127 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:18]
Elapsed: 35.8s; now at t=0.405s; mean rate (Hz): 5.494 (V1 Pyr); 15.628 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:23]
Elapsed: 40.8s; now at t=0.455s; mean rate (Hz): 5.929 (V1 Pyr); 16.703 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:28]
Elapsed: 45.9s; now at t=0.504s; mean rate (Hz): 5.911 (V1 Pyr); 16.700 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:33]
Elapsed: 50.9s; now at t=0.554s; mean rate (Hz): 5.745 (V1 Pyr); 16.419 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:38]
Elapsed: 55.9s; now at t=0.605s; mean rate (Hz): 5.990 (V1 Pyr); 17.149 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:43]
Elapsed: 60.9s; now at t=0.658s; mean rate (Hz): 5.682 (V1 Pyr); 16.527 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:48]
Elapsed: 65.9s; now at t=0.710s; mean rate (Hz): 5.630 (V1 Pyr); 16.450 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:53]
Elapsed: 70.9s; now at t=0.764s; mean rate (Hz): 5.498 (V1 Pyr); 16.165 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:14:58]
Elapsed: 76.0s; now at t=0.821s; mean rate (Hz): 5.195 (V1 Pyr); 15.429 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:15:03]
Elapsed: 81.0s; now at t=0.881s; mean rate (Hz): 5.258 (V1 Pyr); 15.645 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:15:08]
Elapsed: 86.0s; now at t=0.938s; mean rate (Hz): 5.331 (V1 Pyr); 15.832 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:15:13]
Elapsed: 91.0s; now at t=0.993s; mean rate (Hz): 5.404 (V1 Pyr); 15.958 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:15:18]
Elapsed: 96.0s; now at t=1.046s; mean rate (Hz): 5.381 (V1 Pyr); 15.886 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:15:23]


%0.94
Elapsed: 1488.6s; now at t=0.301s; mean rate (Hz): 2.536 (V1 Pyr); 0.650 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [11:51:32]

%Conn256N1200_2020-10-07.mat:
%Elapsed: 1135.7s; now at t=0.202s; mean rate (Hz): 1.639 (V1 Pyr); 0.889 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [14:01:13]

%Conn256N1200_2020-10-26.mat:
%Elapsed: 1094.4s; now at t=0.397s; mean rate (Hz): 6.533 (V1 Pyr); 4.828 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [12:28:46]

%Conn256N1200_2020-10-28.mat:
%Elapsed: 1096.1s; now at t=0.398s; mean rate (Hz): 5.883 (V1 Pyr); 5.115 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [12:10:01]

%Conn256N1200_2020-10-29.mat:
%Elapsed: 1112.5s; now at t=0.399s; mean rate (Hz): 5.725 (V1 Pyr); 4.887 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [11:45:44]

strConnNew = 'Conn256N1200_2021-03-22'; %new
Elapsed: 1458.4s; now at t=0.663s; mean rate (Hz): 7.040 (V1 Pyr); 7.296 (V1 Int); NaN (V2 Pyr); NaN (V2 Int) [19:40:02]


%}