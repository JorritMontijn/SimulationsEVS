%% load data
strConn = 'Conn256N1200_2020-10-26.mat'; %new
strStim = 'Ret256Noise0.0Ori5_x2R1_2020-07-17.mat'; %new
load(['F:\Code\Simulations\SimulationsEVS\Connectivity\' strConn ]);
load(['F:\Code\Simulations\SimulationsEVS\Stimulation\' strStim]);

%% plot response of example cell
intNeuron=25;
dblPrefSF = sConnectivity.vecPrefSF(intNeuron);
indSyn = sConnectivity.matSynConnOFF_to_Cort(:,2)==intNeuron;
vecSourceOff = sConnectivity.matSynConnOFF_to_Cort(indSyn,1);
vecW_LGN = sConnectivity.vecSynWeightOFF_to_Cort(indSyn);
vecC_LGN = sConnectivity.vecSynConductanceOFF_to_Cort(indSyn);

vecW_E = sConnectivity.vecSynWeight(sConnectivity.vecSynExcInh==1);
vecC_E = sConnectivity.vecSynConductance(sConnectivity.vecSynExcInh==1);

vecW_I = sConnectivity.vecSynWeight(sConnectivity.vecSynExcInh==2);
vecC_I = sConnectivity.vecSynConductance(sConnectivity.vecSynExcInh==2);

%% PSP
dblTauPeakE = sConnectivity.vecTauPeakByType(1); %E
dblTauPeakI = sConnectivity.vecTauPeakByType(2); %E

%synapse function
fPSP = @(dblT,vecSpikeTimes,dblTauPeak) sum(max(0,dblT - vecSpikeTimes).*(1/dblTauPeak).*exp(1-(dblT - vecSpikeTimes)/dblTauPeak));
fPSP2 = @(matT,dblTauPeak) sum(max(0,matT).*(1/dblTauPeak).*exp(1-matT/dblTauPeak),2);

%peak of PSP is 1, so we can ignore it
%i=0;
%for t=dblDeltaT:dblDeltaT:0.2
%	i=i+1;
%	x(i)=fPSP(t,0,dblTauPeakE);
%end
%max(x);

%% d_mV per spike
vecCellThresh = sConnectivity.vecCellThresh;
vecCellTauPeak = sConnectivity.vecCellTauPeak;
vecCellV_E = sConnectivity.vecCellV_E;
vecCellV_I = sConnectivity.vecCellV_I;

vec_d_mV_E = vecCellV_E - vecCellThresh;
dblMedian_dmV_LGN = median(vecW_LGN.*vecC_LGN)*median(vec_d_mV_E)
dblMax_dmV_LGN = max(vecW_LGN.*vecC_LGN)*median(vec_d_mV_E)

	
vec_d_mV_E = vecCellV_E - vecCellThresh;
dblMedian_dmV_E = median(vecW_E.*vecC_E)*median(vec_d_mV_E)
dblMax_dmV_E = max(vecW_E.*vecC_E)*median(vec_d_mV_E)


vec_d_mV_I = vecCellV_I - vecCellThresh;
dblMedian_dmV_I = median(vecW_I.*vecC_I)*median(vec_d_mV_I)
dblMax_dmV_I = max(vecW_I.*vecC_I)*median(vec_d_mV_I)

fprintf('LGN dmV, med=%.3f,max=%.3f; E dmV, med=%.3f,max=%.3f; E dmV, med=%.3f,max=%.4f\n',...
	dblMedian_dmV_LGN,dblMax_dmV_LGN,dblMedian_dmV_E,dblMax_dmV_E,dblMedian_dmV_I,dblMax_dmV_I);