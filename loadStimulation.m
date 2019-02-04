function [sStimParams,sStimInputs] = loadStimulation(strStimDir,strStimFile)
	%UNTITLED4 Summary of this function goes here
	%   Detailed explanation goes here
	
	%% load data
	sLoad = load([strStimDir strStimFile]);
	sStimInputs = sLoad.sStimInputs;
	sStimParams = sLoad.sStimParams;
end

