function [sConnParams,sConnectivity] = loadConnectivity_xArea(strConnDir,strConnFile)
	%UNTITLED4 Summary of this function goes here
	%   Detailed explanation goes here
	
	%% load data
	sLoad = load([strConnDir strConnFile]);
	sConnParams = sLoad.sConnParams;
	sConnectivity = sLoad.sConnectivity;
end

