function [matStim,rampGrid] = buildGrating(dblSF,dblStimSizeRetDeg,vecScrPixWidthHeight,vecScrDegWidthHeight,dblPhaseOffset)
	%buildGrating outputs greyscale grating patch
	%   syntax: matStim = buildGrating(SpatFreq,StimSizeRetDeg,ScrPixWidthHeight,ScrDegWidthHeight,phaseOffset)
	%		- Creates a grating enclosed by circular neutral-grey window
	%		  with cosine ramps around the border. For more information on
	%		  the circular window, see buildCircularCosineRamp.m
	%
	%	output: 
	%		- matStim; matrix, range 0-255
	%
	%	inputs:
	%		- SpatFreq; double, spatial frequency (cycles per degree)
	%		- StimSizeRetDeg; double, size of stimulus in retinal degree
	%		- ScrPixWidthHeight; 2-element (integer) vector, first element
	%			is screen width; second is height in pixels
	%		- ScrDegWidthHeight; 2-element (double) vector, first element
	%			is screen width; second is height in retinal degrees
	%		- phaseOffset; double, phase offset in range 0-2pi
	%
	%	Version History:
	%	2012-03-27	Created by Jorrit Montijn
	
	
	%these are some default values that can be used to test the script:
	if ~exist('dblSF','var')
		dblSF = 0.05; %cycles per degree
		dblStimSizeRetDeg = 40;
		ScrW_pix = 1280;
		ScrH_pix = 1024;
		ScrW_deg = 93;
		ScrH_deg = 80;
	else
		ScrH_pix = vecScrPixWidthHeight(1);
		ScrW_pix = vecScrPixWidthHeight(end);
		ScrH_deg = vecScrDegWidthHeight(1);
		ScrW_deg = vecScrDegWidthHeight(end);
	end
	
	if ~exist('dblPhaseOffset','var') || isempty(dblPhaseOffset)
		dblPhaseOffset = 0;
	end
	
	pixPerDeg = ScrH_pix / ScrH_deg; %number of pixels in single retinal degree
	imSize = ScrH_pix;
	
	degreesPerCycle = 1/dblSF; %number of degrees in a single cycle (black-white block)
	pixPerCycle = degreesPerCycle * pixPerDeg; %number of pixels in such a cycle
	
	[grid]=meshgrid(1:imSize); %create a grid with the size of the required image
	
	%build the grating
	modmat = mod(grid-1+pixPerCycle*(dblPhaseOffset/(2*pi)),pixPerCycle);  %every pixPerCycle pixels the grid flips back to 0 with an offset of phaseOffset
	matGrat = (modmat >= pixPerCycle/2)*255; %create logical 1s and 0s to build the black/white grating
	
	%calculate window variables
	windowDiameter = dblStimSizeRetDeg * pixPerDeg; %size of non-ramped stimulus
	rampWidth = 0.5*pixPerDeg; %size of cosine ramp where the stimulus goes from maximal intensities to background grey
	
	%build the circular ramp using another function and also create an inverse mask
	rampGrid = buildCircularCosineRampPix(imSize,windowDiameter,rampWidth);
	rampGridInverse = abs(rampGrid - 1);
	
	stimPart = matGrat .* rampGrid; %multiply the ramped mask with the stimulus
	bgPart = 128 .* rampGridInverse; %and multiple the inverse of that mask with the background
	
	matStim = stimPart + bgPart; %add them together and we're done!
end
