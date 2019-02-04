function [matStim,rampGrid] = buildLine(LineWidthRetDeg,StimSizeRetDeg,ScrPixWidthHeight,ScrDegWidthHeight,phaseOffset)
	%buildLine outputs greyscale line stimulus
	%   syntax: [matStim,rampGrid] = buildLine(LineWidthRetDeg,StimSizeRetDeg,ScrPixWidthHeight,ScrDegWidthHeight,phaseOffset)
	%		- Creates a grating enclosed by circular neutral-grey window
	%		  with cosine ramps around the border. For more information on
	%		  the circular window, see buildCircularCosineRamp.m
	%
	%	output: 
	%		- matStim; matrix, range 0-255
	%
	%	inputs:
	%		- StimSizeRetDeg; double, size of stimulus in retinal degree
	%		- ScrPixWidthHeight; 2-element (integer) vector, first element
	%			is screen width; second is height in pixels
	%		- ScrDegWidthHeight; 2-element (double) vector, first element
	%			is screen width; second is height in retinal degrees
	%		- phaseOffset; double, phase offset in range 0-2pi
	%
	%	Version History:
	%	2017-01-13	Created by Jorrit Montijn
	
	%LineWidthRetDeg = 10
	%StimSizeRetDeg = 16
	%ScrPixWidthHeight = [128 128]
	%ScrDegWidthHeight = [25.6 25.6]
	%phaseOffset = 0.1257
	
	%these are some default values that can be used to test the script:
	ScrH_pix = ScrPixWidthHeight(1);
	ScrW_pix = ScrPixWidthHeight(end);
	ScrH_deg = ScrDegWidthHeight(1);
	ScrW_deg = ScrDegWidthHeight(end);
	
	if ~exist('phaseOffset','var') || isempty(phaseOffset)
		phaseOffset = 0;
	end
	
	pixPerDeg = ScrH_pix / ScrH_deg; %number of pixels in single retinal degree
	imSize = ScrH_pix;
	
	maxDisplacementDeg = LineWidthRetDeg;%number of degrees in a single cycle (black-white block)
	maxDisplacementPix = maxDisplacementDeg * pixPerDeg; %number of pixels in such a cycle
	
	[grid]=meshgrid(1:imSize); %create a grid with the size of the required image
	
	%build the line
	LineWidthPixBorder = (LineWidthRetDeg/2 * pixPerDeg);
	dblPhase = mod(phaseOffset/pi,2);
	modmat = abs(grid-imSize/2-maxDisplacementPix/2+maxDisplacementPix*dblPhase) <= LineWidthPixBorder;
	matLine = modmat*255; %create logical 1s and 0s to build the black/white grating
	
	%calculate window variables
	windowDiameter = StimSizeRetDeg * pixPerDeg; %size of non-ramped stimulus
	rampWidth = 4 * pixPerDeg; %size of cosine ramp where the stimulus goes from maximal intensities to background grey
	
	%build the circular ramp using another function and also create an inverse mask
	rampGrid = buildCircularCosineRamp(imSize,windowDiameter,rampWidth);
	rampGridInverse = abs(rampGrid - 1);
	
	stimPart = matLine .* rampGrid; %multiply the ramped mask with the stimulus
	bgPart = 128 .* rampGridInverse; %and multiple the inverse of that mask with the background
	
	matStim = stimPart + bgPart; %add them together and we're done!
end
