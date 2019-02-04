function cellSpikeTimes=doSpikeIntToDouble(cellSpikeTimes,dblStepSize,dblFirstStamp)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	
	for intCell=1:numel(cellSpikeTimes)
		cellSpikeTimes{intCell} = ((double(cellSpikeTimes{intCell})-1)*dblStepSize)+dblFirstStamp;
	end
end

