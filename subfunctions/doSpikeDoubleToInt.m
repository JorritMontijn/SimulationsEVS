function cellSpikeTimes= doSpikeDoubleToInt(cellSpikeTimes,dblStepSize,dblFirstStamp)
	%UNTITLED Summary of this function goes here
	%   Detailed explanation goes here
	if iscell(cellSpikeTimes)
		for intCell=1:numel(cellSpikeTimes)
			if ~isa(cellSpikeTimes{intCell},'double'),error([mfilename ':NotADouble'],'Supplied cell array contents are not doubles!');end
			cellSpikeTimes{intCell} = int32(ceil(((cellSpikeTimes{intCell}-dblFirstStamp)/dblStepSize)+1));
		end
	else
		if ~isa(cellSpikeTimes,'double'),error([mfilename ':NotADouble'],'Supplied array is not a double!');end
		cellSpikeTimes = int32(ceil(((cellSpikeTimes-dblFirstStamp)/dblStepSize)+1));
	end
end

