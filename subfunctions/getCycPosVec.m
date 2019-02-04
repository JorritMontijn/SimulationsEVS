function vecPos = getCycPosVec(intNewPos,intStoredVals,intCylinderSize)
	if intNewPos <= intStoredVals
		intNewPos = intNewPos + intCylinderSize;
	end
	vecPos = mod((intNewPos-intStoredVals):(intNewPos-1),intCylinderSize);
	vecPos(vecPos==0)=intCylinderSize;
end

