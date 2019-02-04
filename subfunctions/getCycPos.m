function intStartPos = getCycPos(intNewPos,intStoredVals,intCylinderSize)
	intStartPos = mod(intNewPos-intStoredVals,intCylinderSize);
	if intStartPos == 0,intStartPos = intCylinderSize;end
end
