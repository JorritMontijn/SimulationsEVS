function sStimDrive = loadStimDrive(sParams,strDir)
	%UNTITLED4 Summary of this function goes here
	%   Detailed explanation goes here
	
	if ~exist('strDir','var') || isempty(strDir),strDir = 'D:\Simulations\StimDrives\';end
	
	%get all file names
	sDir = dir(strDir);
	cellFiles = struct2cell(sDir);
	cellFiles = cellFiles(strcmp(fieldnames(sDir),'name'),:);
	
	%get files that satisfy criteria
	indSatisFiles = true(size(cellFiles));
	
	%get criteria
	cellParams = fieldnames(sParams);
	for intParam=1:numel(cellParams)
		strParam = cellParams{intParam};
		varVal = sParams.(strParam);
		
		%loop through files
		for intFile=find(indSatisFiles==true)
			strFile = cellFiles{intFile};
			
			strFileVal=getFlankedBy(strFile,['_' strParam],'_');
			if isempty(strFileVal)
				indSatisFiles(intFile) = false;
			else
				cellVals = strsplit(strFileVal,'.');
				if numel(cellVals) == 1
					intPrecision = 0;
				else
					intPrecision = length(cellVals{2});
				end
				if roundi(varVal,intPrecision) ~= roundi(str2double(strFileVal),intPrecision)
					indSatisFiles(intFile) = false;
				end
			end
		end
	end
	
	%check if any files match criteria
	if sum(indSatisFiles) == 0
		fprintf('loadStimDrive:No files were found that match requested criteria [%s]\n',getTime);
		sStimDrive = [];
	else
		if sum(indSatisFiles) > 1
			warning([mfilename ':AmbiguousRequest'],sprintf('Multiple files (%d) matched requested criteria',sum(indSatisFiles)));
		end
		strLoadFile = cellFiles{find(indSatisFiles==true,1)};
		fprintf('loadStimDrive:loading <%s> [%s]\n',strLoadFile,getTime);
		sStimDrive = load([strDir strLoadFile]);
	end
end

