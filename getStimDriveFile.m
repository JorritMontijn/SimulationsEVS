function strFile = getStimDriveFile(sParams)
	%UNTITLED2 Summary of this function goes here
	%   Detailed explanation goes here
	
	strFile = 'StimDrive';
	if isfield(sParams,'ST')
		strFile = [strFile sprintf('_ST%d_',sParams.ST)];
	end
	if isfield(sParams,'SF')
		strFile = [strFile sprintf('_SF%.3f_',sParams.SF)];
	end
	if isfield(sParams,'TF')
		strFile = [strFile sprintf('_TF%.3f_',sParams.TF)];
	end
	if isfield(sParams,'Ori')
		strFile = [strFile sprintf('_Ori%.3f_',sParams.Ori)];
	end	
	
	if isfield(sParams,'FR')
		strFile = [strFile sprintf('_FR%3d_',sParams.FR)];
	end
	
	if isfield(sParams,'SSRD')
		strFile = [strFile sprintf('_SSRD%.1f_',sParams.SSRD(1))];
	end
	if isfield(sParams,'PixWH')
		strFile = [strFile sprintf('_PixWH%d_',sParams.PixWH(1))];
	end
	if isfield(sParams,'DegWH')
		strFile = [strFile sprintf('_DegWH%.1f_',sParams.DegWH(1))];
	end	
	
	if isfield(sParams,'dT')
		strFile = [strFile sprintf('_dT%.4f_',sParams.dT)];
	end
	if isfield(sParams,'SD')
		strFile = [strFile sprintf('_SD%.3f_',sParams.SD)];
	end
	if isfield(sParams,'BD')
		strFile = [strFile sprintf('_BD%.3f_',sParams.BD)];
	end	
	strFile = [strFile '.mat'];
end

