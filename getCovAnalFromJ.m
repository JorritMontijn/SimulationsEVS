function matCov = getCovAnalFromJ(matJ, dblVarIndep, dblVarShared)
% matJ should take form of to-from

if nargin < 2
	dblVarIndep = 76.5;
	dblVarShared = 3.5;
end

intN = size(matJ,1);

matI = eye(intN);
mat1 = ones(intN);
invJ = inv(matJ);

varNumerator = (dblVarIndep*matI + dblVarShared * mat1);
matCov = invJ * (varNumerator / matJ');
