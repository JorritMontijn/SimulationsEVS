function matJ = getBandJ(intPyrWidth,intInhWidth)

intN = 1440;
inh_period = 5;
if ~exist('intPyrWidth','var'),intPyrWidth = 10;end
if ~exist('intInhWidth','var'),intInhWidth = 160;end

dblConnEE = 6/(intPyrWidth-1); %to E from E; 6
dblConnIE = 5.9/(intPyrWidth); %to I from E; 5.9
dblConnEI = -9.5/(intInhWidth); %to E from I; -9.5
dblConnII = -9.4/(intInhWidth-1); %to I from I; -9.4
	
idx_pyr = mod([1:intN],inh_period) ~= 0;
idx_inh = mod([1:intN],inh_period) == 0;

intPyr = sum(idx_pyr);
intInh = sum(idx_inh);

pyr_start = -floor(intPyrWidth/2);
pyr_range = pyr_start:pyr_start+intPyrWidth-1;
inh_start = -floor(intInhWidth/2);
inh_range = inh_start:inh_start+intInhWidth-1;

pyr_pyr_filter = ...
	bsxfun(@plus,...
		mod(bsxfun(@plus,pyr_range,[0:intPyr-1]'),intPyr),...
		(1+intN*[0:intPyr-1]')...
	);

pyr_inh_filter = ...
	bsxfun(@plus,...
		mod(bsxfun(@plus,pyr_range,floor(intPyr*1.0/intInh*[0:intInh-1]')),intPyr),...
		(1+intN*intPyr+intN*[0:intInh-1]')...
	);

inh_inh_filter = ...
	bsxfun(@plus,...
		mod(bsxfun(@plus,inh_range,[0:intInh-1]'),intInh),...
		(1+(intN+1)*intPyr+intN*[0:intInh-1]')...
	);
	
inh_pyr_filter = ...
	bsxfun(@plus,...
		mod(bsxfun(@plus,inh_range,floor(intInh*1.0/intPyr*[0:intPyr-1]')),intInh),...
		(1+intPyr+intN*[0:intPyr-1]')...
	);

matJ = zeros(intN);
matJ(pyr_pyr_filter) = dblConnEE;
matJ(inh_inh_filter) = dblConnII;
matJ(pyr_inh_filter) = dblConnIE;
matJ(inh_pyr_filter) = dblConnEI;
matJ(logical(eye(intN))) = 0;%dblV_reset - dblV_thresh;

matJ = matJ';
