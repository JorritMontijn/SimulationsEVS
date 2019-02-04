function matJ = getBandJ()

intN = 1440;
inh_period = 5;

pyr_width = 160;
inh_width = 120;

dblV_reset = 0;
dblV_thresh = 5;

matConductancesFromTo(1,:) = [1.4643 2.1299]; %from pyramid to [pyr inter]
matConductancesFromTo(2,:) = [-4.0563 -2.7042]; %from inter to [pyr inter]

idx_pyr = mod([1:intN],inh_period) ~= 0;
idx_inh = mod([1:intN],inh_period) == 0;

intPyr = sum(idx_pyr);
intInh = sum(idx_inh);

pyr_start = -floor(pyr_width/2);
pyr_range = pyr_start:pyr_start+pyr_width-1;
inh_start = -floor(inh_width/2);
inh_range = inh_start:inh_start+inh_width-1;

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
matJ(pyr_pyr_filter) = matConductancesFromTo(1,1);
matJ(inh_inh_filter) = matConductancesFromTo(2,2);
matJ(pyr_inh_filter) = matConductancesFromTo(1,2);
matJ(inh_pyr_filter) = matConductancesFromTo(2,1);
matJ(logical(eye(intN))) = dblV_reset - dblV_thresh;
