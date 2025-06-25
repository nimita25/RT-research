function [D95, Dmax, CI, BEDmeanbladder, BED50_bladder, BED20_bladder, BED50_rectum, BED20_rectum, BED10_rectum,...
    BED10_femhead, BED50_penilebulb] = calcparaBED_51477044(MD, BED2, px, c)
% Calculate dosimetric parameters for radiation therapy planning.
%
% Inputs:
%   d - Dose distribution in a linear array.
%   px - Prescription dose.
%   ctv, bladder, rectum, femhead, penilebulb, body - Logical arrays indicating the region of interest.
%
% Outputs:
%   D95, Dmax - Dose covering 95% of the CTV and maximum dose in CTV.
%   CI - Conformity Index for the CTV.
%   Dmean_bladder, D10_bladder - Mean and D10 dose for the bladder.
%   Dmean_rectum, D10_rectum - Mean and D10 dose for the rectum.
%   Dmean_femhead - Mean dose for the femoral head.
%   Dmean_penilebulb - Mean dose for the penile bulb.
%   Dmean_body - Mean dose for the whole body.

% Get row indices of target, body, OARs
ctv = c{1};
bladder = c{3};
rectum = c{4};
femhead = c{5};
penilebulb = c{7};



% Calculate maximum dose and D95 for CTV.
MD_T = MD(ctv);
Dmax = max(MD_T) / px * 100;
sorted_doses_ctv = sort(MD_T, 'descend');
D95 = sorted_doses_ctv(round(0.95 * numel(ctv)));

% Calculate the Conformity Index (CI) for the CTV.
CI = sum(MD_T >= px)^2 / (sum(MD >= px) * numel(ctv));

% Calculate BED constraint values for bladder.
BED3 = BED2(bladder);
BEDmeanbladder = mean(BED3);
sorted_bladder = sort(BED3, 'descend');
BED50_bladder = sorted_bladder(round(0.5*numel(bladder)));
BED20_bladder = sorted_bladder(round(0.2*numel(bladder)));

% Calculate BED constraint values for rectum.
BED3 = BED2(rectum);
sorted_rectum = sort(BED3, 'descend');
BED50_rectum = sorted_rectum(round(0.5*numel(rectum)));
BED20_rectum = sorted_rectum(round(0.2*numel(rectum)));
BED10_rectum = sorted_rectum(round(0.1*numel(rectum)));

% Calculate BED constraint values for femhead.
BED3 = BED2(femhead);
sorted_femhead = sort(BED3, 'descend');
BED10_femhead = sorted_femhead(round(0.1*numel(femhead)));

% Calculate BED constraint values for penilebulb.
BED3 = BED2(penilebulb);
sorted_penilebulb = sort(BED3, 'descend');
BED50_penilebulb = sorted_penilebulb(round(0.5*numel(penilebulb)));

end
