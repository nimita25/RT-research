function [D95, Dmax, CI, Dmean_bladder, D10_bladder, Dmean_rectum, D10_rectum, Dmean_femhead,...
    Dmean_penilebulb, Dmean_body] = calcpara_51477044(d, px, ctv, bladder, rectum, femhead, penilebulb, body)
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

% Flatten dose array to ensure it is a vector.
d = d(:);

% Calculate maximum dose and D95 for CTV.
d2 = d(ctv);
Dmax = max(d2) / px * 100;
sorted_doses_ctv = sort(d2, 'descend');
D95 = sorted_doses_ctv(round(0.95 * numel(ctv)));

% Calculate the Conformity Index (CI) for the CTV.
CI = sum(d2 >= px)^2 / (sum(d >= px) * numel(ctv));

% Calculate mean and D10 dose for the rectum.
d2 = d(rectum);
Dmean_rectum = mean(d2) / px * 100;
sorted_doses_rectum = sort(d2, 'descend');
D10_rectum = sorted_doses_rectum(round(0.1 * numel(sorted_doses_rectum))) * 100 / px;

% Calculate mean and D10 dose for the bladder.
d2 = d(bladder);
Dmean_bladder = mean(d2) / px * 100;
sorted_doses_bladder = sort(d2, 'descend');
D10_bladder = sorted_doses_bladder(round(0.1 * numel(sorted_doses_bladder))) * 100 / px;

% Calculate mean dose for the femoral heads.
d2 = d(femhead);
Dmean_femhead = mean(d2) / px * 100;

% Calculate mean dose for the penile bulb.
d2 = d(penilebulb);
Dmean_penilebulb = mean(d2) / px * 100;

% Calculate mean dose for the entire body.
d2 = d(body);
Dmean_body = mean(d2) / px * 100;

end
