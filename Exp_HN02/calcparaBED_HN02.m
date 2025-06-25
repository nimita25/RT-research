function [D95, Dmax, CI, BEDmax_P, BEDmean_P, BEDmax_OC, BEDmean_OC,...
    BEDmax_OP, BEDmean_OP, BEDmax_L, BEDmean_L] = calcparaBED_HN02(MD, BED2, px, c)
% Calculates various dose parameters for DVH analysis in radiation therapy.
%
% INPUT:
%   d      - Dose vector (typically in Gy).
%   px     - Scaling factor for prescription dose.
%   ctv    - Logical or index array for Clinical Target Volume.
%   oar1 to oar5 - Logical or index arrays for different Organs At Risk.
%   body   - Logical or index array for the whole body.
%
% OUTPUT:
%   D95        - Dose covering 95% of the CTV.
%   Dmax       - Maximum dose received by the CTV.
%   CI         - Conformity index for the CTV.
%   Dmax_oar1  - Maximum dose received by OAR 1.
%   V10_oar1   - Volume of OAR 1 receiving at least 10 Gy.
%   Dmax_oar2  - Maximum dose received by OAR 2.
%   Dmax_oar3  - Maximum dose received by OAR 3.
%   Dmax_oar4  - Maximum dose received by OAR 4.
%   V12_oar5   - Volume of OAR 5 receiving at least 12 Gy.
%   Dmean_body - Mean dose received by the body.

% Get row indices for target, body, OARs
ctv = c{1};
oar1 = c{3};
oar2 = c{4};
oar3 = c{5};
oar5 = c{6};

% Calculate maximum dose and D95 for CTV.
MD_T = MD(ctv);
Dmax = max(MD_T) / px * 100;
sorted_doses_ctv = sort(MD_T, 'descend');
D95 = sorted_doses_ctv(round(0.95 * numel(ctv)));

% Calculate the Conformity Index (CI) for the CTV.
CI = sum(MD_T >= px)^2 / (sum(MD >= px) * numel(ctv));

% Calculate BED constraint values for parotid.
BED3 = BED2(oar1);
BEDmax_P = max(BED3);
BEDmean_P = mean(BED3);

% Calculate BED constraint values for oral cavity.
BED3 = BED2(oar2);
BEDmax_OC = max(BED3);
BEDmean_OC = mean(BED3);

% Calculate BED constraint values for Oropharynx.
BED3 = BED2(oar3);
BEDmax_OP = max(BED3);
BEDmean_OP = mean(BED3);

% Calculate BED constraint values for Larynx.
BED3 = BED2(oar5);
BEDmax_L = max(BED3);
BEDmean_L = mean(BED3);

end
