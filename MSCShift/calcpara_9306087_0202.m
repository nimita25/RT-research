function [D95, Dmax, CI, Dmax_oar1, V10_oar1, Dmax_oar2, Dmax_oar3, Dmax_oar4, V12_oar5, Dmean_body] = calcpara_9306087_0202(nfrac, d, px, c)
% Calculates various dose parameters for DVH analysis in radiation therapy.
%
% INPUT:
%   nfrac  - Number of fractions of therapy.
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
oar4 = c{6};
oar5 = c{7};
body = c{2};

% Ensure the dose vector is in the correct format
d = d(:);

% CTV-related calculations
d_ctv = d(ctv);
Dmax = max(d_ctv) / px * 100;
sorted_ctv_doses = sort(d_ctv, 'descend');
D95 = sorted_ctv_doses(round(0.95 * numel(d_ctv)));

% Conformity index calculation
CI = (sum(d_ctv >= px)^2) / (sum(d >= px) * numel(ctv));

% OAR calculations
Dmax_oar1 = max(d(oar1));
V10_oar1 = sum(d(oar1) >= 15 / nfrac) / numel(oar1) * 100; % Percent volume receiving at least 15 Gy

Dmax_oar2 = max(d(oar2));
Dmax_oar3 = max(d(oar3));
Dmax_oar4 = max(d(oar4));

V12_oar5 = sum(d(oar5) >= 12 / nfrac) / numel(oar5) * 100; % Percent volume receiving at least 12 Gy

% Whole body mean dose calculation
Dmean_body = mean(d(body));

end
