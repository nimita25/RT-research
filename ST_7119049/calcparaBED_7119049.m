function [D95, Dmax, CI, BEDmean_lung, BED30_lung, BEDmean_heart, BED60_heart, ...
    BEDmean_eso] =  calcparaBED_7119049(MD, BED2, px, c)
    
    
% Calculates radiation therapy parameters for specified anatomical regions.
%
% INPUTS:
%   d     - Array of dose values (can be in 2D or 3D matrix form, must be converted to vector).
%   px    - Prescription dose, used for scaling doses.
%   ctv   - Logical index or mask for Clinical Target Volume.
%   lung  - Logical index or mask for lungs.
%   heart - Logical index or mask for heart.
%   esoph - Logical index or mask for esophagus.
%   cord  - Logical index or mask for spinal cord.
%   body  - Logical index or mask for whole body.
%
% OUTPUTS:
%   D95         - Dose covering 95% of the CTV.
%   Dmax        - Maximum dose in the CTV.
%   CI          - Conformity index based on the prescription dose.
%   Dmean_lung  - Mean dose to the lungs.
%   V20_lung    - Percentage of lung volume receiving at least 20% of the prescription dose.
%   Dmean_heart - Mean dose to the heart.
%   V30_heart   - Percentage of heart volume receiving at least 30% of the prescription dose.
%   Dmean_esoph - Mean dose to the esophagus.
%   Dmean_cord  - Mean dose to the spinal cord.
%   Dmean_body  - Mean dose to the whole body.

% Get row indices of target, body, OARs
ctv = c{1};
lung= c{3};
heart= c{4};
esoph= c{5};

% Calculate maximum dose and D95 for CTV.
MD_T = MD(ctv);
sorted_doses_ctv = sort(MD_T, 'descend');
tmp_NF = sorted_doses_ctv(round(0.95 * numel(ctv)))/px;
Dmax = (max(MD_T) / px * 100)/tmp_NF;
D95 = (sorted_doses_ctv(round(0.95 * numel(ctv))))/tmp_NF;

% Calculate the Conformity Index (CI) for the CTV.
CI = sum(MD_T >= px)^2 / (sum(MD >= px) * numel(ctv));

% Calculate BED constraint values for lung.
BED3 = BED2(lung);
sorted_lung = sort(BED3, 'descend');
BED30_lung = sorted_lung(round(0.3*numel(lung)));
BEDmean_lung = mean(BED3);

% Calculate BED constraint values for heart.
BED3 = BED2(heart);
sorted_heart = sort(BED3, 'descend');
BED60_heart = sorted_heart(round(0.3*numel(heart)));
BEDmean_heart = mean(BED3);

% Calculate BED constraint values for esophagus.
BED3 = BED2(esoph);
BEDmean_eso = mean(BED3);

end
