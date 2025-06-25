function [D95, Dmax, CI, Dmean_lung, V20_lung, Dmean_heart, V30_heart, Dmean_esoph,...
    Dmean_body] = calcpara_7119049(d, px, ctv, lung, heart, esoph, body)
    %cord, ...
    
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

% Ensure the dose matrix is in vector form for processing
d = d(:);

% CTV-specific calculations
d_ctv = d(ctv);
Dmax = max(d_ctv) / px * 100; % Calculate maximum dose as a percentage of prescription dose
sorted_doses_ctv = sort(d_ctv, 'descend');
D95 = sorted_doses_ctv(round(0.95 * numel(d_ctv))); % Dose to 95% of CTV

% Conformity Index for the CTV
CI = (sum(d_ctv >= px)^2) / (sum(d >= px) * numel(ctv));

% Lung calculations
d_lung = d(lung);
Dmean_lung = mean(d_lung) / px * 100; % Mean dose as a percentage of prescription dose
V20_lung = sum(d_lung >= px * 20 / 100) / numel(lung) * 100; % Volume receiving at least 20% of prescription dose

% Heart calculations
d_heart = d(heart);
Dmean_heart = mean(d_heart) / px * 100;
V30_heart = sum(d_heart >= px * 30 / 100) / numel(heart) * 100;

% Esophagus calculations
d_esoph = d(esoph);
Dmean_esoph = mean(d_esoph) / px * 100;

% % Spinal cord calculations
% d_cord = d(cord);
% Dmean_cord = mean(d_cord) / px * 100;

% Whole body dose calculation
d_body = d(body);
Dmean_body = mean(d_body) / px * 100;

end
