function ptv = ctv2ptv_080720(ctv, r, cubeDim, cubeRes)
% CTV2PTV_080720 Converts CTV to PTV using expansion.
%
% INPUT:
%   ctv       - Indices related to the CTV.
%   r         - Expansion radius in millimeters.
%   cubeDim   - Dimensions of the cube.
%   cubeRes   - Voxel resolution in millimeters.
%
% OUTPUT:
%   ptv       - Indices related to the PTV.

% Extract cube dimensions and voxel resolutions
nx = cubeDim(1);
ny = cubeDim(2);
nz = cubeDim(3);
dx = cubeRes.x;
dy = cubeRes.y;
dz = cubeRes.z;

% Convert CTV indices to coordinates
[ctvx, ctvy, ctvz] = ind2sub([nx, ny, nz], ctv);

N = numel(ctv);

% Initialize PTV mask
mask_ptv = zeros([nx, ny, nz]);

% Calculate PTV using expansion
for i = 1:N
    xc = max(1, ctvx(i) - ceil(r / dx)):min(ctvx(i) + ceil(r / dx), nx);
    yc = max(1, ctvy(i) - ceil(r / dy)):min(ctvy(i) + ceil(r / dy), ny);
    zc = max(1, ctvz(i) - ceil(r / dz)):min(ctvz(i) + ceil(r / dz), nz);
    [xc, yc, zc] = ndgrid(xc, yc, zc);
    id = find(dx^2 * (ctvx(i) - xc).^2 + dy^2 * (ctvy(i) - yc).^2 + dz^2 * (ctvz(i) - zc).^2 <= r^2);
    id2 = sub2ind([nx, ny, nz], xc(id), yc(id), zc(id));
    mask_ptv(id2) = 1;
end

% Find indices related to the PTV
ptv = find(mask_ptv == 1);
end
