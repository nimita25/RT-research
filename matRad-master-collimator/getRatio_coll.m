function stf = getRatio_coll(stf, pln)
% calculate the ratio for each ray

% prepare structures necessary for particles
fileName = [pln.radiationMode '_' pln.machine];
try
   load([fileparts(mfilename('fullpath')) filesep fileName]);
   SAD = machine.meta.SAD;
end
availableEnergies = [machine.data.energy];

for i = 1 : 1:length(pln.propStf.gantryAngles)
    for j = 1 : 1:stf(i).numOfRays
        stf(i).ray(j).ratio = zeros(size(stf(i).ray(j).energy));
        for k = 1 : stf(i).numOfBixelsPerRay(j)
            id = find(stf(i).ray(j).energy(k) == availableEnergies);
            sigma = machine.data(id).initFocus.sigma(1);
            theta = atan(SAD/sqrt(stf(i).ray(j).rayPos_bev(1)^2 + stf(i).ray(j).rayPos_bev(3)^2));
            center_x = (SAD - pln.coll_iso) * stf(i).ray(j).rayPos_bev(1) / SAD;
            center_y = (SAD - pln.coll_iso) * stf(i).ray(j).rayPos_bev(3) / SAD;
            area = 0;
            for n = 1 : size(pln.collimator, 1)
                ax = pln.collimator(n, 1);
                bx = pln.collimator(n, 2);
                ay = pln.collimator(n, 3);
                by = pln.collimator(n, 4);
                area = area + sqrt(2*pi)*sigma/(2 * sin(theta)) * (erf(sin(theta) * (bx - center_x) / sqrt(2) / sigma) - erf(sin(theta) * (ax - center_x) / sqrt(2) / sigma))...
                    * sqrt(2*pi)*sigma/(2 * sin(theta)) * (erf(sin(theta) * (by - center_y) / sqrt(2) / sigma) - erf(sin(theta) * (ay - center_y) / sqrt(2) / sigma))...
                    * 1 / (2 * pi * sigma^2);
            end
            stf(i).ray(j).ratio(k) = area;
        end
    end
end
end