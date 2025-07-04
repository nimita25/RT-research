ptid='9306087'; 
gridDim = '111';
addpath('../pMBRT');
method=1; %1: pMBRTL-2, 2: pMBRTL-1, 3: pMBRTL-1e, 4: conv
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);



ctv1 = cst{38,4}{1};   % PTV_59.5 Hao, 45Gy in 3 fractions
body = cst{1,4}{1};
N_oar = 5;
oar = cell(N_oar,1);
oar(1) = cst{37,4}(1); % Brainstem Dmax<15Gy; Dmax<54Gy
oar(2) = cst{11,4}(1); % OpticChiasm Dmax<10Gy; Dmax<36Gy
oar(3) = cst{12,4}(1); % OpticNrv_R Dmax<10Gy; Dmax<54Gy
oar(4) = cst{13,4}(1); % OpticNrv_L Dmax<10Gy; Dmax<54Gy
oar(5) = cst{44,4}(1); % Brain V12<5cc; Dmax<60Gy


N_iter=20;

%% Generate peak and valley
if method <= 3

lattice_x = [249;249;249;259;259;259;259;269;269;269;269;259;259;]; %13 peaks with best result for 8 angles
lattice_y = [238;248;268;238;248;258;268;238;248;268;278;228;278];
lattice_z = 79*ones(numel(lattice_x),1);%[79;79;79;79;79;79;79;79;79;79];
dr = 10;
ddr = 1.5;
margin = 5;
else
lattice_x = [259;259];
lattice_y = [238;268];
lattice_z = 79*ones(numel(lattice_x),1);
dr = 30;
ddr = 10;
margin = 5;
end

%peakid = setdiff(peakid,[20482154]);

%% Set indices of CTV, peak, valley according to gridDim. Note default gridDim = 3x3x3
ctv1=ctv2ptv_080720(ctv1,margin,ct.cubeDim,ct.resolution);
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv1) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv1 = VdoseGrid;
if 0%exist('9306087_lattice.mat','file')
    load('9306087_lattice.mat');
else
    %[lattice_x,lattice_y,lattice_z]  = [peakidx,peakidy,peakidz];% ind2sub(doseGrid.cubeDim,peakid);
    N_obj = numel(lattice_x);
    Vlattice=cell(N_obj,1);
    %disp(N_obj)
    for i=1:N_obj
        Vlattice(i)={{generate_id(doseGrid.x(lattice_x(i)),doseGrid.y(lattice_y(i)),doseGrid.z(lattice_z(i)),ddr,doseGrid)}};
    end
    Vpeak=[Vlattice{:,1}];
    Vpeak=unique(vertcat(Vpeak{:}));
    Vvalley=setdiff(ctv1,Vpeak);
    Plattice=[lattice_x,lattice_y,lattice_z];
    save('9306087_lattice.mat', 'Vpeak','Vvalley','Vlattice','Plattice','dr','ddr');
end



%% Set objective weights
ctv=cell(1,1);
ctv{1}=ctv1;
n_oar=zeros(N_oar,1);
for i=1:N_oar
    n_oar(i)=numel(oar{i});
end
c=[{Vpeak};{Vvalley};{body};oar;];
%c=[{Vpeak};{mod_Vvalley};{setdiff(Vvalley,mod_Vvalley)};{body};oar;];

N_c=numel(c);
n_c=zeros([N_c 1]);
for i=1:N_c
    n_c(i)=numel(c{i});
end

px=2;nfrac=10;px0=px*nfrac;pvdr = 5;

N_obj = 13;
w_obj = [1; 100; 1; 1; 1; 100; 0.1; 1; 1; 1; 1; 1; 1;];
s_obj = [pvdr*px; pvdr*px; pvdr*px * 1.1; px; px; px * 1.1; 0; px; [15; 10; 10; 10; 12]/nfrac];
n_obj = round([nan; n_c(1) * 0.95; nan; nan; n_c(2) * 0.95; nan; nan; nan; nan; nan; nan; nan; 5 / 0.3^3;]);
c_obj = [1; 1; 1; 2; 2; 2; 3; 3; 4; 5; 6; 7; 8;];
id_obj = cell(N_obj, 1);
type_obj = [0; 2; 3; 0; 2; 3; 0; 3; 3; 3; 3; 3; 1;];



%% Load dij
%id=[90 180]; % (depend on setup)
id=[0 45 90 135 180 225 270 315];
nN=numel(id);
ctc = [7]; %[3 5 7];


if method == 1

    coll_angle = [0 90];
    Dij = {};
    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            if isempty(Dij)
                Dij = {[dij.physicalDose{1}]};
            else
                Dij = {[Dij{1},dij.physicalDose{1}]};
            end
        end
    end

elseif method == 2
    
    coll_angle = [0];
    Dij = {};
    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
           
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
            if isempty(Dij)
                Dij = {[dij.physicalDose{1}]};
            else
                Dij = {[Dij{1},dij.physicalDose{1}]};
            end
        end
    end
elseif method == 3
    
    Dij = {};
    for idc = 1:numel(id)
        if mod(idc,2) == 0
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(0) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
        else
            load([folder  'dij_' ptid '_' gridDim '_collimator' num2str(90) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat']);
        end
        if isempty(Dij)
            Dij = {[dij.physicalDose{1}]};
        else
            Dij = {[Dij{1},dij.physicalDose{1}]};
        end
    end
else
    Dij = {};
    for idc = 1:numel(id)
        load([folder  'dij_' ptid '_' gridDim '_'  num2str(id(idc)) '.mat']);
        if isempty(Dij)
            Dij = {[dij.physicalDose{1}]};
        else
            Dij = {[Dij{1},dij.physicalDose{1}]};
        end
    end
    
end

[nY,nX]=size(Dij{1});

if method == 1

    coll_angle = [0 90];

    nf = numel(coll_angle)*numel(id);

    pe = cell(nf,1);
    eg = cell(nf,1);
    pe_unique = cell(nf,1);
    ii = 0;

    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
            
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat'],'dij' , 'stf');
            i = (ca-1)*numel(id)+idc;
            disp(i)
                pe{i} = zeros([stf.totalNumOfBixels 1]);
                nr = numel(stf.numOfBixelsPerRay);
                ie = 0;
                for j = 1:nr
                    ne = stf.numOfBixelsPerRay(j);
                    pe{i}(ie+(1:ne)) = stf.ray(j).energy;
                    ie = ie + ne;
                end
                pe_unique{i} = unique(pe{i});
                eg(i) = {cell([numel(pe_unique{i}),1])};
                for j = 1:numel(pe_unique{i})
                    eg{i}{j} = find(pe{i} == pe_unique{i}(j)) + ii;
                end
                ii = ii + stf.totalNumOfBixels;
            
        end
    end

    ne_unique = 0;
    pe2 = [];
    
    for i = 1:nf
        ne_unique = ne_unique+numel(pe_unique{i});
        pe2 = [pe2;pe{i}];
    end
    id_gs = zeros([nX 1]);
    n_gs = zeros([ne_unique 1]);
    ii = 0; jj = 0;
    for i = 1:nf
        for j = 1:numel(pe_unique{i})
            ii = ii + 1;
            n_gs(ii) = numel(eg{i}{j});
            id_gs(jj + (1:n_gs(ii))) = eg{i}{j};
            jj = jj + n_gs(ii);
        end
    end

elseif method == 2
    
    coll_angle = [0];

    nf = numel(coll_angle)*numel(id);

    pe = cell(nf,1);
    eg = cell(nf,1);
    pe_unique = cell(nf,1);
    ii = 0;

    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)
            
            load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(coll_angle(ca)) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat'],'dij' , 'stf');
            i = (ca-1)*numel(id)+idc;
            disp(i)
                pe{i} = zeros([stf.totalNumOfBixels 1]);
                nr = numel(stf.numOfBixelsPerRay);
                ie = 0;
                for j = 1:nr
                    ne = stf.numOfBixelsPerRay(j);
                    pe{i}(ie+(1:ne)) = stf.ray(j).energy;
                    ie = ie + ne;
                end
                pe_unique{i} = unique(pe{i});
                eg(i) = {cell([numel(pe_unique{i}),1])};
                for j = 1:numel(pe_unique{i})
                    eg{i}{j} = find(pe{i} == pe_unique{i}(j)) + ii;
                end
                ii = ii + stf.totalNumOfBixels;
            
        end
    end

    ne_unique = 0;
    pe2 = [];
    
    for i = 1:nf
        ne_unique = ne_unique+numel(pe_unique{i});
        pe2 = [pe2;pe{i}];
    end
    id_gs = zeros([nX 1]);
    n_gs = zeros([ne_unique 1]);
    ii = 0; jj = 0;
    for i = 1:nf
        for j = 1:numel(pe_unique{i})
            ii = ii + 1;
            n_gs(ii) = numel(eg{i}{j});
            id_gs(jj + (1:n_gs(ii))) = eg{i}{j};
            jj = jj + n_gs(ii);
        end
    end

elseif method == 3

    coll_angle = [0];
    
    nf = numel(coll_angle)*numel(id);

    pe = cell(nf,1);
    eg = cell(nf,1);
    pe_unique = cell(nf,1);
    ii = 0;

    for ca = 1:numel(coll_angle)
        for idc = 1:numel(id)

            if mod(idc,2) == 0
                load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(0) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat'],'dij' , 'stf');
            else
                load([folder 'dij_' ptid '_' gridDim '_collimator' num2str(90) '_ctc' num2str(ctc) 'mm_' num2str(id(idc)) '.mat'],'dij' , 'stf');
            end
            
            i = (ca-1)*numel(id)+idc;
            disp(i)
                pe{i} = zeros([stf.totalNumOfBixels 1]);
                nr = numel(stf.numOfBixelsPerRay);
                ie = 0;
                for j = 1:nr
                    ne = stf.numOfBixelsPerRay(j);
                    pe{i}(ie+(1:ne)) = stf.ray(j).energy;
                    ie = ie + ne;
                end
                pe_unique{i} = unique(pe{i});
                eg(i) = {cell([numel(pe_unique{i}),1])};
                for j = 1:numel(pe_unique{i})
                    eg{i}{j} = find(pe{i} == pe_unique{i}(j)) + ii;
                end
                ii = ii + stf.totalNumOfBixels;
            
        end
    end

    ne_unique = 0;
    pe2 = [];
    
    for i = 1:nf
        ne_unique = ne_unique+numel(pe_unique{i});
        pe2 = [pe2;pe{i}];
    end
    id_gs = zeros([nX 1]);
    n_gs = zeros([ne_unique 1]);
    ii = 0; jj = 0;
    for i = 1:nf
        for j = 1:numel(pe_unique{i})
            ii = ii + 1;
            n_gs(ii) = numel(eg{i}{j});
            id_gs(jj + (1:n_gs(ii))) = eg{i}{j};
            jj = jj + n_gs(ii);
        end
    end

    
else
    disp('no energy layer info for conventional')
    % Dij = {};
    % for idc = 1:numel(id)
    %     load([folder ptid '/dij_' ptid '_' gridDim '_'  num2str(id(idc)) '.mat']);
    %     if isempty(Dij)
    %         Dij = {[dij.physicalDose{1}]};
    %     else
    %         Dij = {[Dij{1},dij.physicalDose{1}]};
    %     end
    % end
    
end


if method == 1
    load('output_9306087/res_9306087_8_1_2024-12-03-21-18.mat')
elseif method == 2
    load('output_9306087/res_9306087_8_2_2024-12-03-23-55.mat')
elseif method == 3
    load('output_9306087/res_9306087_8_3_2024-12-04-00-58.mat')
elseif method == 4
    % load('output_9306087/res_9306087_8_4_2024-12-04-04-09.mat')
    load('output_9306087/res_9306087_8_4_2024-12-04-04-25.mat')
end


if method <= 3
    MU_val = sum(x)/5;
    ii=0;
    active_EL=0;
    for i_gs = 1:numel(n_gs)
        if nnz(x(id_gs(ii+(1:n_gs(i_gs))))) > 0
            active_EL=active_EL+1;
        end
        ii= ii+n_gs(i_gs);
    end
else
    MU_val = sum(x);
    active_EL=numel(id)*20;
end
%disp(active_EL)

disp('Total MU value (in Gp) =')
disp(MU_val)
disp('Using dose rate of 400 Giga-proton per second: https://www.sciencedirect.com/science/article/pii/S0360301615001297')
disp('Beam on time (in secs) =')
disp(MU_val/400)
disp('Total gantry setup time = number of angles*30secs')
disp('Total collimator rotation time = 30secs')
disp('Energy layer switching time (1 secs per layer)')
EL_switch =1;
disp(active_EL*EL_switch)
disp('Total time (in secs)=')
if method == 1 || method == 3
    Tot_time =(MU_val/400)+(16*30)+(30)+(active_EL*EL_switch);
else
    Tot_time =(MU_val/400)+(8*30)+(0)+(active_EL*EL_switch);
end
disp(Tot_time)
disp('Total time (in mins) =')
disp(Tot_time/60)

for i=3:N_c
tmpCube=zeros(ct.cubeDim);
tmpCube(c{i}) = 1;
% interpolate cube
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
c{i} = VdoseGrid;
end

for i=1:N_c
    n_c(i)=numel(c{i});
end
body=setdiff(c{3},c{1});
body=setdiff(body,c{2});
for i = [1 2 3 4 8]
disp([max(d(c{i})) mean(d(c{i}))])
end
disp(max(d(body)))