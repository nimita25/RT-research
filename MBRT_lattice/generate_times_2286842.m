ptid='2286842'; 
gridDim = '111'; % set the dimension of the grid. Default = 3x3x3 mm
addpath('../pMBRT');
method=1; %1: M2, 2: M0, 3: M1, 4: conv
% 1. load ct & define oars & optimization parameter (depend on ptid)
folder = ['../' ptid '/'];
load([folder ptid '.mat'], 'ct', 'cst');
%load([folder ptid '/' ptid '.mat'],'cst','ct');
load([folder 'dij_' ptid '_doseGrid' gridDim '.mat']);


ctv1=cst{22,4}{1};   % ptv69.96, 2.12*33
body=cst{21,4}{1};
N_oar=7; % (!)
oar=cell(N_oar,1);
oar(1)=cst{4,4}(1); % right parotid D50<30
oar(2)=cst{5,4}(1); % left parotid Dmean<26
oar(3)=cst{24,4}(1); % larynx Dmean<32 ??
oar(4)=cst{8,4}(1); % esophagus Dmean<22
oar(5)=cst{3,4}(1); % mandible Dmax<73 Dmean<50
oar(6)=cst{1,4}(1); % oral Dmean<35
oar(7)=cst{2,4}(1); % lip Dmean<15

N_iter=20;


%% Generate peak and valley
if method <= 3
% Set peak center acc to info from Dij --> max dose in 250 plane at 252,325
lattice_x = [252;262;272;252;252;262;262;272;272;272];
lattice_y = [325;315;305;305;315;305;325;315;325;335];
lattice_z = 250*ones(numel(lattice_x),1);%[79;79;79;79;79;79;79;79;79;79];
dr = 10; %min distance between vertex center
ddr = 1.5; %diameter of the vertex
margin = 5;
else %this defines peak vertices for conv. To do: set values acc to voxels in 228... case
lattice_x = [257;272];
lattice_y = [310;325];
lattice_z = 250*ones(numel(lattice_x),1);
dr = 30;
ddr = 8;
margin = 5;
end

%% Set indices of CTV, peak, valley according to gridDim. Note default gridDim = 3x3x3
ctv1=ctv2ptv_080720(ctv1,margin,ct.cubeDim,ct.resolution);
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv1) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv1 = VdoseGrid;
if 0%exist('2286842_lattice.mat','file')
    load('2286842_lattice.mat');
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
    Plattice=[lattice_x ,lattice_y,lattice_z];
    save('2286842_lattice.mat', 'Vpeak','Vvalley','Vlattice','Plattice','dr','ddr');
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


px=2.12;nfrac=33;px0=px*nfrac;pvdr = 5;

N_obj=17;
type_obj=[0;2;3;0;2;3;0;3;1;1;1;1;3;1;1;1;0];
w_obj=[1;10;1;1;10;1;0.1;1;1;1;1;1;1;1;1;1;0.1];
s_obj=[px0*pvdr;px0*pvdr;px0*1.1*pvdr;px0;px0;px0*1.1;0;px0;[30;26;32;22;73;50;35;15;0]]*px/px0;
n_obj=round([nan;n_c(1)*0.95;nan;nan;n_c(2)*0.95;nan;nan;nan;n_c(4)*0.5;n_c(5)*0.5;n_c(6)*0.5;n_c(7)*0.5;nan;n_c(8)*0.5;n_c(9)*0.5;n_c(10)*0.5;nan;]);
c_obj=[1;1;1;2;2;2;3;3;4;5;6;7;8;8;9;10;6];
id_obj=cell(N_obj,1);



%% Load dij
%id=[90 180]; % (depend on setup
% id=[0 45 90 135 180 225 270 315];
id = [45 135 225 315];
%id = [45 135];
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
    load('output_2286842/res_2286842_4_1_2024-12-07-11-18.mat')
elseif method == 2
    load('output_2286842/res_2286842_4_2_2024-12-07-07-38.mat')
elseif method == 3
    load('output_2286842/res_2286842_4_3_2024-12-07-09-43.mat')
elseif method == 4
    load('output_2286842/res_2286842_4_4_2024-12-07-23-52.mat')
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
    active_EL=numel(id)*33;
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
if method == 1 
    Tot_time =(MU_val/400)+(8*30)+(30)+(active_EL*EL_switch);
elseif method == 3
    Tot_time =(MU_val/400)+(4*30)+(30)+(active_EL*EL_switch);
else
    Tot_time =(MU_val/400)+(4*30)+(0)+(active_EL*EL_switch);
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
for i = [1 2 3 6 9]
disp([max(d(c{i})) mean(d(c{i}))])
end
disp(max(d(body)))