ptid='9306087';
folder = 'C:\Users\nshinde\Desktop\pMBRT\';
load([folder ptid '\' ptid '.mat'],'cst','ct');
load([folder ptid '\' 'dij_' ptid '_doseGrid111.mat']);
body=cst{1,4}{1};
ctv2=cst{38,4}{1};
margin=5; % 5 mm inward margin
ctv1=ctv2ptv_080720(ctv2,margin,ct.cubeDim,ct.resolution);

% Generate id for 113 grid resolution
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv2) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv2 = VdoseGrid;
tmpCube=zeros(ct.cubeDim);
tmpCube(ctv1) = 1;
VdoseGrid=find(matRad_interp3(ct.x,ct.y,ct.z,tmpCube, ...
                            doseGrid.x,doseGrid.y',doseGrid.z,'nearest'));
ctv1 = VdoseGrid;

dr=10;ddr=3;
[Vpeak,Vvalley,Vlattice,Plattice]=generate_lattice(doseGrid,ctv1,ctv2,dr,ddr,ptid);