clc;clear;close all;
ptid='2286842';
id=[45 135 225 315]; % in this script, only to put tv matrices of all angles together
load([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'],'slice_Dim');
nZ=slice_Dim(1)*slice_Dim(2);
inta=[];
intb=[];
intc=[];
intd=[];
inte=[];
intf=[];
for i=1:nZ
    [id_x,id_z]=ind2sub(slice_Dim,i);
    inta=[inta i];
    intb=[intb i];
    intc=[intc 1];
    if (id_x>1)
        inta=[inta i];
        intb=[intb sub2ind(slice_Dim,id_x-1,id_z)];
        intc=[intc -1];
    end
    intd=[intd i];
    inte=[inte i];
    intf=[intf 1];
    if (id_z>1)
        intd=[intd i];
        inte=[inte sub2ind(slice_Dim,id_x,id_z-1)];
        intf=[intf -1];
    end
end
tvx=sparse(inta,intb,intc,nZ,nZ);
tvz=sparse(intd,inte,intf,nZ,nZ);
% v=ones(numel(id),1);
v=[1;2;2;1];
m=diag(v);
tvx=kron(m,tvx);
tvz=kron(m,tvz);
save([ptid '\' ptid '_' num2str(numel(id)) '_tvm.mat'],'tvx','tvz');