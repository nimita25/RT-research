clc;clear;close all;
ptid='9306087';
cubeDim=[511,511,80];
id=[45 135 225 315];
weight=[1 0.2 0.2 1]; % [1,1] - no weight;
nY=cubeDim(1)*cubeDim(2)*cubeDim(3);
idx0=259;idy0=253; % tumor center
id0=30; % slice to center distance, (!) change with angle
id1=30;
dim1=2*id1+1;
slice_Dim=[dim1,cubeDim(3)];
nZ=slice_Dim(1)*slice_Dim(2);
for i=1:numel(id)
    if (id(i)>0 && id(i)<90)
        idx1=idx0-id0-id1;
        idy1=idy0+id0-id1;
        idx2=idx0-id0+id1;
        idy2=idy0+id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
    elseif (id(i)>90 && id(i)<180)
        idx2=idx0+id0-id1;
        idy2=idy0+id0+id1;
        idx1=idx0+id0+id1;
        idy1=idy0+id0-id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
    elseif (id(i)>180 && id(i)<270)
        idx1=idx0+id0-id1;
        idy1=idy0-id0-id1;
        idx2=idx0+id0+id1;
        idy2=idy0-id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
    elseif (id(i)>270 && id(i)<360)
        idx2=idx0-id0+id1;
        idy2=idy0-id0-id1;
        idx1=idx0-id0-id1;
        idy1=idy0-id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
    end
    inta=[];
    intb=[];
    intc=[];
    for k=1:nZ
        [id_x,id_z]=ind2sub(slice_Dim,k);
        j=sub2ind(cubeDim,ceil(iddx(id_x)),ceil(iddy(id_x)),id_z);
        inta=[inta k];
        intb=[intb j];
        intc=[intc (iddx(id_x)-ceil(iddx(id_x))+1)*(iddy(id_x)-ceil(iddy(id_x))+1)];
        if (ceil(iddy(id_x))>iddy(id_x))
            j=sub2ind(cubeDim,ceil(iddx(id_x)),floor(iddy(id_x)),id_z);
            inta=[inta k];
            intb=[intb j];
            intc=[intc (iddx(id_x)-ceil(iddx(id_x))+1)*(ceil(iddy(id_x))-iddy(id_x))];
        end
        if (ceil(iddx(id_x))>iddx(id_x))
            j=sub2ind(cubeDim,floor(iddx(id_x)),ceil(iddy(id_x)),id_z);
            inta=[inta k];
            intb=[intb j];
            intc=[intc (ceil(iddx(id_x))-iddx(id_x))*(iddy(id_x)-ceil(iddy(id_x))+1)];
        end
        if (ceil(iddy(id_x))>iddy(id_x) && ceil(iddx(id_x))>iddx(id_x))
            j=sub2ind(cubeDim,floor(iddx(id_x)),floor(iddy(id_x)),id_z);
            inta=[inta k];
            intb=[intb j];
            intc=[intc (ceil(iddx(id_x))-iddx(id_x))*(ceil(iddy(id_x))-iddy(id_x))];
        end
    end
    eval(['intp_',num2str(id(i)),'=sparse(inta,intb,intc,slice_Dim(1)*slice_Dim(2),nY);'])
    if (i==1)
        eval(['intp=weight(i)*intp_',num2str(id(i)),';'])
    else
        eval(['intp=[intp;weight(i)*intp_',num2str(id(i)),'];'])
    end
end
save([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'],['intp']);
for i=1:numel(id)
    save([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'],['intp_',num2str(id(i))],'-append');
end
save([ptid '\' ptid '_' num2str(numel(id)) '_intp.mat'],['slice_Dim'],'-append');