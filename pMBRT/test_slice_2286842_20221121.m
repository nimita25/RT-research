id=[45 135 225 315];
px=1;
cubeDim=[622,622,154];
slice_Dim=[622,154];
iddz=1:slice_Dim(2);
load("2286842_c.mat");
if 1
    d=zeros(59580136,1);
    d(c{2})=0.5;
    d(c{1})=1;
%     d(c{4})=0.2;
%     d(c{5})=0.4;
%     d(c{7})=0.6;
%     d(c{8})=0.8;
else
    d=d(:);
end
idx0=265;idy0=316; % tumor center
id1=30;
dim1=2*id1+1;
for i=1:numel(id)
    if (id(i)<90)
        id0=25; % slice to center distance, (!) change with angle
        idx1=idx0-id0-id1;
        idy1=idy0+id0-id1;
        idx2=idx0-id0+id1;
        idy2=idy0+id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
        d(sub2ind(cubeDim,repmat(ceil(iddx),1,cubeDim(3)),repmat(ceil(iddy),1,cubeDim(3)),repelem(iddz,dim1)))=px;
    elseif (id(i)<180)
        id0=40; % slice to center distance, (!) change with angle
        idx2=idx0+id0-id1;
        idy2=idy0+id0+id1;
        idx1=idx0+id0+id1;
        idy1=idy0+id0-id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
        d(sub2ind(cubeDim,repmat(ceil(iddx),1,cubeDim(3)),repmat(ceil(iddy),1,cubeDim(3)),repelem(iddz,dim1)))=px;
    elseif (id(i)<270)
        id0=40; % slice to center distance, (!) change with angle
        idx1=idx0+id0-id1;
        idy1=idy0-id0-id1;
        idx2=idx0+id0+id1;
        idy2=idy0-id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
        d(sub2ind(cubeDim,repmat(ceil(iddx),1,cubeDim(3)),repmat(ceil(iddy),1,cubeDim(3)),repelem(iddz,dim1)))=px;
    else
        id0=25; % slice to center distance, (!) change with angle
        idx2=idx0-id0+id1;
        idy2=idy0-id0-id1;
        idx1=idx0-id0-id1;
        idy1=idy0-id0+id1;
        iddx=linspace(idx1,idx2,dim1);
        iddy=linspace(idy1,idy2,dim1);
        d(sub2ind(cubeDim,repmat(ceil(iddx),1,cubeDim(3)),repmat(ceil(iddy),1,cubeDim(3)),repelem(iddz,dim1)))=px;
    end
end
d=reshape(d,cubeDim);
figure;imshow3D(d,[0 px]);