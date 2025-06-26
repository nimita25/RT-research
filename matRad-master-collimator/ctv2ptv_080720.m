% close all;
function ptv=ctv2ptv_080720(ctv,r,cubeDim,cubeRes)
% nx=ct.cubeDim(1);
% ny=ct.cubeDim(2);
% nz=ct.cubeDim(3);
nx=cubeDim(1);
ny=cubeDim(2);
nz=cubeDim(3);

dx=cubeRes.x;
dy=cubeRes.y;
dz=cubeRes.z;

% dx=3;
% ctv=cst{10,4}{1};
% ptv0=cst{14,4}{1};
% mask=zeros([nx ny nz]);
% mask(ctv)=mask(ctv)+1;
% mask(ptv0)=mask(ptv0)+1;
% figure;imshow3D(mask,[]);

% shift=5;
% r=ceil(shift/dx);

[ctvx,ctvy,ctvz]=ind2sub([nx ny nz],ctv);

N=numel(ctv);
% tic
mask_ptv=zeros([nx ny nz]);
for i=1:N
    xc=max(1,ctvx(i)-ceil(r/dx)):min(ctvx(i)+ceil(r/dx),nx);
    yc=max(1,ctvy(i)-ceil(r/dy)):min(ctvy(i)+ceil(r/dy),ny);
    zc=max(1,ctvz(i)-ceil(r/dz)):min(ctvz(i)+ceil(r/dz),nz);
    [xc,yc,zc]=ndgrid(xc,yc,zc);
    id=find(dx^2*(ctvx(i)-xc).^2+dy^2*(ctvy(i)-yc).^2+dz^2*(ctvz(i)-zc).^2<=r^2);
    id2=sub2ind([nx ny nz],xc(id),yc(id),zc(id));
    mask_ptv(id2)=1;
end
% toc
ptv=find(mask_ptv==1);

% mask_ptv(ctv)=mask_ptv(ctv)+1;
% figure;imshow3D(mask_ptv,[]);
