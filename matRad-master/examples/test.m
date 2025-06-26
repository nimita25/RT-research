D = dij.physicalDose{1};
[nX,nY] = size(D);
d = D*ones(nY,1);
figure()
A = squeeze(max(d, [], 3));
imagesc(A)
d1 = reshape(d,ct.cubeDim);
line1 = d1(:,:,40);
figure()
plot(line1,'r','linewidth',2)