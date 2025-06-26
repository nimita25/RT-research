
close all

% figure,imshowpair(ct.cubeHU{1}(:,:,120),resultGUI.physicalDose(:,:,120),'blend');
% return

ix=130+(-50:50);
iy=130+(-50:50);

x1=ct.cubeHU{1}(ix,iy,120);
x1=x1/max(x1(:));
x2=resultGUI.physicalDose(ix,iy,120);
x2=x2/max(x2(:));

E=x1;
figure, imshow(E), hold on
red = cat(3, ones(size(E)), zeros(size(E)), zeros(size(E)));
h = imshow(red,[-0.1,0.1]);
set(h, 'AlphaData', 0.75*x2);

return
figure;
% imagesc(ct.cubeHU{1}(:,:,120));
imshow(x1,[]);
hold on

h = imagesc(x2);colormap jet
set( h, 'AlphaData', .5 ); % .5 transparency
return

close all
figure;
imshow(ct.cubeHU{1}(:,:,120),[]);
hold on
% imshow(resultGUI.physicalDose(:,:,120),[]);

h = imagesc(resultGUI.physicalDose(:,:,120) ); % show the edge image
set( h, 'AlphaData', .5 ); % .5 transparency
colormap gray

return
hold on;
h = imagesc( edgeImg ); % show the edge image
set( h, 'AlphaData', .5 ); % .5 transparency
colormap gray
