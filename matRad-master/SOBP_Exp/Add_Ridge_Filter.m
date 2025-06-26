clc,clear,close all;
load('PROSTATE.mat');
A = ct.cube{1};
ct_RF = ct;
[nX,nY,nZ] = size(A);
ct_RF.cube{1} = zeros(nX,nY,nZ);
Ty = 10;
Tx = 80;
l = 5;
ntri = 5;
nsize = 9;
cor = zeros((nsize + 1)^2/4,3);
choice = 0;
for i = 1:nZ
    B = zeros(nX,nY);
    id1 = nX - Tx;
    id2 = nY - Ty;
    temp1  = [];
    temp2  = [];
    temp3  = [];
    temp4  = [];
    for k = 1:(nsize + 1)/2
        temp1 = [temp1; ((id1+k):1:(id1+nsize-(k-1)))'];
        temp2 = temp1 + (nsize + 1)/2;
        temp3 = [temp3; repmat(id2 - (k - 1),nsize - 2*(k - 1),1)];
        temp4 = [temp4; repmat(id2 + k,nsize - 2*(k - 1),1)];
    end
    cor(:,1) = temp1;
    if choice == 1
        cor(:,2) = temp2;
    else
        cor(:,2) = temp1;
    end
    cor(:,3) = temp3;
    cor(:,4) = temp4;
    for j = 1:5
        for k = 1:(nsize + 1)^2/4
            B(cor(k,1) - (j - 1)*nsize,cor(k,3)) = 0;
            B(cor(k,2) - (j - 1)*nsize,cor(k,4)) = 0;
        end
    end
    ct_RF.cube{1}(:,:,i) = B;
end
save('CT_SOBP.mat','ct_RF')