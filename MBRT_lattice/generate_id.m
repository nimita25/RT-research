function id = generate_id(x,y,z,ddr,ct)
dim=ct.cubeDim;
id=[];
for i=1:dim(1)
    for j=1:dim(2)
        for k=1:dim(3)
            distance=sqrt((ct.x(i)-x)^2+(ct.y(j)-y)^2+(ct.z(k)-z)^2);
            if (distance<ddr)
                id=[id;sub2ind(dim,i,j,k)];
            end
        end
    end
end
end