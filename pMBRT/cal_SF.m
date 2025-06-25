function sf = cal_SF(d,type)
alpha=-0.211;
beta=-0.068;
if (type==1)
    sf=mean(exp(alpha*d+beta*mlq(0.15*d).*d.^2));
else
    sf=mean(exp(alpha*d+beta*d.^2));
end
end