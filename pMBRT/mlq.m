function r = mlq(d)
d=max(d,10^(-8));
r=2*(d+exp(-d)-1)./d.^2;
end