function plotdvh(d0,var)

n_c=var.n_c;
dmax=var.dmax;

n=100;
t=(0:1/(n-1):1)*dmax;
ptv1=zeros(n,1);
body=zeros(n,1);
% bladder=zeros(n,1);
% rectum=zeros(n,1);
% femhead_lt=zeros(n,1);
% femhead_rt=zeros(n,1);
% penilebulb=zeros(n,1);
for i=1:n
    ptv1(i)=numel(find(d0{1}>=t(i)))/n_c(1);
    body(i)=numel(find(d0{2}>=t(i)))/n_c(2);
%     bladder(i)=numel(find(d0{3}>=t(i)))/n_c(3);w
%     rectum(i)=numel(find(d0{4}>=t(i)))/n_c(4);
%     femhead_lt(i)=numel(find(d0{5}>=t(i)))/n_c(5);
%     femhead_rt(i)=numel(find(d0{6}>=t(i)))/n_c(6);
%     penilebulb(i)=numel(find(d0{7}>=t(i)))/n_c(7);
end

figure;hold on;
plot(t,ptv1,'r--');
plot(t,body,'k--');
% plot(t,bladder,'b');
% plot(t,rectum,'b');
% plot(t,femhead_lt,'g');
% plot(t,femhead_rt,'g');
% plot(t,penilebulb,'y');
% hold off
% drawnow
