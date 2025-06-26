clc;clear;close all;
load protons_Generic_DoNotChange.mat machine
for i=1:114
machine.data(i).Z=machine.data(i).Z/25;
end
figure;hold on;plot(machine.data(14).Z);plot(machine.data(48).Z);plot(machine.data(81).Z);hold off
save protons_Generic.mat machine

