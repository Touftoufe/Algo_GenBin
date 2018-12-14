clear all
close all
clc

load('matlab.mat');
plot(lambda,J,lambda,f(lambda,0.722296, 6562.565300, 10.774000, 0.218124),'r');
figure(2)
plot(lambda,JC,lambda,f(lambda,0.722296, 6562.565300, 10.774000, 0.218124),'r');