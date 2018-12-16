%clear all
close all
clc

load('matlab.mat');
plot(lambda,J,lambda,f(lambda,c(1), c(2), c(3), c(4)),'r');
figure(2)
plot(lambda,JC,lambda,f(lambda,c(1), c(2), c(3), c(4)),'r');