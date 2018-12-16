%clear all
close all
clc

load('matlab.mat');

subplot(1,2,1)
plot(lambda,J,lambda,f(lambda,c(1), c(2), c(3), c(4)),'r');
xlabel("Lambda")
ylabel("J")
title("Approximation of the experimental data")

subplot(1,2,2)
plot(lambda,JC,lambda,f(lambda,c(1), c(2), c(3), c(4)),'r');
xlabel("Lambda")
ylabel("J")
title("Comparition with the approximation data provided in the 'profile.txt' file")
