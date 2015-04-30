% Load data and make bw
clear all;close all; clc; 
set(0,'Defaultfigurewindowstyle','docked')

I = imread('football.jpg');
Ibw = im2bw(I,0.95);
Ibw = not(Ibw);

figure(1);clf
imagesc(Ibw);colormap(gray)

%% Calculate axis and draw

[M N] = size(Ibw);
[X Y] = meshgrid(1:N,1:M);

%Mass and mass center
m = sum(sum(Ibw));
x0 = sum(sum(Ibw.*X))/m;
y0 = sum(sum(Ibw.*Y))/m;

%Covariance matrix elements
Mxx = sum(sum((X-x0).^2.*Ibw))/m;
Myy = sum(sum((Y-y0).^2.*Ibw))/m;
Mxy = sum(sum((Y-y0).*(X-x0).*Ibw))/m;

MM = [Mxx Mxy; Mxy Myy];

[U S V] = svd(MM);

W = V(:,1)/sign(V(1,1)); %Extremal directions (normalized to have first coordinate positive)
H = V(:,2);
W = 2*sqrt(S(1,1))*W; %Scaling of extremal directions to give ellipsis half axis
H = 2*sqrt(S(2,2))*H;

figure(1)
hold on
    plot(x0,y0,'r*');
    quiver(x0,y0,W(1),H(1),'r')
    quiver(x0,y0,W(2),H(2),'r')
hold off