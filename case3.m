% This source code is written to implement the Christmas tree
% Author: wei-fan
% Email: zwf_jaccount@sjtu.edu.cn
% Open Source License: GPL

% shall we go!
clear all
clc
% simulation paraments set up
dt=0.01;
stime=50;
loop=stime/dt;

% flocking paraments set up
d=3;
n=2;

% init state
s=zeros(12,n);
rand = unifrnd(-pi,pi,[1,n]); % init random position on the circle
for i=1:n
   s(1:3,i) = [5*pi*cos(rand(i));5*pi*sin(rand(i));0]; 
end
s(4:6,:)=unifrnd(-0,0,[3,n]);
s(7,:)=unifrnd(-0.0*pi,0.0*pi,[1,n]);
s(8,:)=unifrnd(-0.0*pi,0.0*pi,[1,n]);
s(9,:)=unifrnd(-1*pi,1*pi,[1,n]);
x=s(1);y=s(2);z=s(3);
vx=s(4);vy=s(5);vz=s(6);
phi=s(7);theta=s(8);psi=s(9);
vphi=s(10);vtheta=s(11);vpsi=s(12);

% public virtual leadr init
xl=[5*pi;0;0];
vl=[0;0;0];

%parameters for quadrotor (this can be identified for crazyflie)
para.g=9.8;
para.m=1.2;
para.Iy=0.05;
para.Ix=0.05;
para.Iz=0.1;
para.b=10^-4;
para.l=0.5;
para.d=10^-6;
para.Jr=0.01;
para.k1=0.02;
para.k2=0.02;
para.k3=0.02;
para.k4=0.1;
para.k5=0.1;
para.k6=0.1;
para.omegaMax=330;

% history capture
xyHis=zeros(d,n+1,loop+1);
xyHis(:,:,1)=[xl s(1:3,:)];

%simulation start
hwait=waitbar(0,'Starting>>>>>>>>>>');

height = 1;
beta = pi/1000;
L = n*diag(ones(1,n)) - ones(n,n);
step = 0.2;
alpha = zeros(1,n);
valpha = zeros(1,n);
omegaHis=zeros(4,n,loop);
for t=1:loop
    
    %leader information generator
    if t/loop<0.05 %takeoff
        al=([0;0;height]-vl);
    elseif t/loop<0.99 %circle
        al=([-5*sin(beta*(t-0.05*loop));-5*cos(beta*(t-0.05*loop));0]-vl);
    else 
        al=([0;0;-height]-vl);
    end

    vl=vl+dt*al;
    xl=xl+dt*vl;
    % get x and v by using agreement protocol
    if t/loop<0.05 %takeoff
        x_1 = [s(1:2,1);xl(3)];
        x_2 = [s(1:2,2);xl(3)];
        v_1 = [s(4:5,1);vl(3)];
        v_2 = [s(4:5,2);vl(3)];
    elseif t/loop<0.99 %circle
        alpha(1) = (s(1:3,1)'*xl)/sqrt(sum(xl.^2))/sqrt(sum(s(1:3,1).^2));
        alpha(2) = (s(1:3,2)'*xl)/sqrt(sum(xl.^2))/sqrt(sum(s(1:3,2).^2));
        valpha = -alpha*L;
        alpha = valpha*step + alpha;
        x_1 = [5*pi*cos(alpha(1)-atan(xl(2)/xl(1)));5*pi*sin(alpha(1)-atan(xl(2)/xl(1)));xl(3)];
        x_2 = [5*pi*cos(alpha(2)-atan(xl(2)/xl(1)));5*pi*sin(alpha(2)-atan(xl(2)/xl(1)));xl(3)];
        v_1 = [valpha(1)*5*pi*sin(alpha(1)-atan(xl(2)/xl(1)));valpha(1)*5*pi*cos(alpha(1)-atan(xl(2)/xl(1)));vl(3)];
        v_2 = [valpha(2)*5*pi*sin(alpha(2)-atan(xl(2)/xl(1)));valpha(2)*5*pi*cos(alpha(2)-atan(xl(2)/xl(1)));vl(3)];
    else 
        x_1 = [s(1:2,1);xl(3)];
        x_2 = [s(1:2,2);xl(3)];
        v_1 = [s(4:5,1);vl(3)];
        v_2 = [s(4:5,2);vl(3)];
    end
    
    % get motor speeds form the controller 
    omega1=quadrotor_controller(s(:,1),x_1,v_1,0,para,1,10);
    omega2=quadrotor_controller(s(:,2),x_2,v_2,0,para,1,10);
    
    %record the speeds
    omegaHis(:,1,t)=omega1;
    omegaHis(:,2,t)=omega2;
    
    %send speeds of four motors to quadrotor and get its state
    s(:,1)=quadrotor_dynamic(s(:,1),omega1,para,dt);
    s(:,2)=quadrotor_dynamic(s(:,2),omega2,para,dt);
    
    %recodrd the position of quadrotor at time t/loop*stime
    xyHis(:,:,t+1)=[xl s(1:3,:)];
    
    waitbar(t/loop,hwait,'loading');
end

close(hwait);
%show the animation of the flight process
figure(1)
plotHis3(xyHis,dt,-1,200)
axis equal
grid on

%show changes in motor speeds during the flight
%figure(2)
%plot(omegaHis')%,'linesmooth','on')
%grid on
