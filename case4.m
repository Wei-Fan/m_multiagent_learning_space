% This source code is written to implement the Christmas tree
% Author: wei-fan
% Email: zwf_jaccount@sjtu.edu.cn
% Open Source License: GPL

function [] = case4( n,stime )
clc
%init
dt = 0.01;
loop = stime/dt;
beta = pi/20;
L = n*diag(ones(1,n)) - ones(n,n);
x = zeros(loop+1,n,3);
v = zeros(loop+1,n,3);
alpha = zeros(loop+1,n,1);
dalpha = zeros(loop+1,n,1);
valpha = zeros(loop+1,n,1);
pp = -n/2:n/2;
for i=1:n
   alpha(1,i,1) = unifrnd(0,pi);
   %x(1,i,:) = [5*cos(alpha(1,i,1));5*sin(alpha(1,i,1));0];
   x(1,i,:) = [(10+10/n*pp(i))*cos(alpha(1,i,1));(10+10/n*pp(i))*sin(alpha(1,i,1));0];
   %v(1,i,:) = [-5*beta.*sin(alpha(1,i,1));5*beta.*cos(alpha(1,i,1));0];
   v(1,i,:) = [0,0,0];
end

% init state
s=zeros(12,n);
for i=1:n
   s(1:3,i) = [x(1,i,1);x(1,i,2);x(1,i,3)]; 
   s(4:6,i) = [v(1,i,1);v(1,i,2);v(1,i,3)];
end
s(7,:)=unifrnd(-0.0*pi,0.0*pi,[1,n]);
s(8,:)=unifrnd(-0.0*pi,0.0*pi,[1,n]);
s(9,:)=unifrnd(-1*pi,1*pi,[1,n]);
%x=s(1);y=s(2);z=s(3);
%vx=s(4);vy=s(5);vz=s(6);
%phi=s(7);theta=s(8);psi=s(9);
%vphi=s(10);vtheta=s(11);vpsi=s(12);

%init virtue leader 
alpha_ref = zeros(loop+1,1);
alpha_ref(1,1) = mean(alpha(1,:,1));
dalpha(1,:,1) = alpha(1,:,1) - alpha_ref(1,1);
x_ref = [10*cos(alpha_ref(1,1));10*sin(alpha_ref(1,1));0];
v_ref = [0;0;0];
%a_ref = [0;0;0];
var(alpha(1,:,1))
step = 0.01/var(alpha(1,:,1))/n^2/beta %important!!!converge speed parameter

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
xyHis=zeros(3,n+1,loop+1);
xyHis(:,:,1)=[x_ref s(1:3,:)];

%simulation start
hwait=waitbar(0,'Starting>>>>>>>>>>');

omegaHis=zeros(4,n,loop);
omega = zeros(4,n,loop);
for t=1:loop
    if t/loop < 0.1 %take off
        a_ref=([0;0;2]-v_ref);
        v_ref=v_ref+dt*a_ref;
        x_ref=x_ref+dt*v_ref;
        %pause agreement protocl
        alpha_ref(t+1,1) = alpha_ref(t,1);
        dalpha(t+1,:,1) = dalpha(t,:,1);
        alpha(t+1,:,1) = alpha(t,:,1);
        for i = 1:n
            x(t+1,i,:) = [x(t,i,1);x(t,i,2);x_ref(3)*(1-1/n*pp(i))];
            v(t+1,i,:) = [v(t,i,1);v(t,i,2);v_ref(3)*(1-1/n*pp(i))];
        end
    else %circle
        %generate leader trajectory
        alpha_ref(t+1,1) = alpha_ref(t,1) + beta*dt;
        %a_ref=([-5*beta*sin(alpha_ref(t+1,1));5*beta*cos(alpha_ref(t+1,1));0]-v_ref);
        %v_ref=v_ref+dt*a_ref;
        %x_ref=x_ref+dt*v_ref;
        x_ref(1) = 10*cos(alpha_ref(t+1,1));
        x_ref(2) = 10*sin(alpha_ref(t+1,1));
        %v_ref = [-5*beta*sin(alpha_ref(t+1,1));5*beta*cos(alpha_ref(t+1,1));0];

        % get x and v by using agreement protocol
        dalpha(t,:,1) = alpha(t,:,1) - alpha_ref(t,1);
        valpha(t,:,1) = -dalpha(t,:,1)*L;
        dalpha(t+1,:,1) = valpha(t,:,1)*step + dalpha(t,:,1);

        alpha(t+1,:,1) = dalpha(t+1,:,1) + alpha_ref(t+1,1);
        %x(t+1,:,:) = [5*cos(alpha(t+1,:,1));5*sin(alpha(t+1,:,1));zeros(1,n)];
        for i = 1:n
            x(t+1,i,:) = [(10+10/n*pp(i))*cos(alpha(t+1,i,1));(10+10/n*pp(i))*sin(alpha(t+1,i,1));x_ref(3)*(1-1/n*pp(i))];
            v(t+1,i,:) = [-(10+10/n*pp(i))*(valpha(t+1,i,1)*step/dt+beta).*sin(alpha(t+1,i,1));(10+10/n*pp(i))*(valpha(t+1,i,1)*step/dt+beta).*cos(alpha(t+1,i,1));0];
            %v(t+1,i,:) = [-5*beta*sin(alpha(t+1,i,1));5*beta*cos(alpha(t+1,i,1));0];
       end
        end
    
    % get motor speeds form the controller and record the speeds
    for i = 1:n
        x_tmp = [x(t+1,i,1);x(t+1,i,2);x(t+1,i,3)];
        v_tmp = [v(t+1,i,1);v(t+1,i,2);v(t+1,i,3)];
        omega(:,i,t) = quadrotor_controller(s(:,i),x_tmp,v_tmp,0,para,1,10);
        omegaHis(:,i,t) = omega(:,i,t);
    end
    
    %send speeds of four motors to quadrotor and get its state
    for i = 1:n 
        s(:,i)=quadrotor_dynamic(s(:,i),omega(:,i,t),para,dt);
    end
    
    %recodrd the position of quadrotor at time t/loop*stime
    xyHis(:,:,t+1)=[x_ref s(1:3,:)];
    
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

%figure('Name','dalpha')
%hold on;
%grid on;
%for i=1:n
%    plot(1:loop+1,dalpha(:,i,1));
%end

figure('Name','x')
hold on;
grid on;
for i=1:n
    %plot(1:loop+1,x(:,i,1));
    plot(x(:,i,1),x(:,i,2))
end

%figure('Name','v')
%hold on;
%grid on;
%for i=1:n
%    plot(1:loop+1,v(:,i,1));
%    plot(1:loop+1,v(:,i,2));
%    plot(1:loop+1,v(:,i,3));
    %plot(v(:,i,1),v(:,i,2))
%end

end

