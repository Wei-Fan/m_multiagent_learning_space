% This source code is the core of Christmas tree
% Author: wei-fan
% Email: zwf_jaccount@sjtu.edu.cn
% Open Source License: GPL

function [] = case2(n,stime)
%init
dt = 0.01;
loop = stime/dt;
beta = pi/5;
L = n*diag(ones(1,n)) - ones(n,n);
step = 0.05/n;
x = zeros(loop+1,n,3);
v = zeros(loop+1,n,3);
alpha = zeros(loop+1,n,1);
dalpha = zeros(loop+1,n,1);
valpha = zeros(loop+1,n,1);
for i=1:n
   alpha(1,i,1) = unifrnd(-pi,pi);
   x(1,i,:) = [5*cos(alpha(1,i,1));5*sin(alpha(1,i,1));0];
   v(1,i,:) = [-5*beta.*sin(alpha(1,i,1));5*beta.*cos(alpha(1,i,1));0];
end

%init virtue leader
alpha_ref = zeros(loop+1,1);

alpha_ref(1,1) = mean(alpha(1,:,1));
%dalpha(1,:,1) = alpha(1,:,1) - alpha_ref(1,1);

%start iteration
for t = 1:loop
    %generate leader trajectory
    alpha_ref(t+1,1) = alpha_ref(t,1) + beta*dt;
    
    %agreement protocol
    dalpha(t,:,1) = alpha(t,:,1) - alpha_ref(t,1);
    valpha(t,:,1) = -dalpha(t,:,1)*L;
    dalpha(t+1,:,1) = valpha(t,:,1)*step + dalpha(t,:,1);
    
    alpha(t+1,:,1) = dalpha(t+1,:,1) + alpha_ref(t+1,1);
    %x(t+1,:,:) = [5*cos(alpha(t+1,:,1));5*sin(alpha(t+1,:,1));zeros(1,n)];
    for i = 1:n
        x(t+1,i,:) = [5*cos(alpha(t+1,i,1));5*sin(alpha(t+1,i,1));0];
        v(t+1,i,:) = [-5*(valpha(t+1,i,1)+beta).*sin(alpha(t+1,i,1));5*(valpha(t+1,i,1)+beta).*cos(alpha(t+1,i,1));0];
    end
end

%show results

figure('Name','dalpha')
hold on;
grid on;
for i=1:n
    plot(1:loop+1,dalpha(:,i,1));
end

figure('Name','x')
hold on;
grid on;
for i=1:n
    %plot(1:loop+1,x(:,i,1));
    plot(x(:,i,1),x(:,i,2))
end

figure('Name','v')
hold on;
grid on;
for i=1:n
    plot(1:loop+1,v(:,i,1));
    %plot(v(:,i,1),v(:,i,2))
end

end
