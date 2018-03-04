% This source code is written to implement agreement protocol
% Author: wei-fan
% Email: zwf_jaccount@sjtu.edu.cn
% Open Source License: GPL

function [x,y] = case1( t,step )
%note: step should be around 0.1~0.3 in order to let the curves converged
%and smooth.

%init position
x = zeros(t,4);
y = zeros(t,4);
x(1,:) = normrnd(10,2,1,4);
y(1,:) = normrnd(10,2,1,4);

%init weight graphs for both x,y direction
weight_x = zeros(4,4);
weight_x(1,2) = rand(1);
%weight_x(2,1) = weight_x(1,2);
%weight_x(1,3) = rand(1);
%weight_x(3,1) = weight_x(1,3);
weight_x(1,4) = rand(1);
%weight_x(4,1) = weight_x(1,4);
weight_x(3,2) = rand(1);
%weight_x(2,3) = weight_x(3,2);
weight_x(2,4) = rand(1);
weight_x(4,2) = weight_x(2,4);
weight_x(3,4) = rand(1);
%weight_x(4,3) = weight_x(3,4);
weight_x
weight_y = zeros(4,4);
weight_y(1,2) = rand(1);
weight_y(2,1) = weight_y(1,2);
weight_y(1,3) = rand(1);
weight_y(3,1) = weight_y(1,3);
weight_y(1,4) = rand(1);
weight_y(4,1) = weight_y(1,4);
weight_y(2,3) = rand(1);
weight_y(3,2) = weight_y(2,3);
weight_y(2,4) = rand(1);
weight_y(4,2) = weight_y(2,4);
weight_y(3,4) = rand(1);
weight_y(4,3) = weight_y(3,4);

%calculate delta matrix and Laplacian matrix
delta_x = zeros(1,4);
delta_y = zeros(1,4);
for i=1:4
    for j=1:4
        delta_x(i) = delta_x(i) + weight_x(i,j);
        delta_y(i) = delta_y(i) + weight_y(i,j);
    end
end

delta_x = diag(delta_x);
delta_y = diag(delta_y);
L_x = delta_x - weight_x;
L_y = delta_y - weight_y;

eig(L_x)
%eig(L_y)

%iteration
for i=2:t
    vx = -L_x*x(i-1,:)';
    x(i,:) = vx'*step+ x(i-1,:);
    vy = -L_y*y(i-1,:)';
    y(i,:) = vy'*step+ y(i-1,:);
end

%show results
t_p = 1:t;

figure('Name','x')
hold on;
grid on;
plot(t_p,x(:,1));
plot(t_p,x(:,2));
plot(t_p,x(:,3));
plot(t_p,x(:,4));

figure('Name','y')
hold on;
grid on;
plot(t_p,y(:,1));
plot(t_p,y(:,2));
plot(t_p,y(:,3));
plot(t_p,y(:,4));

figure('Name','trajectory')
hold on;
grid on;
plot(x(:,1),y(:,1));
plot(x(:,2),y(:,2));
plot(x(:,3),y(:,3));
plot(x(:,4),y(:,4));
end
