%% Rotating Body
% Neil McHenry
% AERO 650, March 2018

clc 
close all
clear
%% Problem setup

m=1.2; %kg
R=0.024; %m
g=9.81;
Ia=m*R^2/2;
It=m*R^2/4;
P=[0;R; 0;];  %initial center from o on B frame

%% Case 1
W=[0; 1 ; 0]; % Initial angular rates, principal body axes, rad/s

Eulers=[0 ;deg2rad(89.7); 0];

Ainitial=[1 0 0;
      0 sin(Eulers(2)) 0;
      0 cos(Eulers(2)) 1];
  
eulerDot=Ainitial\W;

C1=[1 0 0;
       0 cos(Eulers(2)) sin(Eulers(2));
       0 sin(Eulers(2)) cos(Eulers(2))];

Pinertial=C1'*P;

ic=[W;Eulers;Pinertial];

tspan=[0 2];

%% Solver
options=odeset('RelTol', 1e-5);
%Solver
[t,y]=ode45(@(t,y) rotatingBodyEqs(t,y,m,R,g,Ia,It) ,tspan,ic,options);

%% Plots

figure(1)
subplot(2,2,1)
plot3(y(:,7),y(:,8),y(:,9))
hold on
plot3(y(2,7),y(2,8),y(2,9),'o')
grid on
axis equal
xlabel('X, m')
ylabel('Y, m')
zlabel('Z, m')
title('Inertial Position of Center of Mass')

subplot(2,2,2)
plot(y(:,7),y(:,8))
grid on
hold on
axis equal
xlabel('X, m')
ylabel('Y, m')
title('X-Y view, Position Center of Mass')

subplot(2,2,3)
hold on
plot(t, y(:,4))
plot(t, y(:,5))
xlabel('Time, s')
ylabel('Angle, rad')
title('Euler Angles of StarTracker')
legend('\phi', '\theta')

subplot(2,2,4)
hold on
plot(t, y(:,1))
plot(t, y(:,2))
plot(t, y(:,3))
xlabel('Time, s')
ylabel('Angle, rad/s')
title('Angular Velocities')
legend('\omega 1', '\omega 2','\omega 3')

