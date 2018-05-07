
function [omegax,omegay,omegaz] = bodyDynamics(tarray)
% Written by Neil McHenry
% AERO 650

%Tarray input is a matrix of time intervals, e.g. [1:10]

Forces = [0; 0; 0];

Moments = [.2; .2; -4];

% Inertia for a solid cuboid of height h, width w, depth d, mass m

h = 0.5;  %0.5m
w = 0.5;  %0.5m
d = 0.5;  %0.5m
m = 1.2;  %1.2kg

Mass = m;

Ih = 1/12*m*(w^2+d^2);
Iw = 1/12*m*(d^2+h^2);
Id = 1/12*m*(w^2+h^2);


Inertia = [Ih,0,0;
          0, Iw ,0;
          0,0, Id];
%tarray = [1:10];

[t,y] = sixdof(Forces, Moments, Mass, Inertia, tarray)

%% Calculating Star Plot Projection due to Dynamics

%Part A: Convert euler rates & body position from Body frame to Inertial frame

%Part B: Calculate Projections based on Rotations and Displacement
    % Component 1: Dislocation in X due to X rotation

    % Component 2: Dislocation in Y due to Y rotation
    
    % Component 3: Dislocation in X & Y due to Z rotation
    
    % Component 4: Sum dislocations
    
 %Part C: Output Xp and Yp over time (projection based on orientation/position)

 omegax = y(:,1);
 omegay = y(:,2);
 omegaz = y(:,3);
 
end