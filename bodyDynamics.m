
function [omegaxyz,pos] = bodyDynamics(tarray,fm,dim,ihwd)
% Written by Neil McHenry
% AERO 650

%Tarray input is a matrix of time intervals, e.g. [1:10]
%fm input is a 3x2 matrix of the forces and moments e.g.:
    % Forces = [.2; 0; .3];
    % Moments = [20; .4; -2];
    % fm = [Forces,Moments];

% dim and Inrt are inputs for Ihwd for a solid cuboid of height h, 
%width w, depth d, mass m, (1x4 and 1x3) e.g.:

    % h = 0.8;  %0.5m
    % w = 0.5;  %0.5m
    % d = 0.5;  %0.5m
    % m = 8.2;  %8.2kg
    % dim = ([h,w,d,m]);

    % Ih = 1/12*m*(w^2+d^2);
    % Iw = 1/12*m*(d^2+h^2);
    % Id = 1/12*m*(w^2+h^2);
    % ihwd = ([Ih,Iw,Id]);

    %Assign Values based on Inputs
    Forces = fm(:,1);
    Moments = fm(:,2);
    Mass = dim(4);
    
    Ih = ihwd(1);
    Iw = ihwd(2);
    Id = ihwd(3);
    Inertia = [Ih,0,0;
              0, Iw ,0;
              0,0, Id];
%tarray = [1:10];

[t,y] = sixdof(Forces, Moments, Mass, Inertia, tarray);

%% Calculating Star Plot Projection due to Dynamics

%Part A: Convert euler rates & body position from Body frame to Inertial frame

%Part B: Calculate Projections based on Rotations and Displacement
    % Component 1: Dislocation in X due to X rotation

    % Component 2: Dislocation in Y due to Y rotation
    
    % Component 3: Dislocation in X & Y due to Z rotation
    
    % Component 4: Sum dislocations
    
 %Part C: Output Xp and Yp over time (projection based on orientation/position)

 omegaxyz = y(:,1:3);
 pos = y(:,7:9);

  
 
end