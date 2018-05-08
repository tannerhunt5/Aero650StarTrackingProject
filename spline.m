clc;
close all;

%Initialization
%fm input is a 3x2 matrix of the forces and moments e.g.:
    Forces = [1; 10; 0];
    Moments = [0; .9; 0];
    fm = [Forces,Moments];
    
    %Number of Steps
    npts = 70;
    
    %Set Optical Axis Pointing Direction relative to body
    opt_i = (.001*npts)*[1,1,0]';    
    
% dim and Inrt are inputs for Ihwd for a solid cuboid of height h, 
%width w, depth d, mass m, (1x4 and 1x3) e.g.:

    h = 0.5;  %0.5m
    w = 0.5;  %0.5m
    d = 0.5;  %0.5m
    m = 10000.2;  %100.2kg
    dim = ([h,w,d,m]);

    Ih = 1/12*m*(w^2+d^2);
    Iw = 1/12*m*(d^2+h^2);
    Id = 1/12*m*(w^2+h^2);
    ihwd = ([Ih,Iw,Id]);
    
    %Vector to store body rotations
    rot = [];
    
%Calculate angular velocites and body position over time
    [omegaxyz,pos] = bodyDynamics(1:npts,fm,dim,ihwd);
    pos = pos';
    
    anglexyz = cumtrapz(1:npts,omegaxyz);
       
%Calculate rotation (rot) at each time interval
    for i = 1:npts
        phi = anglexyz(i,1);   %Roll
        theta = anglexyz(i,2); %Pitch
        psi = anglexyz(i,3);   %Yaw
        
        %Calculate DCM for 321 euler rotation
      C =  [cos(theta)*cos(psi),cos(theta)*sin(psi),-sin(theta)
          
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), ...
            sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), ...
            sin(phi)*cos(theta)
            
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), ...
            cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), ...
            cos(phi)*cos(theta)];

        %Rotate Optical axis by euler angles
        opt_f = C*opt_i;
        rot = [rot,opt_f];
    end
    
%Plot Position Points    
    plot3(pos(1,:),pos(2,:),pos(3,:),'ro','LineWidth',2);
%Plot Point labels
    text(pos(1,:),pos(2,:),pos(3,:),[repmat('  ',npts,1), num2str((1:npts)')])

box on

hold on
% fnplt(cscvn(xyz(:,[1:end 1])),'r',2)
fnplt(cscvn(pos),'r',2)
        xlabel('X Axis'); 
        ylabel('Y axis');
        zlabel('Z axis');  
%Plot optical axes for each point
 for j= 1:npts
    plot3([pos(1,j),rot(1,j)+pos(1,j)], ...
          [pos(2,j),rot(2,j)+pos(2,j)], ...
          [pos(3,j),rot(3,j)+pos(3,j)], 'b-');
 end
hold off