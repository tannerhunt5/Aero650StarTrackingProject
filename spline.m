clc;
close all;

%Initialization
%fm input is a 3x2 matrix of the forces and moments e.g.:
    Forces = [1; 10; 0];
    Moments = [0; .9; 0];
    fm = [Forces,Moments];
    Forces2 = [.6; 10; 0];
    Moments2 = [0; .5; 0];
    fm2 = [Forces2,Moments2];
    %Number of Steps
    npts = 50;
    
    %Set Optical Axis Pointing Direction relative to body
    opt_i = (.05*npts*norm(Forces))*[1,1,0]'; 
    
    %Hold each animation step for this many seconds
    pause_interval = .1;
% dim and Inrt are inputs for Ihwd for a solid cuboid of height h, 
%width w, depth d, mass m, (1x4 and 1x3) e.g.:

    h = 5;  %0.2m
    w = 5;  %0.1m
    d = 5;  %0.1m
    m = 100.2;
    m2 = 1000.2;

    dim = ([h,w,d,m]);
    dim2 = ([h,w,d,m2]);

    Ih = 1/12*m*(w^2+d^2);
    Iw = 1/12*m*(d^2+h^2);
    Id = 1/12*m*(w^2+h^2);
    ihwd = ([Ih,Iw,Id]);
    
    %Vector to store body rotations
    rot = [];
    rot2 = [];
    
%Calculate angular velocites and body position over time
    [omegaxyz,pos] = bodyDynamics(1:npts,fm,dim,ihwd);
    pos = pos';
    [omegaxyz2,pos2] = bodyDynamics(1:npts,fm2,dim2,ihwd);
    pos2 = pos2';   
    
    anglexyz = cumtrapz(1:npts,omegaxyz);
    anglexyz2 = cumtrapz(1:npts,omegaxyz2);    
       
%Calculate rotation (rot) at each time interval (first group)
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
        %Pause for a certain interval to play an animation
    end
    
       
%Calculate rotation (rot) at each time interval (second group)
    for i = 1:npts
        phi = anglexyz2(i,1);   %Roll
        theta = anglexyz2(i,2); %Pitch
        psi = anglexyz2(i,3);   %Yaw
        
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
        rot2 = [rot2,opt_f];
        %Pause for a certain interval to play an animation
    end    
    
%Plot Position Points  (red)  
    plot3(pos(1,:),pos(2,:),pos(3,:),'ro','LineWidth',2);
%Plot Point labels
    text(pos(1,:),pos(2,:),pos(3,:),[repmat('  ',npts,1), num2str((1:npts)')])

box on

hold on

%Plot Position Points2    (green)
    plot3(pos2(1,:),pos2(2,:),pos2(3,:),'go','LineWidth',2);
%Plot Point labels2
    text(pos2(1,:),pos2(2,:),pos2(3,:),[repmat('  ',npts,1), num2str((1:npts)')])


% fnplt(cscvn(xyz(:,[1:end 1])),'r',2)
%Plot Spline 1
fnplt(cscvn(pos),'r',2)
        xlabel('X Axis'); 
        ylabel('Y axis');
        zlabel('Z axis');  

%Plot Spline 2        
fnplt(cscvn(pos2),'g',2)
        xlabel('X Axis'); 
        ylabel('Y axis');
        zlabel('Z axis');          
%Plot optical axes for each point (group 1)
 for j= 1:npts
    plot3([pos(1,j),rot(1,j)+pos(1,j)], ...
          [pos(2,j),rot(2,j)+pos(2,j)], ...
          [pos(3,j),rot(3,j)+pos(3,j)], 'bp-');
    pause(pause_interval);
 end
 
%Plot optical axes for each point (group 2) 
 for j= 1:npts
    plot3([pos2(1,j),rot2(1,j)+pos2(1,j)], ...
          [pos2(2,j),rot2(2,j)+pos2(2,j)], ...
          [pos2(3,j),rot2(3,j)+pos2(3,j)], 'gp-');
    pause(pause_interval);
 end
hold off