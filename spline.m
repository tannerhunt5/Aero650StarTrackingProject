clc;
close all;

%Initialization
%fm input is a 3x2 matrix of the forces and moments e.g.:
    Forces = [.2; 0; .3];
    Moments = [20; .4; -2];
    fm = [Forces,Moments];

% dim and Inrt are inputs for Ihwd for a solid cuboid of height h, 
%width w, depth d, mass m, (1x4 and 1x3) e.g.:

    h = 0.8;  %0.5m
    w = 0.5;  %0.5m
    d = 0.5;  %0.5m
    m = 8.2;  %8.2kg
    dim = ([h,w,d,m]);

    Ih = 1/12*m*(w^2+d^2);
    Iw = 1/12*m*(d^2+h^2);
    Id = 1/12*m*(w^2+h^2);
    ihwd = ([Ih,Iw,Id]);
    
    %Number of Steps
    npts = 10;
[omega1,omega2,omega3,pos] = bodyDynamics(1:npts,fm,dim,ihwd);
pos = pos';
plot3(pos(1,:),pos(2,:),pos(3,:),'ro','LineWidth',2);

text(pos(1,:),pos(2,:),pos(3,:),[repmat('  ',npts,1), num2str((1:npts)')])

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
box on

hold on
% fnplt(cscvn(xyz(:,[1:end 1])),'r',2)
fnplt(cscvn(pos),'r',2)

hold off