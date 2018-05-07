% Opening and Reading Hipparcos star catalog

fid = fopen('hip_main.dat','r');
T = readtable('hip_main.dat','Delimiter','|');
fclose(fid);

RA_deg = table2array(T(:,'Var9'));
Dec_deg = table2array(T(:,'Var10'));
C = table2array(T(:,'Var8'));
C = 1;
Mag = table2array(T(:,'Var6'));
CartCoords = zeros(length(RA_deg(:,1)), 3);
Store = 0;
for i = 1:length(RA_deg(:,1))
    if Mag(i) < 5
    
        X = (C * cos(Dec_deg(i))) * cos(RA_deg(i));
        Y = (C * cos(Dec_deg(i))) * sin(RA_deg(i));
        Z = C * sin(Dec_deg(i));
        CartCoords(i,1) = X;
        CartCoords(i,2) = Y;
        CartCoords(i,3) = Z;
        
        TwoDx(i,1) = Y/X;
        TwoDx(i,2) = Z/X;
    
        
    else
        CartCoords(i,1) = 0;
        CartCoords(i,2) = 0;
        CartCoords(i,3) = 0;
    end

    
end
CartCoords(CartCoords == 0) = [];
CartCoords = reshape(CartCoords, [], 3);

rad = .1;
[x,y,z] = sphere(1000);

figure(1)
set(gcf, 'Position', [300 300 1200 550])
mesh(x*rad,y*rad,z*rad); hold on;
plot3(CartCoords(:,1), CartCoords(:,2), CartCoords(:,3), '.k');
grid on;

% Describing the camera and finding stars within the Camera's FOV
FOVx = 10; %degs
FOVy = 10; %degs
foc = 600; %mm 
OptAxis = [0 0 1];
rho = [2.5, 1.5, -2];
FOVnew = [];
sigma = 20/3/3600*pi/180;

% Selecting stars in FOV
Ctrue = Gibbs2DCM(rho);

% Direction of Z in the inertial frame
zInertnew = Ctrue*[0 0 1]';
quiver3(0,0,0, zInertnew(1), zInertnew(2),zInertnew(3))

for i = 1:length(CartCoords(:,1))

    rstar = CartCoords(i,:);
    ctheta = dot(rstar,zInertnew)/(norm(rstar)*norm(zInertnew));
    Thetanew = acos(ctheta)*180/pi;

    if Thetanew <= FOVx
        FOVnew = [FOVnew; i]; 
    end
 
end
NumStars = length(FOVnew);

% Vectors towards captured stars in FOV
for i = 1:NumStars
    quiver3(0,0,0, CartCoords(FOVnew(i),1),CartCoords(FOVnew(i),2),CartCoords(FOVnew(i),3),'r');
end
hold off;

% Converting Points from world space to camera space
CamCoords = [];
for i = 1:NumStars
    CamCoords(i,:) = CartCoords(FOVnew(i),:)*Ctrue;
end

%%
% Each row represents a star, and the position of that star over a number of time steps is equal to the number of
% columns
xp = [];
yp = [];
for i = 1:NumStars
     
    xp(i) = (-foc/(CamCoords(i,3)))*(CamCoords(i,1)); % Sample star streaks[-10 -5 0 5; 10 15 20 25; -60 -55 -50 -45; 70 75 80 85; 40 45 50 55];
    yp(i) = (-foc/(CamCoords(i,3)))*(CamCoords(i,2));% Sample star streaks[50 55 60 65; 15 20 25 35; -75 -70 -65 -60; -65 -60 -55 -50; 60 65 70 75];
    
end
steps = 3*length(xp);
pdfplotter(xp, yp, steps)

function [] = pdfplotter(xp, yp, steps)
   
    Sigma = [steps 0; 0 steps];
    x1 = -300:2:300; 
    x2 = -300:2:300;
    A = zeros(length(x1),length(x2));

    [X1,X2] = meshgrid(x1,x2);
    NumStars = length(xp(:,1));
    
    for j = 1:NumStars
        
        for i = 1:length(xp(1,:))

            F{i} = mvnpdf([X1(:)+xp(j,i) X2(:)+yp(j,i)], [], Sigma);
            F{i} = reshape(F{i}, length(x2), length(x1));
            A = A + F{i};
            
        end
        
        figure(3)
        surf(x1, x2, A); hold on;
        view(2); title('Simulated CCD');
        set(gcf, 'Position', [0 0 1338 1057]);
        xlabel('X, Pixels'); ylabel('Y, Pixels'); %zlabel('z');
        
        
    end
end





