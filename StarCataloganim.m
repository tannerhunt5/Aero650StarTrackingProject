%Written by Neil McHenry
% April 18th 2018
%StarCataloganim.m

clc;
clear;
close all;

%Initialization 
    % Initial attitude: 
%         rho = [1.27, 0.28, -2.7]';
         rho = [0,0,0]';

        
    % Attitude dynamics: pure spin about z-axis with w = 0.01 rad/s;
        w = 0.05;
    % Set Max simulation Time
        maxTime = 10;
        f = 1/10;
    % Select the visible stars falling in the star tracker FOV.
    % Add multiplicative Gaussian noise to observed stars: 3sigma = 20";
        gridSize = 200;
        
        %Initialize Surfaces and gauss plot histories
        gauss_surf = zeros(gridSize);
        dot_gauss1 = zeros(gridSize);
        dot_gauss2 = zeros(gridSize);
        dot_gauss3 = zeros(gridSize);
        dot_gauss4 = zeros(gridSize);
        dot_gauss5 = zeros(gridSize);
        pause_interval = 0;
        %Sigma is the spread value of the gaussian function
        sigma = 3;
        q = [];
        
        %Integration time/Exposure Steps
        Ig = 5;
        
        %Loop Count
        i=0;
        %Loop Count past each division of 5 timesteps
        i5 = i-5;
        

        d=zeros(2,2,2);
        
        %Observation Matrices
            Xhist = [];
            Yhist = [];
        
        

%% Part A: Determination of Observed Stars
    %Load Stars with magnitude threshold 6
    [VS,R] = catalogsearch([-1.5,6]);
    
    %Determine optical axis pointing direction
        opt_i = [0,0,1]';
    %Rotate Z axis to the current direction with DCMs
    
    for time=0:f:maxTime %This is the time in tenths of seconds

        i = i + 1;
        % Convert to DCM
%         rho(3) = rho(3) + (0.0005*time);
%         rho(2) = rho(2) + (0.001*time);

        Att = rhotoC(rho);
        %Calculate final ptical axis direction
        opt_f = Att*opt_i;
        
        omega = w*4*time;
%         %Calculate DCM for rotation about y axis
%             C =  [cos(omega), 0, sin(omega)
%                    0    , 1,    0  
%                  -sin(omega), 0,cos(omega)];

        %Calculate DCM for rotation about z axis
            C =  [cos(omega), -sin(omega), 0
                  sin(omega), cos(omega), 0
                  0,    0, 1  
                 ];
             
        %Rotate Optical axis vector about y axis
        opt_f = C*opt_f;
        
        
    %Compute vectors orthogonal to Optical Axis
        V = null(opt_f(:).') ;
        V = 0.08715*V;
        o = opt_f;
        


    %Create matrix of star positions in 3D and cut it in half
        %Z Segmentation
        %Turn into zeros based on condition
        Z = R(4,:);
        Z(Z < 0) = 0;
        
        %Eliminate the same elements within other rows
        LOGICZ = Z~=0;
        Y = LOGICZ.*R(3,:);
        X = LOGICZ.*R(2,:);
        
        %Replace zero columns with empty columns
        Z(Z == 0) = [];
        Y(Y == 0) = [];
        X(X == 0) = [];  
        %Combine row vectors of coordinates into a matrix
        Stars3D = [X;Y;Z];
        
        %Rotate whole catalog about Z axis
        Stars3D = C*Stars3D;
        
   %% Part B: Plotting Starmap, Star Tracker View
       
        %Plot Optical Axis (Star Tracker direction)
        subplot(1,3,1);
        Loptical = plot3([0,opt_f(1)],[0,opt_f(2)],[0,opt_f(3)], 'b-');
        hold on
        
%         %Plot Orthogonal vectors
%         plot3([o(1),V(1,2)+o(1)],[o(2),V(2,2)+o(2)],[o(3),V(3,2)+o(3)])
%         plot3([o(1),V(1,1)+o(1)],[o(2),V(2,1)+o(2)],[o(3),V(3,1)+o(3)])
        
        %3D Scatterplot of stars
        cataplot = scatter3(Stars3D(1,:),Stars3D(2,:),Stars3D(3,:),'.k');
        xlabel('X Axis'); 
        ylabel('Y axis');
        zlabel('Z axis');    
    
    %Plot FOV circle in 3D
        center = opt_f';
        radius = 1;
        thetaa=0:0.01:2*pi;
        points=repmat(center',1,size(thetaa,2))+radius* ...
            (V(:,1)*cos(thetaa)+V(:,2)*sin(thetaa));
        FOV_Circle = plot3(points(1,:),points(2,:),points(3,:),'b-');
    
    %Plot of 2D view from star tracker
        subplot(1,3,2);
        %Project Points onto plane defined by otrhogonal basis V
        projection = 11.47446930579460699942627653471*(Stars3D')*V;
        
        %XY Segmentation
            %Turn elements that meet a condition into zeros 
            %Combine into (X^2+Y^2)^(1/2) vector
            Xlin = (projection(:,1))';
            Ylin = (projection(:,2))';
            
            %Save Observations
            Xhist(i) = Xlin(1);
            Yhist(i) = Ylin(1);
            
            XYsq = ((Xlin.*Xlin)+(Ylin.*Ylin)).^(1/2);
            XYsq(XYsq > (0.08715)) = 0;

            %Eliminate the same elements within other rows
            LOGICXY = XYsq~=0;
            Ylin = LOGICXY.*Ylin;
            Xlin = LOGICXY.*Xlin;

            %Replace zero columns with empty columns
            Ylin(Ylin == 0) = [];
            Xlin(Xlin == 0) = [];    
        
                hold on
                %Draw Circle with Radius 0.08715 at origin
                th = 0:pi/50:2*pi;
                xunit = 0.08715 * cos(th) + 0;
                yunit = 0.08715 * sin(th) + 0;
                FOV_Circle2 = plot(xunit, yunit, 'b-');    

            
        %Scatterplot of projection
        scatterstar = scatter(Xlin,Ylin,'r*');
    
        title('Star Tracker View (Projected)');
        xlabel('X-Prime Axis'); 
        ylabel('Y-Prime axis');

        %Gaussian Noise Definition
        subplot(1,3,3);
        gauss = 3*sigma;
        
        %Start Position to plot gaussian center
        for count = 1:length(Xlin)
            Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
            Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

            [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
            dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
            lgaussx = length(dot_gauss) +Sx - 1;
            lgaussy = length(dot_gauss) +Sy - 1;
            
            %Save history of last 5 dot_gauss plots
            if i == 1
                dot_gauss1(Sx:lgaussx,Sy:lgaussy) = dot_gauss1(Sx:lgaussx,Sy:lgaussy) + dot_gauss;
            end
            if i == 2
                dot_gauss2(Sx:lgaussx,Sy:lgaussy) = dot_gauss2(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;                
            end
            if i == 3
                dot_gauss3(Sx:lgaussx,Sy:lgaussy) = dot_gauss3(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;                
            end
            if i == 4
                dot_gauss4(Sx:lgaussx,Sy:lgaussy) = dot_gauss4(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;                
            end
            if i == 5
                dot_gauss5(Sx:lgaussx,Sy:lgaussy) = dot_gauss5(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;                
            end
           
            gauss_surf(Sx:lgaussx,Sy:lgaussy) = gauss_surf(Sx:lgaussx,Sy:lgaussy) + dot_gauss;

        end
        
      %If time is greater than 5, start removing old dot_gauss values
            i5 = i5+ 1;
      
            if i5 == 1
                %Subtract Old Plot
                gauss_surf = gauss_surf -dot_gauss1;
                %Reset dot_gauss with zeros
                dot_gauss1 = zeros(200);
                %Save New Plot
                for count = 1:length(Xlin)
                    Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
                    Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

                    [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
                    dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
                    lgaussx = length(dot_gauss) +Sx - 1;
                    lgaussy = length(dot_gauss) +Sy - 1;                
                    dot_gauss1(Sx:lgaussx,Sy:lgaussy) = dot_gauss1(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;
                end

            end

            if i5 == 2
                gauss_surf = gauss_surf - dot_gauss2;
                dot_gauss2 = zeros(200);
                for count = 1:length(Xlin)
                    Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
                    Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

                    [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
                    dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
                    lgaussx = length(dot_gauss) +Sx - 1;
                    lgaussy = length(dot_gauss) +Sy - 1;                
                    dot_gauss2(Sx:lgaussx,Sy:lgaussy) = dot_gauss2(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;
                end
            end

            if i5 == 3
                gauss_surf = gauss_surf - dot_gauss3;
                dot_gauss3 = zeros(200);
                for count = 1:length(Xlin)
                    Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
                    Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

                    [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
                    dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
                    lgaussx = length(dot_gauss) +Sx - 1;
                    lgaussy = length(dot_gauss) +Sy - 1;                
                    dot_gauss3(Sx:lgaussx,Sy:lgaussy) = dot_gauss3(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;
                end
            end

            if i5 == 4
                gauss_surf = gauss_surf - dot_gauss4;
                dot_gauss4 = zeros(200);
                for count = 1:length(Xlin)
                    Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
                    Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

                    [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
                    dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
                    lgaussx = length(dot_gauss) +Sx - 1;
                    lgaussy = length(dot_gauss) +Sy - 1;                
                    dot_gauss4(Sx:lgaussx,Sy:lgaussy) = dot_gauss4(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;
                end
            end

            if i5 == 5
                gauss_surf = gauss_surf - dot_gauss5;
                dot_gauss5 = zeros(200);
                for count = 1:length(Xlin)
                    Sx = (gridSize - round(1000*(Xlin(count)+0.1))) - gauss;
                    Sy = (round(1000*(Ylin(count)+0.1))) - gauss;

                    [xgrid,ygrid] = meshgrid (-gauss:gauss, -gauss:gauss);
                    dot_gauss = exp(-1/(2*sigma^2) * (xgrid.^2 + ygrid.^2));
                    lgaussx = length(dot_gauss) +Sx - 1;
                    lgaussy = length(dot_gauss) +Sy - 1;                
                    dot_gauss5(Sx:lgaussx,Sy:lgaussy) = dot_gauss5(Sx:lgaussx,Sy:lgaussy)+ dot_gauss;
                end
                %Reset Count
                i5 = 0;
            end                   
        
        
        surfplot = surf(gauss_surf);
        view([-90 90]);
        shading interp
        title('Star Tracker View (Gaussian Noise)');
        ylabel('X-Prime Axis'); 
        xlabel('Y-Prime axis');

        set(gcf, 'Position', [200, 200, 1400, 400])
        hold off
        
        %Interval to update the plot (pause in ms)
        pause(pause_interval);
%             xlim([-1 1])
%             ylim([-1 1])
%             zlim([-1 1])  

        delete(cataplot)
        delete(Loptical)
        delete(FOV_Circle)
        delete(FOV_Circle2)
        delete(scatterstar)
        delete(surfplot)
        
        numStars(i)= length(Xlin);

      hold off   
    end
   
    hold on
        subplot(1,3,1)
        cataplot = scatter3(Stars3D(1,:),Stars3D(2,:),Stars3D(3,:),'.k');
        FOV_Circle = plot3(points(1,:),points(2,:),points(3,:),'b-');
        Loptical = plot3([0,opt_f(1)],[0,opt_f(2)],[0,opt_f(3)], 'b-');
        subplot(1,3,2)
        scatterstar = scatter(Xlin,Ylin,'r*');
        FOV_Circle2 = plot(xunit, yunit, 'b-'); 
        subplot(1,3,3)
        surfplot = surf(gauss_surf);
        view([-90 90]);
        shading interp
    hold off

