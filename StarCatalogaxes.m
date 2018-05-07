%Written by Neil McHenry
% April 18th 2018
%StarCataloganim.m

clc;
clear;
close all;

%Initialization 
    % Initial attitude: 
%         rho = [1.27, 0.28, -2.7]';
         rho = [0.3,0.1,0]';

        
    % Attitude dynamics: pure spin about y-axis with w = 0.01 rad/s;
        w = 0.01;
    % Set Max simulation Time
        maxTime = 10;
        f = 1/10;
    % Select the visible stars falling in the star tracker FOV.
    % Add multiplicative Gaussian noise to observed stars: 3sigma = 20";
        gridSize = 200;
        gauss_surf = zeros(gridSize);
%         gauss_hist = [];
        pause_interval = 0;
        sigma = 2;
        q = [];
        i=0;
        
        %Observation Matrices
            Xhist = [];
            Yhist = [];
        
        

%% Part A: Determination of Observed Stars
    %Load Stars with magnitude threshold 6
    [VS,R] = catalogsearch([-1.5,6]);
    

    %Rotate Z axis to the current direction with DCMs
    
    for time=0:f:maxTime %This is the time in tenths of seconds

        i = i + 1;
        % Convert to DCM
%         rho(2) = rho(2) + (0.00008*time);
%         rho(3) = rho(3) + (0.004*time);
%         rho(1) = rho(1) + (0.001*time);
        %Calculate final ptical axis direction
    %Determine optical axis pointing direction
        opt_i = [0,0,1]';
    %Rotate Optical Axis to Initial Orientation     
         Att = rhotoC(rho);
    %Calculate New Optical Axis
         opt_f = Att*opt_i;
        
        
        
%         if time < 5
            omega = w*10*time;
%         else
%             omega = w*time;
%         end

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
        
   %% Part B: Plotting Starmap, Star Tracker View
       
        %Plot Optical Axis (Star Tracker direction)
        subplot(1,3,1);
        
        Loptical = plot3([0,opt_f(1)],[0,opt_f(2)],[0,opt_f(3)], 'b-');
        hold on
        title('3D Hemispherical Plot of Catalog');

%         %Plot Orthogonal vectors
        line1 = plot3([o(1),V(1,2)+o(1)],[o(2),V(2,2)+o(2)],[o(3),V(3,2)+o(3)], 'r');
        line2 = plot3([o(1),V(1,1)+o(1)],[o(2),V(2,1)+o(2)],[o(3),V(3,1)+o(3)], 'g');
        
        set(line1,'LineWidth',3);
        set(line2,'LineWidth',3);
        
        %3D Scatterplot of stars
        cataplot = scatter3(Stars3D(1,:),Stars3D(2,:),Stars3D(3,:),'.k');
        view([-180 70]);

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
        
        line1p = plot([0,0],[0,.08175], 'r');
                set(line1p,'LineWidth',3);
        line2p = plot([-.08175,0],[0,0], 'g');
                set(line2p,'LineWidth',3);
        
                

       
        
        %XY Segmentation
            %Turn elements that meet a condition into zeros 
            %Combine into (X^2+Y^2)^(1/2) vector
            Xlin = (projection(:,1))';
            Xproj = Xlin;

            Ylin = (projection(:,2))';
            Yproj = Ylin;

            
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
%         scatterstar = scatter(Xlin,Ylin,'r*');
        scatterstar = scatter(-Xproj,Yproj,'.k');

    
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


            gauss_surf(Sx:lgaussx,Sy:lgaussy) = gauss_surf(Sx:lgaussx,Sy:lgaussy) + dot_gauss;
%             gauss_hist(i) = gauss_surf;
        end
        surfplot = surf(gauss_surf);
        view([-90 90]);
        shading interp
        title('Star Tracker View (Gaussian Noise)');
        ylabel('X-Prime Axis'); 
        xlabel('Y-Prime axis');

        set(gcf, 'Position', [200, 200, 1600, 400])
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
        delete(line1)
        delete(line2)
        
        numStars(i)= length(Xlin);

      hold off   
    end
   
    hold on
        subplot(1,3,1)
        cataplot = scatter3(Stars3D(1,:),Stars3D(2,:),Stars3D(3,:),'.k');
        view([-180 70]);
        FOV_Circle = plot3(points(1,:),points(2,:),points(3,:),'b-');
        Loptical = plot3([0,opt_f(1)],[0,opt_f(2)],[0,opt_f(3)], 'b-');
        line1 = plot3([o(1),V(1,2)+o(1)],[o(2),V(2,2)+o(2)],[o(3),V(3,2)+o(3)],'r');
        line2 = plot3([o(1),V(1,1)+o(1)],[o(2),V(2,1)+o(2)],[o(3),V(3,1)+o(3)],'g');
        subplot(1,3,2)
%         scatterstar = scatter(Xlin,Ylin,'r*');
        scatterstar = scatter(-Xproj,Yproj,'.k');
        FOV_Circle2 = plot(xunit, yunit, 'b-'); 
        subplot(1,3,3)
        surfplot = surf(gauss_surf);
        view([-90 90]);
        shading interp
    hold off

