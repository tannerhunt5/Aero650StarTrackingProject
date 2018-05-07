function dy=rotatingBodyEqs(~,y,m,R,g,Ia,It)

dy=zeros(9,1);


%Euler angles
dy(5)=y(1);
dy(4)=y(2)/sin(y(5));
dy(6)=y(3)-y(2)*cot(y(5));

%Omegas
dy(1)=-((g*m*R*cos(y(5)) + y(2)*((Ia - It + m*R^2)*y(3) + It*dy(6)))/(It + m*R^2));
dy(2)=(y(1)*((Ia - It)*y(3) + It*dy(6)))/It;
dy(3)=(m*R^2*y(1)*y(2))/(Ia + m*R^2);

%Velocity
v1=-R*(y(3));
v2=0;
v3=R*y(1);

v=[v1;v2;v3];

%Position vector in inertial coordinates

C31=[cos(y(4)) sin(y(4)) 0;
       -cos(y(5))*sin(y(4)) cos(y(4))*cos(y(5)) sin(y(5));
       sin(y(4))*sin(y(5)) -cos(y(4))*sin(y(5)) cos(y(5))];

dy(7:9)=C31'*v;



