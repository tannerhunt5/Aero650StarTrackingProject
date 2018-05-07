function [C] = Gibbs2DCM(rho)
% Converts the Gibbs/Rodrigues vector to the Direction Cosine Matrix
%   Detailed explanation goes here

Const = (1/(1+ rho(1)^2 + rho(2)^2 + rho(3)^2));

C(1,1) = Const*(1 + rho(1)^2 - rho(2)^2 - rho(3)^2);
C(1,2) = Const*(2*(rho(1)*rho(2) + rho(3)));
C(1,3) = Const*(2*(rho(1)*rho(3) - rho(2)));
C(2,1) = Const*(2*(rho(2)*rho(1) - rho(3)));
C(2,2) = Const*(1 - rho(1)^2 + rho(2)^2 - rho(3)^2);
C(2,3) = Const*(2*(rho(2)*rho(3) + rho(1)));
C(3,1) = Const*(2*(rho(3)*rho(1) + rho(2)));
C(3,2) = Const*(2*(rho(3)*rho(2) - rho(1)));
C(3,3) = Const*(1 - rho(1)^2 - rho(2)^2 + rho(3)^2);

end

