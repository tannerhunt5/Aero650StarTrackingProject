%Written by Neil McHenry
% February 13th 2018

% 3.MATLAB Code Suite moving from/to the following attitude
% parameterizations:
%
%   a) Direction Cosine Matrix (C):
%   d) Gibbs/Rodrigues vector (rho):

function [C, flag] = rhotoC(rho)

%define the Identity matrix
I = eye(3);
G = [ 0,      -rho(3),  rho(2);
      rho(3),       0, -rho(1);
     -rho(2),  rho(1),      0];

C = ((1-(rho.'*rho))*I + 2*(rho*rho.') - 2*G)/( 1 + (rho.'*rho));

flag = 'success';

