%Written by Neil McHenry
% February 13th 2018

% 3.MATLAB Code Suite moving from/to the following attitude
% parameterizations:
%
%   a) Direction Cosine Matrix (C):
%   b) quaternion (q):


function [q, flag] = Ctoq(C)

%Calculate quaternion from DCM

q4 = (1/2)*(1+trace(C))^(1/2);

qv = (1/(4*q4))*[C(2,3) - C(3,2)
                 C(3,1) - C(1,3)
                 C(1,2) - C(2,1)];
q = [qv;q4];

flag = 'Success';

