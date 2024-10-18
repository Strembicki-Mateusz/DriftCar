function [ICR] = ICRDiagram(t,q)

eta1 = q(4);
eta2 = q(5);

global a b c x_ICR m g r I T_max;

eta = [eta1 ; eta2];

% Kinematics & Dynamics model
B = [1,-c ; 1,c ; 0,-x_ICR+b; 0,-x_ICR-a];
A = B * eta;

v1x = A(1,1);
v2x = A(1,1);
v3x = A(2,1);
v4x = A(2,1);
v2y = A(3,1); 
v3y = A(3,1); 
v1y = A(4,1); 
v4y = A(4,1); 

v_Y = v1y - v2y - v3y + v4y;
v_X = v1x + v2x + v3x + v4x;

% x_ICR = -v_Y/eta2;
% if x_ICR > 0.038
%     x_ICR = 0.038;
% elseif x_ICR < -0.038
%     x_ICR = -0.038;
% end

ICR = x_ICR;