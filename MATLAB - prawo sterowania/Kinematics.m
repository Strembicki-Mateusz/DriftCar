function [final] = Kinematics(t,q)

x = q(1);
y = q(2);
teta = q(3);
eta1 = q(4);
eta2 = q(5);
z_d_0 = q(6);


global tau1 tau2 tau3 tau4;

global a b c x_ICR m g r I T_max e;

global u_lc u_sc;

tau1 = trapez(0.05,t,T_max);
tau2 = trapez(0.02,t,T_max);
tau3 = trapez(0.02,t,T_max);
tau4 = trapez(0.05,t,T_max);



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

x_ICR = -v_Y*eta2/(eta2^2 + e);
% if x_ICR > 0.039
%     x_ICR = 0.039;
% elseif x_ICR < -0.039
%     x_ICR = -0.039;
% end

F_l1 = u_lc * m * g * sgn(v1y);
F_l2 = u_lc * m * g * sgn(v2y);
F_l3 = u_lc * m * g * sgn(v3y);
F_l4 = u_lc * m * g * sgn(v4y);

F_s1 = u_sc * m * g * sgn(v1x);
F_s2 = u_sc * m * g * sgn(v2x);
F_s3 = u_sc * m * g * sgn(v3x);
F_s4 = u_sc * m * g * sgn(v4x);




S = [cos(teta),  x_ICR*sin(teta);sin(teta), -x_ICR*cos(teta); 0,1];

%F_rx = cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4);
%F_ry = cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4);
%M_r = b*(F_l2 + F_l3) - a*(F_l1 + F_l2 + F_l3 + F_l4) - c*(F_s1 + F_s2 - F_s3 - F_s4);
%R = [cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4) ; cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4) ; b*(F_l2 + F_l3) - a*(F_l1 + F_l2 + F_l3 + F_l4) - c*(F_s1 + F_s2 - F_s3 - F_s4)];
%M =[m, 0,   0;0, m,   0;0, 0, I];
%B = [cos(teta)/r, cos(teta)/r;sin(teta)/r, sin(teta)/r;       -c/r,         c/r];


C1 = [0,       m*x_ICR*eta2 ; -m*x_ICR* eta2, m*0];     %x_ICR]
M1 = [m*cos(teta)^2 + m*sin(teta)^2,0;0, m*x_ICR^2*cos(teta)^2 + m*x_ICR^2*sin(teta)^2 + I];
R1 = [cos(teta)*(cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4)) + sin(teta)*(cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4));b*(F_l2 + F_l3) - a*(F_l1 + F_l2 + F_l3 + F_l4) - c*(F_s1 + F_s2 - F_s3 - F_s4) - x_ICR*cos(teta)*(cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4)) + x_ICR*sin(teta)*(cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4))];
B1 = [cos(teta)^2/r + sin(teta)^2/r, cos(teta)^2/r + sin(teta)^2/r;-c/r,c/r];
tau =[tau1 + tau2;tau3 + tau4];

q2 = inv(M1)*(B1*tau - (C1*eta + R1));
q1 = S*eta;

final = [q1;q2];









