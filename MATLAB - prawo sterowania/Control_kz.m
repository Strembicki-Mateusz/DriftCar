function [final] = Control(t,q)

x = q(1);
y = q(2);
teta = q(3);
%eta1 = q(4);
%eta2 = q(5);
z_d1 = q(4);
z_d2 = q(5);
u1 = q(6);
u2 = q(7);
q_f = [q(8);q(9);q(10)];

qS = [x;y;teta];


global a b c I x_ICR r m g ro a0 a1 e1 e2 k1 k2 k3 eta;
global V_0;
global u_lc u_sc;

%eta = [eta1 ; eta2];
z_d = [z_d1 ; z_d2];
u = [u1;u2];

if t<1e-10
    eta=[0;0];
end

eta_disp=eta;

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

%v_Y = v1y - v2y - v3y + v4y;
%v_X = v1x + v2x + v3x + v4x;

F_l1 = u_lc * m * g * sgn(v1y);
F_l2 = u_lc * m * g * sgn(v2y);
F_l3 = u_lc * m * g * sgn(v3y);
F_l4 = u_lc * m * g * sgn(v4y);

F_s1 = u_sc * m * g * sgn(v1x);
F_s2 = u_sc * m * g * sgn(v2x);
F_s3 = u_sc * m * g * sgn(v3x);
F_s4 = u_sc * m * g * sgn(v4x);

% x_ICR = -v_Y*eta2/(eta2^2 + e);




C1 = [0,       m*x_ICR*eta(2); -m*x_ICR* eta(2), m*0];     %x_ICR]
M1 = [m,0;0, m*x_ICR^2 + I];
R1 = [cos(teta)*(cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4)) + sin(teta)*(cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4));b*(F_l2 + F_l3) - a*(F_l1 + F_l2 + F_l3 + F_l4) - c*(F_s1 + F_s2 - F_s3 - F_s4) - x_ICR*cos(teta)*(cos(teta)*(F_l1 + F_l2 + F_l3 + F_l4) + sin(teta)*(F_s1 + F_s2 + F_s3 + F_s4)) + x_ICR*sin(teta)*(cos(teta)*(F_s1 + F_s2 + F_s3 + F_s4) - sin(teta)*(F_l1 + F_l2 + F_l3 + F_l4))];
B1 = [1/r, 1/r;-c/r,c/r];
%tau =[tau1 + tau2;tau3 + tau4];
S = [cos(teta),  x_ICR*sin(teta);sin(teta), -x_ICR*cos(teta); 0,1];
%eta = [eta1 ; eta2];
q1 = S*eta;


global q_r q_r_k q_r_kk;


x_r = q_r(1);
y_r = q_r(2);
teta_r = q_r(3);
q1_r = [q_r_k(1); q_r_k(2); q_r_k(3)];
q1_r_k = [q_r_kk(1); q_r_kk(2); q_r_kk(3)];

q_f = qS - q_r;

S_2_2_r = [cos(teta_r),  x_ICR*sin(teta_r);sin(teta_r), -x_ICR*cos(teta_r)];
eta_r = inv(S_2_2_r)*[q1_r(1);q1_r(2)];

%S_r = [cos(teta_r),  x_ICR*sin(teta_r);sin(teta_r), -x_ICR*cos(teta_r); 0,1];
%q1_r_liczone = S_r * eta_r; % q ref. 1 pochodna


S_2_2_k_r = [-sin(teta_r),cos(teta_r);(1/x_ICR)*cos(teta_r), (1/x_ICR)*sin(teta_r)];
%S_2_2_k_r = [-sin(teta_r), x_ICR*cos(teta_r) ; cos(teta_r), x_ICR*sin(teta_r)]*eta_r(2);



eta_r_k = S_2_2_k_r*[q1_r(1);q1_r(2)]*eta_r(2) + inv(S_2_2_r)*[q1_r_k(1);q1_r_k(2)];


q_f_k = q1 - q1_r;

J=[0,-1;1,0];

z1 = q_f(3);
z2 = q_f(1)*cos(q(3)) + q_f(2)*sin(q(3));
z = [z1;z2];
w = -z(1)*z(2)+2*(q_f(1)*sin(q(3)) - q_f(2)*cos(q(3)) - x_ICR*q_f(3));
z1_k = q_f_k(3);
z2_k = q_f_k(1)*cos(q(3)) - q_f(1)*sin(q(3))*q1(3) + q_f_k(2)*sin(q(3)) + q_f(2)*cos(q(3))*q1(3);
z_k = [z1_k;z2_k];


%Z = [w ; z1 ; z2];

f = 2*(-eta_r(1)*sin(z(1)) + eta_r(2)*(x_ICR + z(2) - x_ICR*cos(z(1))));
f_k = 2*(-eta_r_k(1)*sin(z(1)) - eta_r(1)*cos(z(1))*z_k(1) + eta_r_k(2)*x_ICR + eta_r_k(2)*z(2)+eta_r(2)*z_k(2) + x_ICR*(-eta_r_k(2)*cos(z(1)) + eta_r(2)*sin(z(1))*z_k(1))); 
w_k = z_k(1)*z(2)-z(1)*z_k(2) + f;



l = q_f(1)*sin(q(3)) - q_f(2)*cos(q(3));
l_k = q_f_k(1)*sin(q(3)) + q_f(1)*q1(3)*cos(q(3)) - (q_f_k(2)*cos(q(3)) - sin(q(3))*q1(3)*q_f(2));
T = [l,1;1,0];
%T__1 = inv(T);
T_k = [l_k,0;0,0];
II = [eta_r(1)*cos(z(1)) + eta_r(2)*(-x_ICR*sin(z(1))+l); eta_r(2)];
II_k = [eta_r_k(1)*cos(z(1)) - eta_r(1)*sin(z(1))*z_k(1) + eta_r_k(2)*l + eta_r(2)*l_k - x_ICR*(eta_r_k(2)*sin(z(1)) + eta_r(2)*cos(z(1))*z_k(1)) ; eta_r_k(2)];






delta_d = a0 * exp(-a1 * t) + e1;
delta_d_k = a0 * -a1 * exp(-a1*t);
delta_d_2k = a0 * a1^2 * exp(-a1*t);
Omega1 = k2 + delta_d_k/delta_d + w*(k1*w+f)/(delta_d^2); 
Omega1_k = (delta_d_2k*delta_d - delta_d_k^2)/(delta_d^2) + k1*((2*w*w_k*delta_d^2 - w^2*2*delta_d*delta_d_k)/(delta_d^4)) + (f_k * delta_d^2 - f * 2*delta_d*delta_d_k)/(delta_d^4); 


z_d_k = (delta_d_k/delta_d) * z_d+((k1*w+f)/(delta_d^2) +w*Omega1) * J * z_d;


u_a_k = (k1*w+f)/(delta_d^2) * J * z_d + Omega1 * z_d;
u_a_k1 = (((k1*w_k+f_k)*delta_d^2 - (k1*w+f)*2*delta_d*delta_d_k)/(delta_d^4))*J *z_d + ((k1*w+f)/delta_d^2)*J*z_d_k+Omega1_k*z_d+Omega1*z_d_k;


u_d = u_a_k - k2*z;
%u = T__1*(eta-II);
u_d_k = u_a_k1-k2*z_k;
u_f = u_d - u;


%dobrać alfa_0 (symulacje)

C2 = T.'*(C1*T+M1*T_k); 
M2 = T.'*M1*T;
R2 = T.'*(C1*II + M1*II_k+R1);
B2 = T.'*B1;

SUMx = g*(sgn(v1x)+sgn(v2x)+sgn(v3x)+sgn(v4x));
SUMy = g*(sgn(v1y)+sgn(v2y)+sgn(v3y)+sgn(v4y));

A = [0, q1(3) ; -q1(3),0];
Y1 = T.'*([1,0;0,x_ICR^2]*(T*u_d_k + T_k*u_d + II_k) + x_ICR*A*(T*u_d+II));
Y2 = T.'*[0,0;0,1]*(T*u_d_k+T_k*u_d + II_k);
%Y2 = T.'*[0,0;0,1]*(T*(u_d_k+u_d) + II_k);
Y3 = T.'*[cos(teta) * SUMx ; x_ICR*sin(teta)*SUMx+c*((g*sgn(v3x))+(g*sgn(v4x)) - (g*sgn(v1x))+(g*sgn(v2x)))];
Y4 = T.'*[-sin(teta) * SUMy ; x_ICR*cos(teta)*SUMy-a*((g*sgn(v1y))+(g*sgn(v4y))) + b*((g*sgn(v3y))+(g*sgn(v2y)))];
Y_d = [Y1,Y2,Y3,Y4];

tau_a = Y_d * (ro^2*Y_d.'*u_f)/(norm(Y_d.'*u_f)*ro+e2);
z_f = z_d - z;

tau = inv(B2)*(w*J*z + z_f + Y_d*V_0 + tau_a + k3*u_f);
u_k = inv(M2)*(B2*tau - (C2*u + R2));
eta = T*u+II;
q1 = S*eta;


final = [q1;z_d_k;u_k;q_f_k];