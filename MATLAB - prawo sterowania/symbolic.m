clear all;
close all;
hold on;

syms t;
syms eta1;
syms B1;
syms R1;
syms M1;
syms a;
syms b;
syms c;
syms C1 tau B M R M_r F_ry F_rx S;
syms eta2;
syms teta x_ICR;
syms tau1;
syms tau2;
syms tau3;
syms tau4;
syms r;
syms F_s1;
syms F_s2;
syms F_s3;
syms F_s4;
syms F_l1;
syms F_l2;
syms F_l3;
syms F_l4;
syms F1 F2 F3 F4;
syms q1 q2 q4 q5 q3;
syms m;

S = [cos(teta), x_ICR*sin(teta) ; sin(teta),-x_ICR*cos(teta) ; 0,1];
F_rx = cos(teta) * (F_s1 + F_s2 + F_s3 + F_s4) - sin(teta) * (F_l1 + F_l2 + F_l3 + F_l4);
F_ry = sin(teta) * (F_s1 + F_s2 + F_s3 + F_s4) + cos(teta) * (F_l1 + F_l2 + F_l3 + F_l4);

M_r = -a*(F_l1 + F_l2 + F_l3 + F_l4) + b* (F_l2 + F_l3) + c*(F_s3 + F_s4 - (F_s1 + F_s2));

R = [F_rx ; F_ry ; M_r];
M = [m,0,0 ; 0,m,0 ; 0,0,0.5];
B = (1/r).*[cos(teta), cos(teta) ; sin(teta),sin(teta) ; -c,c];

C1 = m*x_ICR*diff([0 , teta ; -teta , x_ICR]);
M1 = S.'*M*S;
R1 = S.'*R;
B1 = S.'*B;
tau = [tau1+tau2 ; tau3+tau4];




q_r = [(1+.2*sin(.9*t)) * cos(.15*t) - 1 ; (1+.2*sin(.9*t)) * sin(.15*t)]

q_r_k = diff(q_r,t)
q_r_kk = diff(q_r_k,t)








