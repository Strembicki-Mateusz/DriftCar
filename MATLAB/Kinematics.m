function [q1] = Kinematics(t,q)

% Variables
global u1;
global u2;
u = [u1;u2];
x = q(1);
y = q(2);
teta = q(3);


% Kinematics model
S = [cos(teta),0 ; sin(teta),0 ; 0,1];
q1 = S*u;

q1(3)
