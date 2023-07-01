clear all;
close all;
hold on;

% Variables
x_0 = 0;
y_0 = 0;
teta_0 = deg2rad(0);
q_0 = [x_0,y_0,teta_0];
T = [0 20];

% Controls
global u1;
global u2;
u1 = 1;
u2 = .5;




% Calculations
[t,q] = ode45(@Kinematics,T,q_0);

% Drawing
plot(q(:,1), q(:,2));
title("Diagram of the robot's movement with given controls on the XY plane");
xlabel("X");
ylabel("Y");