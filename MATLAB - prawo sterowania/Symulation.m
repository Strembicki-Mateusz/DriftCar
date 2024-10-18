clear all;
close all;
hold on;


global e a b c I x_ICR r m g ro a0 a1 e1 e2 k1 k2 k3 k4 eta;

% Variables
x_0 = 0;
y_0 = 1;
teta_0 = 0;
a0=2;
phi = -7/12*pi;
z_d_01 = a0*cos(phi);
z_d_02 = a0*sin(phi);
a1 = 0.3;
e1 = 0.01;
e2 = 0.02;
k1 = .5;
k2 = 1;
k3 = 5;
k4 = 2;
u10 = 0;
u20 = 0;
ex = x_0 - 0;
ey = y_0 - 0;
eteta = teta_0 - deg2rad(0);
q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
global T_max;
T_max = 20;
T = [0 T_max];



ro = 2;
a = 0.039;
b = 0.039;
c = 0.034;
I = 0.0036;
m = 1;
g = 9.81;
r = 0.0265;
e = 0.0001;
x_ICR = -0.015;
%eta1=0;
%eta2=0;


%

global u_lc;
global u_sc;

u_lc = 0.2;
u_sc = 0.02;



% Dynamics
global V V_0;
m_0 = 1.2;
I_0 = 0.0054;
u_s0 = .1;
u_l0 = .5;


V = [m ; I ; u_sc*m ; u_lc*m];
V_0 = [m_0 ; I_0 ; u_s0*m_0 ; u_l0*m_0];



global q_r q_r_k q_r_kk;

q_r = [0;0;0]; % Pozycja do dojazdu

% Wartości dla sterowania P2P
q_r_k = [0;0;0];
q_r_kk = [0;0;0];

% Wartości dla sterowania śledzenia trajektorii
% q_r_k = [cos(t);sin(t);0];
% q_r_kk = [-sin(t);cos(t);0];

% Calculations
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,q] = ode45(@Control_kz, T, q_0, opts);

% Drawing

figure(1);hold on; grid on;
plot(q(:,1), q(:,2),'r- ','LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('X [m]','FontSize',20);
ylabel('Y [m]','FontSize',20);
% 
% 
% figure(2);hold on; grid on;
% plot(t, q(:,1),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('X [m]','FontSize',20);
% 
% figure(3);hold on; grid on;
% plot(t, q(:,2),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]', 'FontSize',20);
% ylabel('Y [m]', 'FontSize',20);
% 
% 
% figure(4); hold on; grid on;
% plot(t, q(:,3),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('θ [rad]','FontSize',20);
% 
% 
% figure(5);hold on; grid on;
% plot(t, q(:,4), 'r- ', 'LineWidth',3);
% plot(t, q(:,5), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('z_d','FontSize',20);
% legend('z_{d1}', 'z_{d2}', 'FontSize',16);
% 
% 
% figure(6);hold on; grid on;
% plot(t, q(:,6), 'r- ', 'LineWidth',3);
% plot(t, q(:,7), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('u','FontSize',20);
% legend('u_1', 'u_2', 'FontSize',16);
% 
% figure(7);hold on; grid on;
% plot(t, q(:,8), 'r- ', 'LineWidth',3);
% plot(t, q(:,9), 'b- ', 'LineWidth',3);
% plot(t, q(:,10), 'g- ', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('Błędy [m]','FontSize',20);
% legend('e_{x}', 'e_{y}', 'e_{\theta}', 'FontSize',16);
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0 = 0;
% y_0 = 1;
% teta_0 = -pi;
% a0=3.2;
% phi = pi;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% a1 = 0.5;
% e1 = 0.01;
% e2 = 0.02;
% k1 = 1;
% k2 = 1;
% k3 = 5;
% k4 = 2;
% u10 = 0;
% u20 = 0;
% ex = x_0 - 0;
% ey = y_0 - 0;
% eteta = teta_0 - deg2rad(0);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% 
% figure(11);hold on; grid on;
% plot(q(:,1), q(:,2),'b- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% 
% figure(12);hold on; grid on;
% plot(t, q(:,1),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('X [m]','FontSize',20);
% 
% figure(13);hold on; grid on;
% plot(t, q(:,2),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]', 'FontSize',20);
% ylabel('Y [m]', 'FontSize',20);
% 
% 
% figure(14); hold on; grid on;
% plot(t, q(:,3),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('θ [rad]','FontSize',20);
% 
% 
% figure(15);hold on; grid on;
% plot(t, q(:,4), 'r- ', 'LineWidth',3);
% plot(t, q(:,5), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('z_d','FontSize',20);
% legend('z_{d1}', 'z_{d2}', 'FontSize',16);
% 
% 
% figure(16);hold on; grid on;
% plot(t, q(:,6), 'r- ', 'LineWidth',3);
% plot(t, q(:,7), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('u','FontSize',20);
% legend('u_1', 'u_2', 'FontSize',16);
% 
% figure(17);hold on; grid on;
% plot(t, q(:,8), 'r- ', 'LineWidth',3);
% plot(t, q(:,9), 'b- ', 'LineWidth',3);
% plot(t, q(:,10), 'g- ', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('Błędy [m]','FontSize',20);
% legend('e_{x}', 'e_{y}', 'e_{\theta}', 'FontSize',16);
% 
% 
% % 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0 = .5;
% y_0 = .5;
% teta_0 = -pi/2;
% a0=1.7;
% phi = pi/12;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% a1 = 0.3;
% e1 = 0.01;
% e2 = 0.02;
% k1 = 2;
% k2 = 1;
% k3 = 5;
% k4 = 2;
% u10 = 0;
% u20 = 0;
% ex = x_0 - 0;
% ey = y_0 - 0;
% eteta = teta_0 - deg2rad(0);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% 
% figure(21);hold on; grid on;
% plot(q(:,1), q(:,2),'b- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% 
% figure(22);hold on; grid on;
% plot(t, q(:,1),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('X [m]','FontSize',20);
% 
% figure(23);hold on; grid on;
% plot(t, q(:,2),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]', 'FontSize',20);
% ylabel('Y [m]', 'FontSize',20);
% 
% 
% figure(24); hold on; grid on;
% plot(t, q(:,4),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('θ [rad]','FontSize',20);
% 
% 
% figure(25);hold on; grid on;
% plot(t, q(:,4), 'r- ', 'LineWidth',3);
% plot(t, q(:,5), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('z_d','FontSize',20);
% legend('z_{d1}', 'z_{d2}', 'FontSize',16);
% 
% 
% figure(26);hold on; grid on;
% plot(t, q(:,6), 'r- ', 'LineWidth',3);
% plot(t, q(:,7), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('u','FontSize',20);
% legend('u_1', 'u_2', 'FontSize',16);
% 
% figure(27);hold on; grid on;
% plot(t, q(:,8), 'r- ', 'LineWidth',3);
% plot(t, q(:,7), 'b- ', 'LineWidth',3);
% plot(t, q(:,9), 'b- ', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('Błędy [m]','FontSize',20);
% legend('e_{x}', 'e_{y}', 'e_{\theta}', 'FontSize',16);














% 
% a0=1;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'g- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% a0=1.25;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'b--','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% a0=1.5;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'b- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);

a0=2;
z_d_01 = a0*cos(phi);
z_d_02 = a0*sin(phi);
q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'g--','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% a0 = 2;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'m- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% a0=2.25;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'k- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% a0=2.5;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'c- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% a0=3;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% 
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'b-','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% a0=4;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,ex,ey,eteta];
% opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,q] = ode45(@Control_kz, T, q_0, opts);
% 
% % Drawing
% figure(1);hold on; grid on;
% plot(q(:,1), q(:,2),'g-','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% 
% 
% 
% figure(1);
% legend('\delta_d(0) = 1,5','\delta_d(0) = 1,75','\delta_d(0) = 2','\delta_d(0) = 2,25','\delta_d(0) = 2,5', '\delta_d(0) = 3','\delta_d(0) = 4', 'FontSize',20);