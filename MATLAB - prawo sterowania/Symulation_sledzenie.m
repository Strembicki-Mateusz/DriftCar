clear all;
close all;
hold on;

global iteratorDisp

iteratorDisp = 0;

global e a b c I x_ICR r m g ro a0 a1 e1 e2 k1 k2 k3 k4 tau_vector t_vector;

% Variables
x_0 = .5;
y_0 = .5;
%teta_0 = -pi/2;
teta_0 = -pi/2;
%eta1_0 = 0;
%eta2_0 = 0.0;
% Obliczyć (pytanie)
%a0 = 1.5;
a0 = 3.2;
e1 = 0.5;
z_d_01 = (a0+e1)*cos(-1*pi/12);
z_d_02 = (a0+e1)*sin(-1*pi/12);

%z_d_01 = a0*cos(1/2*pi);
%z_d_02 = a0*sin(1/2*pi);
u10 = 0;
u20 = 0;
ex = x_0 - 0;
ey = y_0 - 0;
eteta = teta_0 - deg2rad(0);
q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,0,ex,ey,eteta,0,0];
global T_max;
T_max = 100;
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

a1 = 0.3;
e1 = 0.1;
%DP
e2 = 0.05;
k1 = 2;
k2 = 1;
k3 = 80;
k4 = 50;

% Dynamics
global V V_0;
m_0 = 1.2;
I_0 = 0.0054;
u_s0 = .1;
u_l0 = .5;


V = [m ; I ; u_sc*m ; u_lc*m];
V_0 = [m_0 ; I_0 ; u_s0*m_0 ; u_l0*m_0];



global q_r q_r_k q_r_kk;

%q_r = [0;0;0]; % Pozycja do dojazdu

% Wartości dla sterowania P2P
%q_r_k = [0;0;0];
%q_r_kk = [0;0;0];

% Wartości dla sterowania śledzenia trajektorii
% q_r_k = [cos(t);sin(t);0];
% q_r_kk = [-sin(t);cos(t);0];

% Calculations
opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,q] = ode23(@Control_sledzenie, T, q_0, opts);

% Drawing

figure(1);hold on; grid on;
plot(q(:,1), q(:,2),'r-- ','LineWidth',3);
plot(q(:,12), q(:,13) , 'g--','LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('X [m]','FontSize',20);
ylabel('Y [m]','FontSize',20);
legend(' Zakreślona \newline trajektoria \newline przez robota', 'Trajektoria', 'FontSize',12);
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
% plot(t, q(:,9), 'r- ', 'LineWidth',3);
% plot(t, q(:,10), 'b- ', 'LineWidth',3);
% plot(t, q(:,11), 'g- ', 'LineWidth',3);
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
% figure(8);hold on; grid on;
% plot(t_vector, tau_vector(:,1), 'r- ', 'LineWidth',3);
% plot(t_vector, tau_vector(:,2), 'g--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\tau [Nm]','FontSize',20);
% legend("\tau_L = \tau_1 + \tau_2", "\tau_P = \tau_3 + \tau_4", 'FontSize',16);