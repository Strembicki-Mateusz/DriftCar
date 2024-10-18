clear all;
close all;
hold on;


global e a b c I x_ICR r m g ro a0 a1 e1 e2 k1 k2 k3 k4 eta;
global tau delta0;
% Variables
x_0 = 0;
y_0 = 1;
teta_0 = 0;
%eta1_0 = 0;
%eta2_0 = 0.0;
% Obliczyć (pytanie)
a0 = 1;
delta0 = ;
phi = 1/2*pi;
z_d_01 = a0*cos(phi);
z_d_02 = a0*sin(phi);
u10 = 0;
u20 = 0;
q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,0,0];
global T_max;
T_max = 2;
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
x_ICR = 0.015;
%eta1=0;
%eta2=0;


%

global u_lc;
global u_sc;

u_lc = 0.2;
u_sc = 0.02;

a1 = 0.5;
e1 = 0.01;
e2 = 0.02;
k1 = 1;
k2 = 1;
k3 = 5;
k4 = 2;

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



% Calculations
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
[t,q] = ode45(@Control, T, q_0, opts);

% Drawing

%set(gca,'FontSize','25')
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



figure(2);hold on; grid on;
plot(t, q(:,1),'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]','FontSize',20);
ylabel('X [m]','FontSize',20);

figure(3);hold on; grid on;
plot(t, q(:,2),'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]', 'FontSize',20);
ylabel('Y [m]', 'FontSize',20);


figure(4); hold on; grid on;
plot(t, q(:,4),'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]','FontSize',20);
ylabel('θ [rad]','FontSize',20);


figure(5);hold on; grid on;
plot(t, q(:,4), 'r- ', 'LineWidth',3);
plot(t, q(:,5), 'b--', 'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]','FontSize',20);
ylabel('z_d','FontSize',20);
legend('z_{d1}', 'z_{d2}', 'FontSize',16);


figure(6);hold on; grid on;
plot(t, q(:,6), 'r- ', 'LineWidth',3);
plot(t, q(:,7), 'b--', 'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]','FontSize',20);
ylabel('u','FontSize',20);
legend('u_1', 'u_2', 'FontSize',16);

figure(7);hold on; grid on;
plot(t, q(:,8), 'r- ', 'LineWidth',3);
plot(t, q(:,9), 'b--', 'LineWidth',3);
ax = gca; % current axes
ax.FontSize = 15;
%ax.TickDir = 'out';
ax.TickLength = [0.02 0.02];
get(gca)
%get(p)
xlabel('Czas [s]','FontSize',20);
ylabel('\eta','FontSize',20);
legend('ϑ_x', 'ω', 'FontSize',16);




% global q_f;
% 
% rob=1;
% 
% tau_11 = tau(1);
% tau_22 = tau(2);
% x_f = q_f(1);
% y_f = q_f(2);
% teta_f = q_f(3);
% 
% while (rob==1)
%     rob=0;
% for zlicz=2:1:length(t)
% if (t(zlicz-1)>t(zlicz))
%     temp(1)=t(zlicz-1);
%     t(zlicz-1)=t(zlicz);
%     t(zlicz)=temp(1);
% 
%     % tau
%     temp(2)=tau_11(zlicz-1);   
%     tau_11(zlicz-1)=tau_11(zlicz);
%     tau_11(zlicz)=temp(2);
% 
%     temp(2)=tau_22(zlicz-1);   
%     tau_22(zlicz-1)=tau_22(zlicz);
%     tau_22(zlicz)=temp(2);
% 
%     % q_f
%     temp(2)=x_f(zlicz-1);   
%     x_f(zlicz-1)=x_f(zlicz);
%     x_f(zlicz)=temp(2);
% 
%     temp(2)=y_f(zlicz-1);   
%     y_f(zlicz-1)=y_f(zlicz);
%     y_f(zlicz)=temp(2);
% 
%     temp(2)=teta_f(zlicz-1);   
%     teta_f(zlicz-1)=teta_f(zlicz);
%     teta_f(zlicz)=temp(2);
% 
%     rob=1;
% end
% end
% end
% 
% 
% figure(8);hold on; grid on;
% plot(t, tau_11, 'r- ', 'LineWidth',3);
% plot(t, tau_22, 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\tau','FontSize',20);
% legend('\tau_L', '\tau_P', 'FontSize',16);
% 
% 
% figure(9);hold on; grid on;
% plot(t, x_f, 'r- ', 'LineWidth',3);
% plot(t, y_f, 'b-', 'LineWidth',3);
% plot(t, teta_f, 'g-', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('q_f','FontSize',20);
% legend('e_x', 'e_y', 'e_{θ}', 'FontSize',16);
% 
% 
% 
% 
% 
% % inne phi
% x_0 = 0;
% y_0 = 1;
% teta_0 = 0;
% %eta1_0 = 0;
% %eta2_0 = 0.0;
% % Obliczyć (pytanie)
% a0 = 1;
% phi = 1/3*pi;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,0,0];
% 
% opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
% [t1,q1] = ode45(@Control, T, q_0, opts);
% 
% 
% %set(gca,'FontSize','25')
% figure(1);hold on; grid on;
% plot(q1(:,1), q1(:,2),'b- ','LineWidth',3);
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
% 
% 
% 
% 
% figure(10);hold on; grid on;
% plot(t1, q1(:,1),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('X [m]','FontSize',20);
% 
% figure(11);hold on; grid on;
% plot(t1, q1(:,2),'LineWidth',3);
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
% figure(12); hold on; grid on;
% plot(t1, q1(:,4),'LineWidth',3);
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
% figure(13);hold on; grid on;
% plot(t1, q1(:,4), 'r- ', 'LineWidth',3);
% plot(t1, q1(:,5), 'b--', 'LineWidth',3);
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
% figure(14);hold on; grid on;
% plot(t1, q1(:,6), 'r- ', 'LineWidth',3);
% plot(t1, q1(:,7), 'b--', 'LineWidth',3);
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
% figure(15);hold on; grid on;
% plot(t1, q1(:,8), 'r- ', 'LineWidth',3);
% plot(t1, q1(:,9), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\eta','FontSize',20);
% legend('ϑ_x', 'ω', 'FontSize',16);
% 
% 
% 
% 
% global q_f;
% 
% rob=1;
% 
% tau_11 = tau(1);
% tau_22 = tau(2);
% x_f = q_f(1);
% y_f = q_f(2);
% teta_f = q_f(3);
% 
% while (rob==1)
%     rob=0;
% for zlicz=2:1:length(t)
% if (t(zlicz-1)>t(zlicz))
%     temp(1)=t(zlicz-1);
%     t(zlicz-1)=t(zlicz);
%     t(zlicz)=temp(1);
% 
%     % tau
%     temp(2)=tau_11(zlicz-1);   
%     tau_11(zlicz-1)=tau_11(zlicz);
%     tau_11(zlicz)=temp(2);
% 
%     temp(2)=tau_22(zlicz-1);   
%     tau_22(zlicz-1)=tau_22(zlicz);
%     tau_22(zlicz)=temp(2);
% 
%     % q_f
%     temp(2)=x_f(zlicz-1);   
%     x_f(zlicz-1)=x_f(zlicz);
%     x_f(zlicz)=temp(2);
% 
%     temp(2)=y_f(zlicz-1);   
%     y_f(zlicz-1)=y_f(zlicz);
%     y_f(zlicz)=temp(2);
% 
%     temp(2)=teta_f(zlicz-1);   
%     teta_f(zlicz-1)=teta_f(zlicz);
%     teta_f(zlicz)=temp(2);
% 
%     rob=1;
% end
% end
% end
% 
% 
% figure(16);hold on; grid on;
% plot(t1, tau_11, 'r- ', 'LineWidth',3);
% plot(t1, tau_22, 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\tau','FontSize',20);
% legend('\tau_L', '\tau_P', 'FontSize',16);
% 
% 
% figure(17);hold on; grid on;
% plot(t1, x_f, 'r- ', 'LineWidth',3);
% plot(t1, y_f, 'b-', 'LineWidth',3);
% plot(t1, teta_f, 'g-', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('q_f','FontSize',20);
% legend('e_x', 'e_y', 'e_{θ}', 'FontSize',16);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % inne phi
% x_0 = 0;
% y_0 = 1;
% teta_0 = 0;
% %eta1_0 = 0;
% %eta2_0 = 0.0;
% % Obliczyć (pytanie)
% a0 = 1;
% phi = 2/2*pi;
% z_d_01 = a0*cos(phi);
% z_d_02 = a0*sin(phi);
% q_0 = [x_0,y_0,teta_0,z_d_01,z_d_02,u10,u20,0,0];
% 
% opts = odeset('RelTol',1e-5,'AbsTol',1e-5);
% [t2,q2] = ode45(@Control, T, q_0, opts);
% 
% 
% %set(gca,'FontSize','25')
% figure(1);hold on; grid on;
% plot(q2(:,1), q2(:,2),'b- ','LineWidth',3);
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
% 
% 
% 
% 
% figure(20);hold on; grid on;
% plot(t2, q2(:,1),'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('X [m]','FontSize',20);
% 
% figure(21);hold on; grid on;
% plot(t2, q2(:,2),'LineWidth',3);
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
% figure(22); hold on; grid on;
% plot(t2, q2(:,4),'LineWidth',3);
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
% figure(23);hold on; grid on;
% plot(t2, q2(:,4), 'r- ', 'LineWidth',3);
% plot(t2, q2(:,5), 'b--', 'LineWidth',3);
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
% figure(24);hold on; grid on;
% plot(t2, q2(:,6), 'r- ', 'LineWidth',3);
% plot(t2, q2(:,7), 'b--', 'LineWidth',3);
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
% figure(25);hold on; grid on;
% plot(t2, q2(:,8), 'r- ', 'LineWidth',3);
% plot(t2, q2(:,9), 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\eta','FontSize',20);
% legend('ϑ_x', 'ω', 'FontSize',16);
% 
% 
% 
% 
% global q_f;
% 
% rob=1;
% 
% tau_11 = tau(1);
% tau_22 = tau(2);
% x_f = q_f(1);
% y_f = q_f(2);
% teta_f = q_f(3);
% 
% 
% while (rob==1)
%     rob=0;
% for zlicz=2:1:length(t)
% if (t(zlicz-1)>t(zlicz))
%     temp(1)=t(zlicz-1);
%     t(zlicz-1)=t(zlicz);
%     t(zlicz)=temp(1);
% 
%     % tau
%     temp(2)=tau_11(zlicz-1);   
%     tau_11(zlicz-1)=tau_11(zlicz);
%     tau_11(zlicz)=temp(2);
% 
%     temp(2)=tau_22(zlicz-1);   
%     tau_22(zlicz-1)=tau_22(zlicz);
%     tau_22(zlicz)=temp(2);
% 
%     % q_f
%     temp(2)=x_f(zlicz-1);   
%     x_f(zlicz-1)=x_f(zlicz);
%     x_f(zlicz)=temp(2);
% 
%     temp(2)=y_f(zlicz-1);   
%     y_f(zlicz-1)=y_f(zlicz);
%     y_f(zlicz)=temp(2);
% 
%     temp(2)=teta_f(zlicz-1);   
%     teta_f(zlicz-1)=teta_f(zlicz);
%     teta_f(zlicz)=temp(2);
% 
%     rob=1;
% end
% end
% end
% 
% 
% figure(26);hold on; grid on;
% plot(t2, tau_11, 'r- ', 'LineWidth',3);
% plot(t2, tau_22, 'b--', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('\tau','FontSize',20);
% legend('\tau_L', '\tau_P', 'FontSize',16);
% 
% 
% figure(27);hold on; grid on;
% plot(t2, x_f, 'r- ', 'LineWidth',3);
% plot(t2, y_f, 'b-', 'LineWidth',3);
% plot(t2, teta_f, 'g-', 'LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('Czas [s]','FontSize',20);
% ylabel('q_f','FontSize',20);
% legend('e_x', 'e_y', 'e_{θ}', 'FontSize',16);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %set(gca,'FontSize','25')
% figure(30);hold on; grid on;
% plot(q(:,1), q(:,2),'r- ','LineWidth',3);
% plot(q1(:,1), q1(:,2),'b- ','LineWidth',3);
% plot(q2(:,1), q2(:,2),'g- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% legend('$\varphi = \frac{1}{2}\pi$', '$\varphi = \frac{1}{3}\pi$', '$\varphi = \pi$', 'Interpreter','latex');
% 
% figure(31);hold on; grid on;
% plot(q(:,1), q(:,2),'r- ','LineWidth',3);
% plot(q1(:,1), q1(:,2),'b- ','LineWidth',3);
% plot(q2(:,1), q2(:,2),'g- ','LineWidth',3);
% ax = gca; % current axes
% ax.FontSize = 15;
% %ax.TickDir = 'out';
% ax.TickLength = [0.02 0.02];
% get(gca)
% %get(p)
% xlabel('X [m]','FontSize',20);
% ylabel('Y [m]','FontSize',20);
% legend('φ = \pi * ^1/_2', 'φ = \pi * ^1/_3', 'φ = \pi', 'Interpreter','latex');
