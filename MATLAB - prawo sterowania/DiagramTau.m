function [tau] = DiagramTau(t)

global tau1 tau2 tau3 tau4;
global T_max;
global wzm1 wzm2 wzm3 wzm4;

    tau1 = trapez(wzm1,t,T_max);
    tau2 = trapez(wzm2,t,T_max);
    tau3 = trapez(wzm3,t,T_max);
    tau4 = trapez(wzm4,t,T_max);
    


tau = [tau1+tau2,tau3+tau4];