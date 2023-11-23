% ================ TRABALHO 2 - CONTROLE DE UM CSTR ================
clc;clear all;close all;
% PARAMETROS NOMINAIS
F0 = 0.1;
T0 = 350;
c0 = 1;
ro = 1000;
Cp = 0.239;
r = 0.219;
k0 = 7.2e10;
E_R = 8750;
U = 54.94;
delta_H = -5e4;

% Pontos Linearizados - NOMINAL
x1s = 0.878;
x2s = 324.5;
x3s = 0.659;
u1s = 300;
u2s = 0.1;


% ================== Equações de Linearização ==================

% ===> NO PONTO 1

% Pontos Linearizados - OPERAÇÃO 1
x1 = 0.791;
x2 = 332.399;
x3 = 0.659;
u1 = 302.79;
u2 = 0.1;

% ===> NO PONTO 2

% Pontos Linearizados - OPERAÇÃO 2
x1 = 0.966;
x2 = 308.755;
x3 = 0.659;
u1 = 285.556;
u2 = 0.1;

x0 = [0.96 300 0.659]';
xk = x0;

% Equacoes Diferenciais

Ts = 1; %Periodo de amostragem
Kt = 0;

f1 = @(t,y)[F0*(c0-y(1))/(pi*(r^2)*y(3))-k0*y(1)*exp((-E_R)/y(2));
            F0*(T0-y(2))/(pi*(r^2)*y(3))-(delta_H*y(1)/(ro*Cp))*exp((-E_R)/y(2))+2*(U/(r*ro*Cp)*(u1-y(2)));
            F0*(T0-u2)/(pi*(r^2))];
        
[T,Z] = ode45(f1,[Kt*Ts (Kt+1)*Ts],x0)   

x0 = Z(size(Z,1),:)'
xk = x0 - [x1s x2s x3s]'






