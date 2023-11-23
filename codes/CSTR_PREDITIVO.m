% ============== TRABALHO 3 - CONTROLE PREDITIVO DE UM CSTR ===============
clc;clear all;close all;

% PARÂMETROS NOMINAIS
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

% Pontos de Linearização NOMINAL e REFERÊNCIA DO SISTEMA
x1s = 0.878;
x2s = 324.5;
x3s = 0.659;
u1s = 300;
u2s = 0.1;

% ================== Equações de Linearização ==================

% ===> NO PONTO 1

% Pontos de Linearização - PONTO DE OPERAÇÃO 1
x1 = 0.791;
x2 = 332.399;
x3 = 0.659;
u1 = 302.79;
u2 = 0.1;

% MATRIZ A
h1_x1 = -F0/(pi*r^2*x3) - k0*exp(-E_R/x2);
h1_x2 = -(k0*x1*E_R/(x2^2))*exp(-E_R/x2);
h1_x3 = -F0*(c0-x1)/(pi*r^2*x3^2);
h2_x1 = ((-delta_H*k0)/(ro*Cp))*(exp(-E_R/x2));
h2_x2 = (-F0)/(pi*r^2*x3) - (((delta_H*k0*x1*E_R)/(ro*Cp*x2^2)))*(exp(-E_R/x2)) - (2*U)/(r*ro*Cp);
h2_x3 = -F0*(T0-x2)/(pi*r^2*x3^2);
h3 = 0;

A1 = [h1_x1 h1_x2 h1_x3;h2_x1 h2_x2 h2_x3;h3 h3 h3];

% MATRIZ B2

h2_u1 = (2*U)/(r*ro*Cp);
h3_u2 = -1/(pi*r^2);

B2 = [0 0;h2_u1 0;0 h3_u2];

% MATRIZ C
C2 = [1 0 0;0 0 1];

% MATRIZ D
D22 = zeros([size(C2,1) size(B2,2)]);


% ===> NO PONTO 2

% Pontos de Linearização - PONTO DE OPERAÇÃO 2
x1 = 0.966;
x2 = 308.755;
x3 = 0.659;
u1 = 285.556;
u2 = 0.1;

% MATRIZ A
h1_x1 = -F0/(pi*r^2*x3) - k0*exp(-E_R/x2);
h1_x2 = -(k0*x1*E_R/(x2^2))*exp(-E_R/x2);
h1_x3 = -F0*(c0-x1)/(pi*r^2*x3^2);
h2_x1 = ((-delta_H*k0)/(ro*Cp))*(exp(-E_R/x2));
h2_x2 = (-F0)/(pi*r^2*x3) - (((delta_H*k0*x1*E_R)/(ro*Cp*x2^2)))*(exp(-E_R/x2)) - (2*U)/(r*ro*Cp);
h2_x3 = -F0*(T0-x2)/(pi*r^2*x3^2);
h3 = 0;

A2 = [h1_x1 h1_x2 h1_x3;h2_x1 h2_x2 h2_x3;h3 h3 h3];

% ===> Discretização
Ts = 1; %Tempo de amostragem de 1 min

% Sistema no ponto de Operação 1
ssG1 = ss(A1,B2,C2,D22);
ssG1_d = c2d(ssG1,Ts,'zoh');
[A1 B21 C2 D22] = ssdata(ssG1_d);

% Sistema no ponto de Operação 2
ssG2 = ss(A2,B2,C2,D22);
ssG2_d = c2d(ssG2,Ts,'zoh');
[A2 B22 C2 D22] = ssdata(ssG2_d);


% ============= CONTROLE PREDITIVO COM INCERTEZAS DE POLITOPO =============

% Matrizes de Ponderação
Q = (C2')*C2;
R = 100000*eye(2);
Q_2 = sqrt(Q);
R_2 = sqrt(R);
x = []';

%Vértices do Politopo
A = [];
A(:,:,1) = A1;
A(:,:,2) = A2;    

B2 = [];
B2(:,:,1) = B21;
B2(:,:,2) = B22;

% Dimensões das Matrizes
[n,p,ra] = size(A);
[p,m,rb] = size(B2);

% CONDIÇÃO INICIAL DO SISTEMA
x0 = [0.96 300 0.659]';

% Diferenças por conta da linearização
xk = x0 - [x1s x2s x3s]';
q = 1;

% Tempo de Simulação
tempo = 100; 

for Kt = 0:tempo

% ===> PREDITIVO por LMI  
setlmis([]);

% ===> Declaração das variaveis
gama = lmivar(1,[1 1]);
[W1,n_W1] = lmivar(1,[n 1]);
[W2,n_W2] = lmivar(2,[m,n]);

% ===> Declaração das LMIs com for()

%Sujeição do Preditivo
 
lmiterm([-q 1 1 0],1);
lmiterm([-q 2 1 0],xk); 
lmiterm([-q 2 2 W1],1,1);

%Declaração Principal
q = q+1;
for i = 1:ra  
    for ii = 1:rb
       lmiterm([-q 1 1 W1],1,1);
       lmiterm([-q 2 1 W1],A(:,:,i),1);
       lmiterm([-q 2 1 W2],B2(:,:,ii),1);
       lmiterm([-q 3 1 W1],Q_2,1);
       lmiterm([-q 4 1 W2],R_2,1);       
       lmiterm([-q 2 2 W1],1,1);
       lmiterm([-q 3 2 0],0);
       lmiterm([-q 4 2 0],0);
       lmiterm([-q 3 3 gama],1,1);
       lmiterm([-q 4 3 0],0);
       lmiterm([-q 4 4 gama],1,1);       
       q = q+1;     
    end   
end

%Declaração de Gama>0
q = q + 1;
lmiterm([-q 1 1 gama],1,1);

%Declaração de X>0
q = q + 1;
lmiterm([-q 1 1 W1],1,1);

% Fim das declarações
preditivo = getlmis;

%Projeto com mincx
c = [1 zeros(1,n_W2-1)]';
options = [1e-8,3000,1e9,1e-4,1];

[c_opt,x_opt] = mincx(preditivo,c,options);

gama = dec2mat(preditivo,x_opt,gama);
W1_opt = dec2mat(preditivo,x_opt,W1);
W2_opt = dec2mat(preditivo,x_opt,W2);

%O Ganho do Controlador
K = W2_opt*inv(W1_opt);

% Sinal de Controle
if Kt == 0    
    u = K*xk + [u1s u2s]';
    u0 = u;
else
    u = K*xk + [u1s u2s]';    
end

% Resolução do CSTR em Equações Diferenciais
f1 = @(t,y)[(F0*(c0-y(1))/(pi*(r^2)*y(3)))-k0*y(1)*exp((-E_R)/y(2));
            (F0*(T0-y(2))/(pi*(r^2)*y(3)))-(delta_H*k0*y(1)/(ro*Cp))*exp((-E_R)/y(2))+2*(U/(r*ro*Cp)*(u(1,:)-y(2)));
            (F0-u(2,:))/(pi*(r^2))];
        
[T,Y] = ode45(f1,[Kt*Ts (Kt+1)*Ts],x0');   

% Criação do matriz y que contém todos valores de Y
if Kt ==0
   y = Y;
else
   y = vertcat(y,Y);    
end    

% Estados Previstos na estratégia de horizonte deslizante
x(:,Kt+1) = Y(size(Y,1),:)';
x0 = x(:,Kt+1);
xk = x(:,Kt+1) - [x1s x2s x3s]';
end

% Tempo total de simulação
t = linspace(0,tempo,(tempo + 1)*size(Y,1)); % Tamanho do vetor Y ao longo do tempo

% ===============> PLOTAGEM DE GRAFICOS <===============

% SAÍDA DE CONCENTRAÇÃO - Estado X1
figure;set(gcf,'color','w');
plot(t,y(:,1), t, ones(4141,1)*0.878, 'r');
ylim([0.8776 0.879])
% ESTE COMANDO YLINE() SÓ ESTA DISPONÍVEL NAS VERSOES A PARTIR DA VERSÃO 2018 DO MATLAB
%yline(0,0.878,'--r','LineWidth',2);
legend('Saída de Concentração','Concentração de 0,878')
grid minor
title("SAÍDA DE CONCENTRAÇÃO - Estado X1",'FontSize',20)
xlabel("Tempo (min)",'FontSize',18)
ylabel("Concentração (k.mol/m3)",'FontSize',18)

% TEMPERATURA - Estado X2
figure;set(gcf,'color','w');
plot(t,y(:,2), t, ones(4141,1)*324.5, 'r');
ylim([324.4 324.65])
grid minor
% ESTE COMANDO YLINE() SÃ“ ESTA DISPONIVEL NAS VERSOES A PARTIR DA 2018 DO MATLAB
%yline(324.5,'--r','LineWidth',2);
legend('Estado X2 - Temperatura','Temperatura de 324,5 K')
title("TEMPERATURA - Estado X2", 'FontSize',20)
xlabel("Tempo (min)")
ylabel("Temperatura (K)")

% SAIDA DE NÍVEL - Estado X3
figure;set(gcf,'color','w');
plot(t,y(:,3), t, ones(4141,1)*0.659, 'r');
ylim([0.6589 0.6591])
grid minor
legend('Saída de Nível', "Nível de 0,659m")
title("SAÍDA DE NÍVEL - Estado X3",'FontSize',20)
xlabel("Tempo (min)",'FontSize',18)
ylabel("Nível (m)",'FontSize',18)

return
