clear all
clc
close all

Fo=0.1;
To=350;
Co=1;
ro=1000;
Cp=0.239;
r=0.219;
Ko=7.2*10^(10);
ER=8750;
U=54.94;
dH=-5*10^4;

% alfa=0.45;
% beta=8;
% teta=53*(pi/180);


%Operação Nominal  %Para usar no Simulink

x1s=0.878;
x2s=324.5;
x3s=0.659;
u1s=300;
u2s=0.1;

x11=(-Fo/(pi*r^2*x3s))-Ko*(exp(-ER/x2s));
x12=-(Ko*x1s*ER*exp(-ER/x2s))/(x2s^2);
x13=-(Fo*Co-Fo*x1s)/(pi*r^2*x3s^2);

x21=(-dH*Ko*exp(-ER/x2s))/(ro*Cp);
x22=(-Fo/(pi*r^2*x3s))-2*U/(r*ro*Cp)-(dH*Ko*x1s*ER*exp(-ER/x2s))/(ro*Cp*x2s^2);
x23=-(Fo*To-Fo*x2s)/(pi*r^2*x3s^2);

x31=0;
x32=0;
x33=0;

Anominal=[x11 x12 x13;
   x21 x22 x23;
   x31 x32 x33];

B221=2*U/(r*ro*Cp);
B232=-1/(pi*r^2);

B2nominal=[0 0;
   B221 0;
   0 B232];

B111=(Co-x1s)/(pi*r^2*x3s);
B121=(To-x2s)/(pi*r^2*x3s);
B131=1/(pi*r^2);

B1nominal=[B111;B121;B131];

C2nominal=[1 0 0;0 0 1];

D11nominal=zeros(2,1);
D12nominal=zeros(2,2);

%Discretização

A=Anominal;
B=B2nominal;
C=C2nominal;
D=D12nominal;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 

B2d=sd.B;
D12d=sd.D;



 
A=Anominal;
B=B1nominal;
C=C2nominal;
D=D11nominal;

 sys=ss(A,B,C,D);
 T=1; 
 sd = c2d (sys,T,'zoh'); 
 
Ad=sd.A;
B1d=sd.B;
C2d=sd.C;
D11d=sd.D;

disp('Matrizes Discretizadas Operação Nominal:')

Ad
B1d
B2d
C2d
D11d
D12d


disp('Solução EDO Operação Nominal:')

u1=.1;
u2=.2;

cstr=@(t,z)[(Fo*Co/(pi*r^2*z(3)))-(Fo*z(1)/(pi*r^2*z(3)))-(Ko*z(1)*exp(-ER/z(2)));
    (Fo*To/(pi*r^2*z(3)))-(Fo*z(2)/(pi*r^2*z(3)))-((dH/(ro*Cp))*z(1)*exp(-ER/z(2)))+((2*U/(ro*Cp))*(u1-z(2)));
    Fo/(pi*r^2)-Fo*u2/(pi*r^2)]

[T,Z]=ode45 (cstr,[0 1],[0.6 300 0.5]')


%return

 

% Operação 1

x1s=0.791;
x2s=332.399;
x3s=0.659;
u1s=302.79;
u2s=0.1;

%Eq1
%Derivadas parciais em relação a x1
x11=(-Fo/(pi*r^2*x3s))-Ko*(exp(-ER/x2s)); %d/dx1
%Derivadas parciais em relação a x2
x12=-(Ko*x1s*ER*exp(-ER/x2s))/(x2s^2); %d/dx2
%Derivadas parciais em relação a x3
x13=-(Fo*Co-Fo*x1s)/(pi*r^2*x3s^2); %d/dx3

%Eq 2
%Derivadas parciais em relação a x1
x21=(-dH*Ko*exp(-ER/x2s))/(ro*Cp);%d/dx1
%Derivadas parciais em relação a x1
x22=(-Fo/(pi*r^2*x3s))-2*U/(r*ro*Cp)-(dH*Ko*x1s*ER*exp(-ER/x2s))/(ro*Cp*x2s^2); %d/dx2
%Derivadas parciais em relação a x1
x23=-(Fo*To-Fo*x2s)/(pi*r^2*x3s^2); %d/dx3

%Eq3
x31=0;
x32=0;
x33=0;

A1=[x11 x12 x13;
   x21 x22 x23;
   x31 x32 x33];

%Derivadas parciais em relação a u1
B221=2*U/(r*ro*Cp);
%Derivadas parciais em relação a u2
B232=-1/(pi*r^2);
%O resto é zero
B2=[0 0;
   B221 0;
   0 B232];

%Derivadas parciais em relação a w1 = Fo das Eq 1, 2 e 3
B111=(Co-x1s)/(pi*r^2*x3s);
B121=(To-x2s)/(pi*r^2*x3s);
B131=1/(pi*r^2);

B11=[B111;B121;B131];


C2=[1 0 0;0 0 1];

D21=zeros(2,1);
D22=zeros(2,2);


C1 = C2;
D11 = D21;
D12=D22;

%Discretização 

A=A1;
B=B2;
% B=B2nominal;
C=C2;
D=D12;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 

B2d=sd.B;
D12d=sd.D;

A=A1;
B=B11;
C=C2;
D=D11;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 
 
A1d=sd.A;
B11d=sd.B;
C2d=sd.C;
D11d=sd.D;

disp('Matrizes Discretizadas Ponto de Operação 1:')

A1d
B11d
B2d
C2d
D11d
D12d



%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Operação 2 

x1s=0.966;
x2s=308.755;
x3s=0.659;
u1s=285.556;
u2s=0.1;

x11=(-Fo/(pi*r^2*x3s))-Ko*(exp(-ER/x2s));
x12=-(Ko*x1s*ER*exp(-ER/x2s))/(x2s^2);
x13=-(Fo*Co-Fo*x1s)/(pi*r^2*x3s^2);

x21=(-dH*Ko*exp(-ER/x2s))/(ro*Cp);
x22=(-Fo/(pi*r^2*x3s))-2*U/(r*ro*Cp)-(dH*Ko*x1s*ER*exp(-ER/x2s))/(ro*Cp*x2s^2);
x23=-(Fo*To-Fo*x2s)/(pi*r^2*x3s^2);

x31=0;
x32=0;
x33=0;

A2=[x11 x12 x13;
   x21 x22 x23;
   x31 x32 x33];

B221=2*U/(r*ro*Cp);
B232=-1/(pi*r^2);

B2=[0 0;
   B221 0;
   0 B232];

B111=(Co-x1s)/(pi*r^2*x3s);
B121=(To-x2s)/(pi*r^2*x3s);
B131=1/(pi*r^2);

B12=[B111;B121;B131];

C2=[1 0 0;0 0 1];

D21=zeros(2,1);
D22=zeros(2,2);


% C1 = [C2; zeros(2,3)];
% D11 = [D21 ; zeros(2,1)];
% D12=[D22 0.04*eye(2)]';

C1 = C2;
D11 = D21;
D12=D22;

%Discretização 

A=A2;
B=B2;
% B=B2nominal;
C=C2;
D=D12;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 

B2d=sd.B;
D12d=sd.D;

A=A2;
B=B12;
C=C2;
D=D11;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 
 
A2d=sd.A;
B12d=sd.B;
C2d=sd.C;
D11d=sd.D;

disp('Matrizes Discretizadas Ponto de Operação 2:')

A2d
B12d
B2d
C2d
D11d
D12d