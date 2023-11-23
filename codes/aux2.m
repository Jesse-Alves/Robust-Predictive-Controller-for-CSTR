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

%Opera��o Nominal  %Para usar no Simulink

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

%Discretiza��o

A=Anominal;
B=B2nominal;
C=C2nominal;
D=D12nominal;

sys=ss(A,B,C,D);
T=1;
Ts=T;
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

disp('Matrizes Discretizadas Opera��o Nominal:')

Ad
B1d
B2d
C2d
D11d
D12d

Adn=Ad;
B1dn=B1d;
B2dn=B2d;

% disp('Solu��o EDO Opera��o Nominal:')
% 
% u1=.1;
% u2=.2;
% 
% fcstr=@(t,z)[(Fo*Co/(pi*r^2*z(3)))-(Fo*z(1)/(pi*r^2*z(3)))-(Ko*z(1)*exp(-ER/z(2)));
%     (Fo*To/(pi*r^2*z(3)))-(Fo*z(2)/(pi*r^2*z(3)))-((dH/(ro*Cp))*z(1)*exp(-ER/z(2)))+((2*U/(ro*Cp))*(u1-z(2)));
%     Fo/(pi*r^2)-Fo*u2/(pi*r^2)]
% 
% [T,Z]=ode45 (fcstr,[0 1],[0.96 300 0.659]')


% return

 

% Opera��o 1

x1s=0.791;
x2s=332.399;
x3s=0.659;
u1s=302.79;
u2s=0.1;

%Eq1
%Derivadas parciais em rela��o a x1
x11=(-Fo/(pi*r^2*x3s))-Ko*(exp(-ER/x2s)); %d/dx1
%Derivadas parciais em rela��o a x2
x12=-(Ko*x1s*ER*exp(-ER/x2s))/(x2s^2); %d/dx2
%Derivadas parciais em rela��o a x3
x13=-(Fo*Co-Fo*x1s)/(pi*r^2*x3s^2); %d/dx3

%Eq 2
%Derivadas parciais em rela��o a x1
x21=(-dH*Ko*exp(-ER/x2s))/(ro*Cp);%d/dx1
%Derivadas parciais em rela��o a x1
x22=(-Fo/(pi*r^2*x3s))-2*U/(r*ro*Cp)-(dH*Ko*x1s*ER*exp(-ER/x2s))/(ro*Cp*x2s^2); %d/dx2
%Derivadas parciais em rela��o a x1
x23=-(Fo*To-Fo*x2s)/(pi*r^2*x3s^2); %d/dx3

%Eq3
x31=0;
x32=0;
x33=0;

A1=[x11 x12 x13;
   x21 x22 x23;
   x31 x32 x33];

%Derivadas parciais em rela��o a u1
B221=2*U/(r*ro*Cp);
%Derivadas parciais em rela��o a u2
B232=-1/(pi*r^2);
%O resto � zero
B2=[0 0;
   B221 0;
   0 B232];

%Derivadas parciais em rela��o a w1 = Fo das Eq 1, 2 e 3
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

%Discretiza��o 

A=A1;
B=B2;
C=C2;
D=D12;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 

B21d=sd.B;
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
C1d=C2d;

disp('Matrizes Discretizadas Ponto de Opera��o 1:')

A1d
B11d
B21d
C2d
D11d
D12d
C1d=C2d;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opera��o 2 

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

%Discretiza��o 

A=A2;
B=B2;
C=C2;
D=D12;

sys=ss(A,B,C,D);
T=1; 
sd = c2d (sys,T,'zoh'); 

B22d=sd.B;
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
C1d=C2d;

disp('Matrizes Discretizadas Ponto de Opera��o 2:')

A2d
B12d
B22d
C2d
D11d
D12d






%V�rtices do Politopo

Ad=[];
Ad(:,:,1)=A1d;
Ad(:,:,2)=A2d;

B2d=[];
B2d(:,:,1)=B21d;
B2d(:,:,2)=B22d;

[n,nn,na]=size(Ad);
[nn,m,nb]=size(B2d);
x0=[0.96 300 0.659]';


%Par�metros

% aloc=0;
x=[]';
y=[];
u=[];
Q=C2'*C2;
R=100000*eye(2);

disp('Sistema Nominal:')

Anominal
B2nominal
C2nominal
N=40; %40
% u_max=2;

%  rng('default')

%Condi��es Nominais
x1s=0.878;
x2s=324.5;
x3s=0.659;
u1s=300;
u2s=0.1;
j=1;

% rng('default')
xk=x0-[x1s x2s x3s]';
return

for k=0:N
 
  
    setlmis([]);
    
    %Decalra��o de Vari�veis
    
    gama=lmivar(1,[1 1]);
    
    [W1,n_W1]=lmivar(1,[n 1]);
    [W2,n_W2]=lmivar(2,[m n]);
    
    %1� LMI 
    
    %j=1;
    lmiterm([-j 1 1 0],1);
    lmiterm([-j 2 1 0],xk);
    lmiterm([-j 2 2 W1],1,1);
    j=j+1;
    
    %2� LMI
    
    for i=1:na
        for l=1:nb
        lmiterm([-j 1 1 W1],1,1);
        lmiterm([-j 2 1 W1],Ad(:,:,i),1);
        lmiterm([-j 2 1 W2],B2d(:,:,l),1);
        lmiterm([-j 2 2 W1],1,1);
        lmiterm([-j 3 1 W1],Q^0.5,1);
        lmiterm([-j 3 2 0],0);
        lmiterm([-j 3 3 gama],1,1);
        lmiterm([-j 4 1 W2],R^0.5,1);
        lmiterm([-j 4 2 0],0);
        lmiterm([-j 4 3 0],0);
        lmiterm([-j 4 4 gama],1,1);
        j=j+1;
        end
    end
    
    %N�o tem restri��es de aloca��o
    
    %2� LMI GAMA > 0
    j=j+1; 
    lmiterm([-j 1 1 gama],1,1);

    %3� LMI W1>0
    j=j+1;
    lmiterm([-j 1 1 W1],1,1);
    
    preditivo=getlmis;
    
    
    %Norma Hinf
    
    %c= mat2dec(LQR_lmi,0,0,0,1);
    c=[1 zeros(1,n_W2-1)]';
    options=[1e-8,3000,1e9,1e-4,1];
    
    [c_opt,x_opt]=mincx(preditivo,c,options);
    %[c_opt,x_opt]=feasp(preditivo,options);
    gama_opt=dec2mat(preditivo,x_opt,gama);
    W1_opt=dec2mat(preditivo,x_opt,W1);
    W2_opt=dec2mat(preditivo,x_opt,W2);
%
% %Crit�rio
% % 
     disp('Norma via LMI:')
     gama_opt
%      W2_=W2_opt
%      W1_=W1_opt
% 
    disp('Ganho Controlador Hinf via LMI:')
    K=W2_opt*inv(W1_opt)



%Aplicar o Controle ao Modelo Nominal 

if k==0
    
    u=K*(xk)+[u1s u2s]';
    u0=u;
 
else
    u=K*(xk)+[u1s u2s]';
    
end

  fcstr=@(t,z)[(Fo*Co/(pi*r^2*z(3)))-(Fo*z(1)/(pi*r^2*z(3)))-(Ko*z(1)*exp(-ER/z(2)));
    (Fo*To/(pi*r^2*z(3)))-(Fo*z(2)/(pi*r^2*z(3)))-((dH/(ro*Cp))*z(1)*exp(-ER/z(2)))+((2*U/(ro*Cp))*(u(1,:)-z(2)));
    Fo/(pi*r^2)-Fo*u(2,:)/(pi*r^2)]

[T,Z]=ode45 (fcstr,[k*Ts (k+1)*Ts],[x0]') 

x(:,k+1)=Z(41,:);%Absoluto 
x0=x(:,k+1);
xk=x(:,k+1)-[x1s x2s x3s]'; % Nas LMIs � relativo


% 
%  vet_x=[x0 x];
%  vet_u=[u0 u];
 


end


% plot(T,vet_x(1,:))
% figure;
% plot(T,vet_x(2,:))
% figure;
% plot(T,vet_x(3,:))
%  
% %  vet_y=[y0 y];


% return
%%

% figure;
% plot(T,Z)
% figure
% plot(T,x)
% figure;
% plot(T,vet_u)
% figure;
% plot(T,vet_y)


return

