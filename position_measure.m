% x_k = Ax_(x-1) + Bu_(k-1) + w_(k_1)
% Z_k = Hx_k + v_k
%Pomiar pozycji naprzyk³ad z odbiornika GPS
%tutaj wystepuje sterowanie -> obiekt siê porusza
%S - pozycja
%V - predkosc 
% Sk=S_(k-1)+[V_(k-1)+w_(k-1)dt
% Vk=Vk+vk
dt=1;
N=100;

A=[1 dt; 0 1];
%x=[S;V];
C=[1 1];

R=0.01;
q = 0.1;               %wariacja procesu
Q=diag([1,1])*q;        %macierz kowariancji
% u  - pomiar pozycji

%% DANE
u = zeros(1,N-1);
for i=1:N-1
    if i<20
        u(i+1)=u(i)+rand(1)*15-rand(1)*6;
    elseif i<40
        u(i+1)=u(i)-rand(1)*15+rand(1)*7;
    elseif i<80
        u(i+1)=5*rand(1)*20;
    else
        u(i+1)=u(i)+rand(1)*10;
    end
end

%WEKTORY WYNIKOW
S = zeros(1,N);
V = zeros(1,N);

%% KALMAN 
%inicjacja:
x0=[0; 0];
P0=[1 0; 0 1]*q;
xpost=x0;
Ppost=P0;

for i=1:N
    %predykcja
    xpri = A*xpost;
    Ppri = A*Ppost*A'+Q;
    %Korekcja
    Kk = Ppri*C'*((C*Ppri*C'+R)^(-1));
    xpost = xpri+Kk*(u(i)-C*xpri);
    Ppost = (1-Kk*C)*Ppri;
    S(i)=xpost(1);
    V(i)= xpost(2);
end

%% Wykres
x = 0:1:N-1;
plot(x,u);
hold on
grid on
plot(x,S);
plot(x,V);
title("Pomiar przebytej odleg³oœci ")
xlabel("Próbki [s]")
ylabel("Pozycja [m]")
legend(["Pomiar pozycji", "Kalman-po³o¿enie","Prêdkoœæ"])
