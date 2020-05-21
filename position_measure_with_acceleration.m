% x_k = Ax_(x-1) + Bu_(k-1) + w_(k_1)
% Y_k = Cx_k + v_k
%Pomiar pozycji naprzyk³ad z odbiornika GPS
%tutaj wystepuje sterowanie -> obiekt siê porusza
%S - pozycja
%V - predkosc 
% Sk=S_(k-1)+[V_(k-1)+w_(k-1)]dt +ak*dt^2/2
% Vk=Vk+vk + ak*dt
% ak = ak
dt=1;
N=100;

A=[1 dt dt^2/2; 0 1 dt; 0 0 1];
%x=[S;V;a];
C=[1 1 1];

R=2000;
q = 0.01;               %wariacja procesu
%Q=diag([1,1,1])*q;        %macierz kowariancji
%Poszczególne wariacje wyliczone:
Q = [dt^5/20 dt^4/8 dt^3/6; dt^4/8 dt^3/3 dt^2*2; dt^3/6 dt^2/2 dt]*q;

% zk  - pomiar pozycji

%% DANE
u = zeros(1,N-1);
for i=1:N-1
    if i<20
        u(i+1)=u(i)+rand(1)*15-rand(1)*6;
    elseif i<40
        u(i+1)=u(i)-rand(1)*15+rand(1)*7;
    elseif i<80
        u(i+1)=5*rand(1)*10;
    else
        u(i+1)=u(i)+rand(1)*10;
    end
end

%WEKTORY WYNIKOW
S = zeros(1,N);
V = zeros(1,N);

%% KALMAN 
%inicjacja:
x0=[0; 0;0 ];
P0=[1 0 0; 0 1 0; 0 0 1]*q;
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
    S(i) = xpost(1);
    V(i) = xpost(2);
    a(i) = xpost(3);
end

%% Wykres
x = 0:1:N-1;
plot(x,u);
hold on
grid on
plot(x,S);
plot(x,V);
plot(x,a);
title("Pomiar przebytej odleg³oœci ")
xlabel("Próbki [s]")
ylabel("Pozycja [m]")
legend(["Pomiar pozycji", "Kalman-po³o¿enie", "Predkoœæ", "Przyœpieszenie"])
