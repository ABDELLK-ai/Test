%% Comparaison de vitesse du Modèle linéaire en régime circulaire uniforme
clc,clear all,close all,
Fig=1;
Mt = 1759;
Mf = 1319;
Mr = 440;
Iz = 2638.5;
Lf = 0.71;
Lr = 2.13;
cyf = 94446;
cyr = 48699;
lambda = 16;
L = Lf+Lr;
Dx=200;
Amplitude=[3.5/2.42 3.5/2.30 3.5/2.10 3.5/1.84 3.5/1.58 3.5/1.34 3.5/1.14];
Vx0=[10 30 50 70 90 110 130]/3.6;
T_t=Dx./Vx0;
% Coefficient modèle circulaire uniforme
e1 = Mt*(Lr*cyr-Lf*cyf)/(2*L*cyr*cyf) + L./(Vx0.^2);

% Modele d'etat X_point = A*X + B*U
A = [0 0;1 0];
B1 = [1./e1(1) ; 0]; % Vx0= 10 km/h
B2 = [1./e1(2) ; 0]; % Vx0= 30 km/h
B3 = [1./e1(3) ; 0]; % Vx0= 50 km/h
B4 = [1./e1(4) ; 0]; % Vx0= 70 km/h
B5 = [1./e1(5) ; 0]; % Vx0= 90 km/h
B6 = [1./e1(6) ; 0]; % Vx0= 110 km/h
B7 = [1./e1(7) ; 0]; % Vx0= 130 km/h

C = [0 1];
D = 0;
n = size(A);
% f(X) = AX + BU = X_point
N = 1000;
T = 100;
h = T/N;
t = 0:h:T;
% Condition Initilale
X1=[0;0];
X2=[0;0];
X3=[0;0];
X4=[0;0];
X5=[0;0];
X6=[0;0];
X7=[0;0];

K0=(2.*cyf.*cyr.*Vx0.*L)./(lambda.*2.*cyf.*cyr.*L^2-lambda.*Mt.*Vx0.^2.*(Lf.*cyf-Lr.*cyr));


%%
% Euler trapèze
for k=1:N
    X1(:,k+1) = (eye(n) + h*A)*X1(:,k) + h*B1.*(1/lambda)*beta(k*(T/N),Vx0(1),Dx,T_t(1),Amplitude(1));
    Y1(k+1) = C*X1(:,k+1);
end

for k=1:N
    X2(:,k+1) = (eye(n) + h*A)*X2(:,k) + h*B2.*(1/lambda)*beta(k*(T/N),Vx0(2),Dx,T_t(2),Amplitude(2));
    Y2(k+1) = C*X2(:,k+1);
end

for k=1:N
    X3(:,k+1) = (eye(n) + h*A)*X3(:,k) + h*B3.*(1/lambda)*beta(k*(T/N),Vx0(3),Dx,T_t(3),Amplitude(3));
    Y3(k+1) = C*X3(:,k+1);
end

for k=1:N
    X4(:,k+1) = (eye(n) + h*A)*X4(:,k) + h*B4.*(1/lambda)*beta(k*(T/N),Vx0(4),Dx,T_t(4),Amplitude(4));
    Y4(k+1) = C*X4(:,k+1);
end

for k=1:N
    X5(:,k+1) = (eye(n) + h*A)*X5(:,k) + h*B5.*(1/lambda)*beta(k*(T/N),Vx0(5),Dx,T_t(5),Amplitude(5));
    Y5(k+1) = C*X5(:,k+1);
end

for k=1:N
    X6(:,k+1) = (eye(n) + h*A)*X6(:,k) + h*B6.*(1/lambda)*beta(k*(T/N),Vx0(6),Dx,T_t(6),Amplitude(6));
    Y6(k+1) = C*X6(:,k+1);
end

for k=1:N
    X7(:,k+1) = (eye(n) + h*A)*X7(:,k) + h*B7.*(1/lambda)*beta(k*(T/N),Vx0(7),Dx,T_t(7),Amplitude(7));
    Y7(k+1) = C*X7(:,k+1);
end

for i=1:length(t)
    a1(i)= beta(t(i),Vx0(1),Dx,T_t(1),Amplitude(1))/(lambda*e1(1));
    a2(i)= beta(t(i),Vx0(2),Dx,T_t(2),Amplitude(2))/(lambda*e1(2));
    a3(i)= beta(t(i),Vx0(3),Dx,T_t(3),Amplitude(3))/(lambda*e1(3));
    a4(i)= beta(t(i),Vx0(4),Dx,T_t(4),Amplitude(4))/(lambda*e1(4));
    a5(i)= beta(t(i),Vx0(5),Dx,T_t(5),Amplitude(5))/(lambda*e1(5));
    a6(i)= beta(t(i),Vx0(6),Dx,T_t(6),Amplitude(6))/(lambda*e1(6));
    a7(i)= beta(t(i),Vx0(7),Dx,T_t(7),Amplitude(7))/(lambda*e1(7));
end

Xi1=a1./Vx0(1);
Xi2=a2./Vx0(2);
Xi3=a3./Vx0(3);
Xi4=a4./Vx0(4);
Xi5=a5./Vx0(5);
Xi6=a6./Vx0(6);
Xi7=a7./Vx0(7);

%% Affichage pour plusieurs vitesses
figure(Fig)
Fig=Fig+1;
plot(t.*Vx0(1),Y1(:))
hold on
plot(t.*Vx0(2),Y2(:))
hold on
plot(t.*Vx0(3),Y3(:))
hold on
plot(t.*Vx0(4),Y4(:))
hold on
plot(t.*Vx0(5),Y5(:))
hold on
plot(t.*Vx0(6),Y6(:))
hold on
plot(t.*Vx0(7),Y7(:))
hold off
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title("Modèle linéaire en régime circulaire uniforme: Position latérale Yg")
xlabel('X_G(m)')
ylabel('Y_G(m)')
axis([0 210 0 4])
grid on



figure(Fig)
Fig=Fig+1;
plot(t*Vx0(1),a1)
hold on
plot(t*Vx0(2),a2)
hold on
plot(t*Vx0(3),a3)
hold on
plot(t*Vx0(4),a4)
hold on
plot(t*Vx0(5),a5)
hold on
plot(t*Vx0(6),a6)
hold on
plot(t*Vx0(7),a7)
hold off
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title("Modèle linéaire en régime circulaire uniforme : Accélération latérale")
xlabel('X_G(m)')
ylabel('\Gamma(t)')
axis([0 210 -0.8 0.8])
grid on



figure(Fig)
Fig=Fig+1;
plot(t*Vx0(1),(180/pi)*Xi1)
hold on
plot(t*Vx0(2),(180/pi)*Xi2)
hold on
plot(t*Vx0(3),(180/pi)*Xi3)
hold on
plot(t*Vx0(4),(180/pi)*Xi4)
hold on
plot(t*Vx0(5),(180/pi)*Xi5)
hold on
plot(t*Vx0(6),(180/pi)*Xi6)
hold on
plot(t*Vx0(7),(180/pi)*Xi7)
hold off
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title("Modèle linéaire en régime circulaire uniforme : Vitesse de lacet")
xlabel('X_G(m)')
ylabel('\Psi(°/t)')
axis([0 210 -1.2 1.2])
grid on
