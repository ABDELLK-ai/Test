%% Modélisation Linéaire

clear all, close all, clc,
Fig=1;
% Paramatre nominaux de la dynamique latarale du véhicule
Mt=1759; Iz= 2638.5; Lf=0.71; Lr= 2.13; c_yf=94446; c_yr=48699; lmbd= 16; Dx=200;L=Lr+Lf;
Vx0=[10 30 50 70 90 110 130]/3.6;
T_t=Dx./Vx0;
Amplitude=[3.5/2.42 3.5/2.30 3.5/2.10 3.5/1.84 3.5/1.58 3.5/1.34 3.5/1.14];
% Modélisation d'état du système
% Ma variable d'état pour le modèle linéaire : y=[Xi,Xi',Vy,Yg] => y'=[Xi',Xi'',Vy',Yg']

% t=0:0.1:10;

a22= (-2*(Lf^2*c_yf+Lr^2*c_yr)/Iz)./Vx0;
a23= (2*(-Lf*c_yf+Lr*c_yr)/Iz)./Vx0;
a32= 2*(-c_yf*Lf+c_yr*Lr)./(Mt.*Vx0)-Vx0;
a33= -2*(c_yf+c_yr)./(Mt.*Vx0);
b2= 2*Lf*c_yf/(Iz*lmbd);
b3= 2*c_yf/(Mt*lmbd);

A1= [0 1 0 0;0 a22(1) a23(1) 0;0 a32(1) a33(1) 0;Vx0(1) 0 1 0];% vitesse Vx0=10 km/h
A2= [0 1 0 0;0 a22(2) a23(2) 0;0 a32(2) a33(2) 0;Vx0(2) 0 1 0];% vitesse Vx0=30 km/h
A3= [0 1 0 0;0 a22(3) a23(3) 0;0 a32(3) a33(3) 0;Vx0(3) 0 1 0];% vitesse Vx0=50 km/h
A4= [0 1 0 0;0 a22(4) a23(4) 0;0 a32(4) a33(4) 0;Vx0(4) 0 1 0];% vitesse Vx0=70 km/h
A5= [0 1 0 0;0 a22(5) a23(5) 0;0 a32(5) a33(5) 0;Vx0(5) 0 1 0];% vitesse Vx0=90 km/h
A6= [0 1 0 0;0 a22(6) a23(6) 0;0 a32(6) a33(6) 0;Vx0(6) 0 1 0];% vitesse Vx0=110 km/h
A7= [0 1 0 0;0 a22(7) a23(7) 0;0 a32(7) a33(7) 0;Vx0(7) 0 1 0];% vitesse Vx0=130 km/h
B= [0;b2;b3;0];
C= [0 1 0 0;0 0 0 1];
D=0;
% u=beta(t,Vx0(5),Dx,T_t(5),1); % Entrée du sytème (angle du volant)
% figure(Fig)
% Fig=Fig+1;
% plot(t,(180/pi)*u)
% title('Angle du volant')
% xlabel('t(s)')
% ylabel('teta(°)')
% grid on
%%
Y1=struct('Y1ode',[],'Y2ode',[],'Y3ode',[],'Y4ode',[],'Y5ode',[],'Y6ode',[],'Y7ode',[]);
t1=struct('t1ode',[],'t2ode',[],'t3ode',[],'t4ode',[],'t5ode',[],'t6ode',[],'t7ode',[]);
a1=struct('a1',[],'a2',[],'a3',[],'a4',[],'a5',[],'a6',[],'a7',[]);
u=struct('u1',[],'u2',[],'u3',[],'u4',[],'u5',[],'u6',[],'u7',[]);


%%
% Résolution par la méthode Runge-Kutta d'ordre 3
y_ode=[0;0;0;0];% Initalisation de mes variable d'état à zéro
% tspan=[0 10];
option = odeset('InitialStep',1e-3);
[t1.t1ode,Y1.Y1ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(1), Dx,T_t(1),Amplitude(1)), [0 100], y_ode,option);
[t1.t2ode,Y1.Y2ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(2), Dx,T_t(2),Amplitude(2)), [0 100], y_ode,option);
[t1.t3ode,Y1.Y3ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(3), Dx,T_t(3),Amplitude(3)), [0 100], y_ode,option);
[t1.t4ode,Y1.Y4ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(4), Dx,T_t(4),Amplitude(4)), [0 100], y_ode,option);
[t1.t5ode,Y1.Y5ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(5), Dx,T_t(5),Amplitude(5)), [0 100], y_ode,option);
[t1.t6ode,Y1.Y6ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(6), Dx,T_t(6),Amplitude(6)), [0 100], y_ode,option);
[t1.t7ode,Y1.Y7ode]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0(7), Dx,T_t(7),Amplitude(7)), [0 100], y_ode,option);

% Entrée de mon système (angle du volant)
u.u1=beta(t1.t1ode,Vx0(1),Dx,T_t(1),3.5/2.42);
u.u2=beta(t1.t2ode,Vx0(2),Dx,T_t(2),3.5/2.30);
u.u3=beta(t1.t3ode,Vx0(3),Dx,T_t(3),3.5/2.10);
u.u4=beta(t1.t4ode,Vx0(4),Dx,T_t(4),3.5/1.84);
u.u5=beta(t1.t5ode,Vx0(5),Dx,T_t(5),3.5/1.58);
u.u6=beta(t1.t6ode,Vx0(6),Dx,T_t(6),3.5/1.34);
u.u7=beta(t1.t7ode,Vx0(7),Dx,T_t(7),3.5/1.14);

%%
% Calcul de l'accélération latérale
for i=1:length(t1.t1ode)
    a1.a1(i)= (2*c_yf*u.u1(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(1)*Mt))*Y1.Y1ode(i,2)- 2*((c_yf+c_yr)/(Vx0(1)*Mt))*Y1.Y1ode(i,3)) ;
end

for i=1:length(t1.t2ode)
    a1.a2(i)= (2*c_yf*u.u2(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(2)*Mt))*Y1.Y2ode(i,2)- 2*((c_yf+c_yr)/(Vx0(2)*Mt))*Y1.Y2ode(i,3)) ;
end

for i=1:length(t1.t3ode)
    a1.a3(i)= (2*c_yf*u.u3(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(3)*Mt))*Y1.Y3ode(i,2)- 2*((c_yf+c_yr)/(Vx0(3)*Mt))*Y1.Y3ode(i,3)) ;
end

for i=1:length(t1.t4ode)
    a1.a4(i)= (2*c_yf*u.u4(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(4)*Mt))*Y1.Y4ode(i,2)- 2*((c_yf+c_yr)/(Vx0(4)*Mt))*Y1.Y4ode(i,3)) ;
end

for i=1:length(t1.t5ode)
    a1.a5(i)= (2*c_yf*u.u5(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(5)*Mt))*Y1.Y5ode(i,2)- 2*((c_yf+c_yr)/(Vx0(5)*Mt))*Y1.Y5ode(i,3)) ;
end

for i=1:length(t1.t6ode)
    a1.a6(i)= (2*c_yf*u.u6(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(6)*Mt))*Y1.Y6ode(i,2)- 2*((c_yf+c_yr)/(Vx0(6)*Mt))*Y1.Y6ode(i,3)) ;
end

for i=1:length(t1.t7ode)
    a1.a7(i)= (2*c_yf*u.u7(i)/(lmbd*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0(7)*Mt))*Y1.Y7ode(i,2)- 2*((c_yf+c_yr)/(Vx0(7)*Mt))*Y1.Y7ode(i,3)) ;
end
%% V1
% Figure de la vitesse du lacet
figure(1)
plot(t1.t1ode*Vx0(1),(180/pi)*Y1.Y1ode(:,2))
hold on % Ajouter cette ligne pour garder la figure ouverte
title('Modèle linéaire : Vitesse de lacet')
xlabel('t(s)')
ylabel('V_lacet(°/s)')
xlim([5 250]) % Définir la limite de l'axe x de 5 à 15
grid on

% Figure de la position latérale
figure(2)
plot(t1.t1ode*Vx0(1),Y1.Y1ode(:,4))
hold on % Ajouter cette ligne pour garder la figure ouverte
title('Modèle linéaire : Position latérale')
xlabel('t(s)')
ylabel('Y_G(m)')
xlim([5 250]) % Définir la limite de l'axe x de 5 à 15
grid on

% Figure de l'accélération latérale
figure(3)
plot(t1.t1ode*Vx0(1),a1.a1)
hold on % Ajouter cette ligne pour garder la figure ouverte
title('Modèle linéaire : Accélération latérale')
xlabel('t(s)')
ylabel('Acc latérale(t)')
xlim([5 250]) % Définir la limite de l'axe x de 5 à 15
grid on

%% V2
% Ajouter les courbes pour V2 aux figures existantes
figure(1)
plot(t1.t2ode*Vx0(2),(180/pi)*Y1.Y2ode(:,2))
figure(2)
plot(t1.t2ode*Vx0(2),Y1.Y2ode(:,4))
figure(3)
plot(t1.t2ode*Vx0(2),a1.a2)

%% V3
% Ajouter les courbes pour V3 aux figures existantes
figure(1)
plot(t1.t3ode*Vx0(3),(180/pi)*Y1.Y3ode(:,2))
figure(2)
plot(t1.t3ode*Vx0(3),Y1.Y3ode(:,4))
figure(3)
plot(t1.t3ode*Vx0(3),a1.a3)

%% V4
% Ajouter les courbes pour V4 aux figures existantes
figure(1)
plot(t1.t4ode*Vx0(4),(180/pi)*Y1.Y4ode(:,2))
figure(2)
plot(t1.t4ode*Vx0(4),Y1.Y4ode(:,4))
figure(3)
plot(t1.t4ode*Vx0(4),a1.a4)

%% V5
% Ajouter les courbes pour V5 aux figures existantes
figure(1)
plot(t1.t5ode*Vx0(5),(180/pi)*Y1.Y5ode(:,2))
figure(2)
plot(t1.t5ode*Vx0(5),Y1.Y5ode(:,4))
figure(3)
plot(t1.t5ode*Vx0(5),a1.a5)

%% V6
% Ajouter les courbes pour V6 aux figures existantes
figure(1)
plot(t1.t6ode*Vx0(6),(180/pi)*Y1.Y6ode(:,2))
figure(2)
plot(t1.t6ode*Vx0(6),Y1.Y6ode(:,4))
figure(3)
plot(t1.t6ode*Vx0(6),a1.a6)

%% V7
% Ajouter les courbes pour V7 aux figures existantes
figure(1)
plot(t1.t7ode*Vx0(7),(180/pi)*Y1.Y7ode(:,2))
figure(2)
plot(t1.t7ode*Vx0(7),Y1.Y7ode(:,4))
figure(3)
plot(t1.t7ode*Vx0(7),a1.a7)

% Ajouter une légende pour chaque sous-figure
figure(1)
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
figure(2)
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
figure(3)
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
%% Représentation fréquentielle
G11=ss(A1,B,C(1,:),D);
G21=ss(A1,B,C(2,:),D);

G12=ss(A2,B,C(1,:),D);
G22=ss(A2,B,C(2,:),D);

G13=ss(A3,B,C(1,:),D);
G23=ss(A3,B,C(2,:),D);

G14=ss(A4,B,C(1,:),D);
G24=ss(A4,B,C(2,:),D);

G15=ss(A5,B,C(1,:),D);
G25=ss(A5,B,C(2,:),D);

G16=ss(A6,B,C(1,:),D);
G26=ss(A6,B,C(2,:),D);

G17=ss(A7,B,C(1,:),D);
G27=ss(A7,B,C(2,:),D);

figure(Fig)
Fig=Fig+1;
bode(G11,G12,G13,G14,G15,G16,G17,{1,1000})
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title('Diagramme de Bode entre la vitesse de lacet et l''angle du volant')
grid on

figure(Fig)
Fig=Fig+1;
bode(G21,G22,G23,G24,G25,G26,G27,{10^-4,10^4})
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title('Digramme de Bode entre position latérale et l''angle du volant')
grid on


function dY = ModeLin(t,Y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0, Dx,T_t,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i)  -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    else
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i) +2*Lf*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lmbd*Iz);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i) -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i) + 2*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lmbd*Mt);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    end
end
end

