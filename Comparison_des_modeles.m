%% Comparaison des trois modèles
clear all, close all, clc,
Fig=1;
% Paramatre nominaux de la dynamique latarale du véhicule
Mt=1759; Iz= 2638.5; Lf=0.71; Lr= 2.13; c_yf=94446; c_yr=48699; lambda= 16; Dx=200;L=Lr+Lf;
Vx0=90/3.6;
T_t=Dx./Vx0;
Amplitude=2.2;

%% Modèle linéaire
% partie temporelle
y_ode=[0;0;0;0];
[tlin,Ylin]=ode23(@(t,y) ModeLin(t,y,Mt, Iz, Lf, Lr, c_yf, c_yr, lambda, Vx0, Dx,T_t,Amplitude), [0 100], y_ode);
u=beta(tlin,Vx0,Dx,T_t,Amplitude);
for i=1:length(tlin)
    alin(i)= (2*c_yf*u(i)/(lambda*Mt)+ (2*(-c_yf*Lf +c_yr*Lr)/(Vx0*Mt))*Ylin(i,2)- 2*((c_yf+c_yr)/(Vx0*Mt))*Ylin(i,3)) ;
end
% Partie fréquentielle
a22= (-2*(Lf^2*c_yf+Lr^2*c_yr)/Iz)./Vx0;
a23= (2*(-Lf*c_yf+Lr*c_yr)/Iz)./Vx0;
a32= 2*(-c_yf*Lf+c_yr*Lr)./(Mt.*Vx0)-Vx0;
a33= -2*(c_yf+c_yr)./(Mt.*Vx0);
b2= 2*Lf*c_yf/(Iz*lambda);
b3= 2*c_yf/(Mt*lambda);
Alin= [0 1 0 0;0 a22 a23 0;0 a32 a33 0;Vx0 0 1 0];% vitesse Vx0=90 km/h
Blin= [0;b2;b3;0];
Clin= [0 0 0 1];
C2lin=[0 1 0 0];
D=0;
Slin=ss(Alin,Blin,Clin,D);
S2lin=ss(Alin,Blin,C2lin,D);
%% Modèle circulaire uniforme
e = Mt*(Lr*c_yr-Lf*c_yf)/(2*L*c_yr*c_yf) + L./(Vx0.^2);

% Modele d'etat X_point = A*X + B*U
Acir = [0 0;1 0];
Bcir = [1./e ; 0]; % Vx0= 90 km/h
Ccir = [0 1];D = 0;
n = size(Acir);
% f(X) = AX + BU = X_point
N = 1000;T = 100;h = T/N; tcir = 0:h:T; Xcir=[0;0];
for k=1:N
    Xcir(:,k+1) = (eye(n) + h*Acir)*Xcir(:,k) + h*Bcir.*(1/lambda)*beta(k*(T/N),Vx0,Dx,T_t,Amplitude);
    Ycir(k+1) = Ccir*Xcir(:,k+1);
end
for i=1:length(tcir)
    acir(i)= beta(tcir(i),Vx0,Dx,T_t,Amplitude)/(lambda*e);
end
dXicir=acir./Vx0;

% partie fréquentielle
K0=(2*c_yf*c_yr*Vx0*L)/(lambda*2*c_yf*c_yr*L^2-lambda*Mt*Vx0^2*(Lf*c_yf-Lr*c_yr));
Scir=ss(Acir,Bcir,Ccir,D);
S2cir=tf([K0],[1]);
%% Modèle cinématique
yode=[0;0];
opts = odeset('RelTol',1e-10);
[tcin,Ycin]=ode23(@(t,y) ModeCin(t,y,Vx0, Dx,T_t,L,Amplitude),[0 10],yode);

for i=1:length(tcin)
    dXicin(i)= (Vx0/L)*tan(beta(tcin(i),Vx0,Dx,T_t,Amplitude)/lambda);
end

Scin=tf([Vx0^2],[lambda*L 0 0]);
S2cin=tf([Vx0],[lambda*L]);
% affichage des positions latérales

figure(Fig)
Fig=Fig+1;
plot(tlin,Ylin(:,4),'cyan')
hold on
plot(tcir,Ycir,'b--')
hold on
plot(tcin,Ycin(:,2),'green--')
hold off
legend('Modèle linéaire','Modèle circulaire uniforme','Modèle cinématique')
title('Position latérale à une vitesse initiale de 90km/h')
xlabel('t(s)')
ylabel('Y_G(m)')
axis([0 10 0 6])
grid on

%% Affichage de l'accélération latérale
figure(Fig)
Fig=Fig+1;
plot(tlin,alin)
hold on
plot(tcir,acir)
hold off
legend('Modèle linéaire','Modèle circulaire uniforme')
title('Accélération latérale à une vitesse initiale de 90km/h')
xlabel('t(s)')
ylabel('Acc_y(m/s²)')
axis([0 10 -0.4 0.4])
grid on
%% Affichage de la vitesse de lacet
figure(Fig)
Fig=Fig+1;
plot(tlin,(180/pi)*Ylin(:,2))
hold on
plot(tcir,(180/pi)*dXicir)
hold on
plot(tcin,(180/pi)*dXicin)
hold off
legend('Modèle linéaire','Modèle circulaire uniforme','Modèle cinématique')
title('Vitesse de lacet à une vitesse initiale de 90km/h')
xlabel('t(s)')
ylabel('dXi(°/s)')
axis([0 10 -1.5 1.5])
grid on
% %% Diagramme de bode
figure(Fig)
 Fig=Fig+1;
 bode(Slin,':',Scir,'--',Scin,'r',{10^-4,10^4})
 title('Digramme de Bode entre position latérale et l''angle du volant')
 legend('Modèle linéaire','Modèle circulaire uniforme','Modèle cinématique')
 grid on
 
 figure(Fig)
 Fig=Fig+1;
 bode(S2lin,':',S2cir,'--',S2cin,'r',{10^-4,10^4})
 title('Digramme de Bode entre la vitesse de lacet et l''angle du volant')
 legend('Modèle linéaire','Modèle circulaire uniforme','Modèle cinématique')
 grid on





function Y= beta(t,Vx0,Dx,T_t,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        Y(i)=0;
    else
        Y(i)= (pi/180)*A*sin(2*pi*(t(i)-(5/Vx0))/T_t);
    end
end
end


function dY=ModeCin(t,Y,Vx0, Dx,T_t,L,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= 0;
        dY(2,i)= Vx0*sin(pi*Y(1,i)/180);

    else
        dY(1,i)= Vx0*tan(beta(t(i),Vx0,Dx,T_t,A)/16)/L;
        dY(2,i)= Vx0*sin(Y(1,i));

    end
end

end


function dY = ModeLin(t,Y,Mt, Iz, Lf, Lr, c_yf, c_yr, lambda, Vx0, Dx,T_t,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i)  -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    else
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i) +2*Lf*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lambda*Iz);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i) -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i) + 2*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lambda*Mt);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    end
end
end

function dY = ModeLinCir(t,Y,Mt,Lf,Lr,c_yf,c_yr,lambda,Vx0,Dx,T_t,A)

L=Lr+Lf;
e=Mt*((Lr*c_yr-Lf*c_yf)/(2*L*c_yr*c_yf))+(L/(Vx0*Vx0));

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= 0;
        dY(2,i)= Y(1,i);
    else
        dY(1,i)= beta(t(i),Vx0,Dx,T_t,A)/(lambda*e);
        dY(2,i)= Y(1,i);
    end
end
end

