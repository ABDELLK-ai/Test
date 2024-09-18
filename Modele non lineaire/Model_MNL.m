%%
clear
close all
clc

%% Chargement des paramètres de simulation
Init_MNL


%%
%date finale de simulation
tspan = [0 10];

y0 = [0 0 0 0];      %condition initiale

% Vitesse de lacet y1 et postion latérale y2
[t,y] = ode45(@(t,y) Model_MNL(t, y), tspan, y0);


%% accélération linéaire latérale gammat
%commande u
u = zeros(length(t),1);
for i=1:length(t)
    if ((t(i) >= 0.2) && (t(i) <= 8.2))
        u(i) = Amp*sin(2*pi*(t(i)-0.2)/T);
        u(i) = (1/lambda)*u(i);
    end
end
%accélération gamma_t
gamma_t = 2*cyf/Mt*deg2rad(u) -2*(cyf+cyr)/(Mt*Vx0)*y(:,3) + (2*(-cyf*Lf+cyr*Lr)/(Mt*Vx0))*y(:,2);

%Plot the results
figure (Name='Résultats de simulation du modèle non linéaire')

plot(t, y(:,4))
xlabel('t(s)')
ylabel('Y_g (m)')
title('position latérale modèle non linéaire')

plot(t, gamma_t)
xlabel('t(s)')
ylabel('Gamma_t (m/s²)')
title('accélération latérale modèle non linéaire')

%% Définition du modèle
function dydt = Model_MNL(t, y)
    Init_MNL

    if (t >= 0.2) && (t <= 8.2)
        u = Amp * sin(2*pi*(t-0.2)/T);
        u = (1/lambda) * u;
    else
        u = 0;
    end

    dydt = [0; 0; 0; 0];

    % y(1) = psi        % angle de lacet
    % y(2) = psidot     % vitesse d'angle de lacet
    % y(3) = Vy         % vitesse latérale
    % y(4) = yG         % position d'angle de lacet

    dydt(1) = y(2);

    % dydt(2) = C_sys/Iz;
    dydt(2) = (2*Lf*cyf/Iz) * deg2rad(u) + 2*(-Lf*cyf+Lr*cyr)/(Vx0*Iz) * y(3) - 2*(Lf*Lf*cyf+Lr*Lr*cyr)/(Vx0*Iz) * y(2);

    % dydt(3) = 1/Mt*F_sys-Vx0*y(2);
    dydt(3) = 2*cyf/Mt * deg2rad(u) - 2*(cyf+cyr)/(Mt*Vx0) * y(3) + (2*(-cyf*Lf+cyr*Lr)/(Mt*Vx0)-Vx0) * y(2);

    dydt(4) = Vx0 * sin(y(1)) + y(3) * cos(y(1));

    dp11 = u - atan((y(3) + Lf*y(2)) / (Vx - Lf*y(2)));
    dp12 = u - atan((y(3) + Lf*y(2)) / (Vx + Lf*y(2)));
    dp21 = -atan((y(3) - Lr*y(2)) / (Vx - Lf*y(2)));
    dp22 = -atan((y(3) - Lr*y(2)) / (Vx - Lf*y(2)));

    By11 = 0;
    Dy11 = 0;
    Cy11 = 0;
    Dy11 = 0;
    Ey11 = 0;

    Fz = Mt * g;
    D = a1*Fz*Fz + a2*Fz;
    C = a0;
    BCD = a3*sin(2*atan(Fz/a4))*(1 - a5*abs(y(1)));
    By = BCD / C / D;
    E = min((a6*Fz + a7), 1);
    Sh = a8*y(1) + a9*Fz + a10;
    Sv = a12*Fz + a13 + (a112*Fz^2 + a111*Fz) * y(1);
end

%% Exécutable de définition des valeurs des paramètres du modèle non linéaire MNL

% Paramètres nominaux du véhicule
Mt     = 1759;     % kg, masse totale du véhicule
Mf     = 1319;     % kg, masse de la partie avant du véhicule
Mr     = 440;      % kg, masse de la partie arrière du véhicule
Iz     = 2638.5;   % kg.m², moment d'inertie de l'axe vertical
Lf     = 0.71;     % m, longueur de l'empattement avant
Lr     = 2.13;     % m, longueur de l'empattement arrière
L      = Lf + Lr;
cyf    = 94446;    % N/rad, coefficient de rigidité de dérive des pneumatiques avant
cyr    = 48699;    % N/rad, coefficient de rigidité de dérive des pneumatiques arrière
lambda = 16;       % facteur de démultiplication entre les angles volant et roues
Vx     = 90/3.6;   % m/s
Vx0    = Vx;
Dx     = 200;      % m
Amp    = 2.25;     % Amplitude du signal sinus d'entrée
T      = Dx/Vx;    % période du signal sinus (angle volant), calculé par T = Dx/Vx avec Dx = 200 m et Vx = 90 km/h

g = 9.81;          % m/s²

% Coefficients des efforts pneumatiques du modèle de Pacejka
a0   = 1.998;
a1   = -33.85;
a2   = 1198;
a3   = 2258;
a4   = 10.74;
a5   = 0.01399;
a6   = -0.1693;
a7   = 1;
a8   = -0.03009;
a9   = -0.009786;
a10  = -0.1149;
a111 = -10.85;
a112 = 3.225;
a12  = 3.225;
a13  = 34.78;

