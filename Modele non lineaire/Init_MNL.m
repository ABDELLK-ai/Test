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