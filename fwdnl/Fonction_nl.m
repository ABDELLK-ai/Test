% Définir les paramètres du véhicule
Mt = 1759;   % Masse totale du véhicule
Mf = 1319;   % Masse avant du véhicule
Mr = 440;    % Masse arrière du véhicule
Iz = 2638.5; % Moment d'inertie
Lf = 0.71;   % Distance entre le centre de gravité et l'essieu avant
Lr = 2.13;   % Distance entre le centre de gravité et l'essieu arrière
cyf = 94446; % Coefficient de raideur latérale de l'essieu avant
cyr = 48699; % Coefficient de raideur latérale de l'essieu arrière
L = 16;      % Constante lambda
mu = 100;    % Coefficient de frottement
f0y = 0;     % Force latérale externe
l = 1; % Largeur du véhicule
Vx0 = 90 / 3.6;

% Coefficients
a0 = 1.998;
a1 = -33.85;
a2 = 1198;
a3 = 2258;
a4 = 10.74;
a5 = 0.01399;
a6 = -0.1693;
a7 = 1;
a8 = -0.03009;
a9 = -0.009786;
a10 = -0.1149;
a111 = -10.85;
a112 = -0.1834;
a12 = 3.225;
a13 = 34.78;


Fz = 0.87; 
gamma = 0.1; 
Delta = 15;
D = a1 * (Fz^2) + a2 * Fz;
C = a0;
BCD = (a3 * sin(2 * atan((Fz / a4))))*(1 - a5 * abs(gamma));
B = BCD/(C*D);
E = min(a6 * Fz + a7, 1);
Sh = a8 * gamma + a9 * Fz + a10;
Sv = a12 * Fz + a13 + (a112 * (Fz^2) + a111 * Fz)*gamma;




