% Nettoyer l'environnement de travail
clc; % Effacer la fenêtre de commande
clear; % Supprimer toutes les variables de l'espace de travail
close all; % Fermer toutes les figures

% Paramètres initiaux
t = 0:0.01:10; % vecteur de temps
Lf = 0.71; Lr = 2.13; L = Lr + Lf; % les paramètres du véhicule. Lf est l’empattement avant du véhicule, Lr est l’empattement arrière du véhicule, et L est la longueur totale du véhicule (la somme de Lf et Lr).
f = 1/8; lambda = 16; Vx = 90; % f est la fréquence de l’angle du volant, lambda est un facteur de conversion pour l’angle de braquage des roues avant, et Vx est la vitesse longitudinale constante du véhicule.

% Calcul de l'angle du volant
thetav = 2.2*sin(2*pi*f*(t-0.25)); % angle du volant en fonction du temps
thetav(t<0.25 | t>8.25) = 0; % angle du volant est nul pour t<0.25 et t>8.25

% Calcul de la vitesse de lacet et de la position Yg
Bv = thetav/lambda; % angle de braquage des roues avant[^1^][1]
[t, Vlacet, VYg] = Model_MC(Vx, L, Bv, t); % appel de la fonction du modèle cinématique MC
Yg = cumtrapz(t, VYg); % calcul de la position Yg par intégration de la vitesse latérale

% Affichage des résultats
figure(1);
plot(t, Vlacet,'r'); % tracé de la vitesse de lacet en fonction du temps
xlabel('Temps (s)'); ylabel('Vitesse de lacet (m/s)'); % étiquettes des axes
title('Vitesse de lacet en modéle cinématique')
grid on; % grille

figure(2); 
plot(t, Yg,'g'); % tracé de la position Yg en fonction du temps
xlabel('Temps (s)'); ylabel('Position Yg (m)'); % étiquettes des axes
title('Position latéral en modéle cinématique')
grid on; % grille

figure(3); 
plot(t, thetav,'m'); % tracé de l'angle du volant en fonction du temps
xlabel('Temps (s)'); ylabel('Angle volant (°)'); % étiquettes des axes
grid on; % grille
title('Angle du volant en modéle cinématique')

% Fonction du modèle cinématique MC
function [t, Psi_dot_deg, VYg] = Model_MC(Vx, L, Bv, t)
    Bv = deg2rad(Bv); % conversion de l'angle de braquage en radians
    Psi_dot = ((Vx/3.6)/L)*(tan(Bv)); % calcul de la vitesse de lacet
    Psi_dot_deg = rad2deg(Psi_dot); % conversion de la vitesse de lacet en degrés
    Psi = cumtrapz(t, Psi_dot); % calcul de l'angle de lacet par intégration de la vitesse de lacet
    VYg = (Vx/3.6)*sin(Psi); % calcul de la vitesse latérale
end
