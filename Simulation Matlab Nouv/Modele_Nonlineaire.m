% Nettoyer l'environnement de travail
clc; % Effacer la fenêtre de commande
clear; % Supprimer toutes les variables de l'espace de travail
close all; % Fermer toutes les figures


% Paramètres du véhicule
Mt = 1759; % Masse totale du véhicule en kg
Mf = 1319; % Masse avant du véhicule en kg
Mr = 440; % Masse arrière du véhicule en kg
Iz = 2638.5; % Moment d'inertie du véhicule en kg.m^2
Lf = 0.71; % Distance entre le centre de gravité et l'essieu avant en m
Lr = 2.13; % Distance entre le centre de gravité et l'essieu arrière en m
L = Lf+Lr
cyf = 94446; % Rigidité latérale de l'essieu avant en N/rad
cyr = 48699; % Rigidité latérale de l'essieu arrière en N/rad
lambda = 16; % Rapport de réduction de la direction

% Vitesse longitudinale
Vxt = 90/3.6; % Conversion de la vitesse longitudinale de km/h en m/s
Vx0 = Vxt; % Vitesse longitudinale initiale en m/s

% Vecteur de temps
t = 0:0.001:10; % Vecteur de temps de 0 à 10 secondes avec un pas de 0.001 seconde

% Fréquence
f = 1/8;

% Initialisation de l'angle de volant
thetav = zeros(1,length(t));

% Calcul de l'angle de volant en fonction du temps
for i = 1:length(t)
    if (0.25<t(i)) && (t(i)<8.25)
        thetav(i) = 1.94*sin(2*pi*f*(t(i)-0.25));
    else
        thetav(i) = 0;
    end
end

% Calcul de l'angle de braquage des roues avant en fonction de l'angle de volant
Bvt = thetav/lambda;

% Calcul de l'accélération latérale en utilisant la deuxième hypothèse
Theta_t= Bvt/((Mt*(Lr*cyr-Lf*cyf)/(2*L*cyr*cyf)+L/(Vx0^2)));

% Tracé de la vitesse de lacet
figure(1);
plot(t,Theta_t/(Vx0),'r');
xlabel("Temps (s)");
ylabel("Vitesse de lacet (rad/s)");
title("Vitesse de lacet");
grid on
% Tracé de l'accélération latérale
figure(2);
plot(t,Theta_t/50,'g');
xlabel("Temps (s)");
ylabel("Accélération latérale (m/s²)");
title("Accélération latérale");
grid on
% Le pas d'intégration
h = t(2) - t(1); 

% Initialisation du vecteur vitesse
vitesse = zeros(size(t)); 

% Calcul de la vitesse en utilisant la méthode des trapèzes
for i = 2:length(t)
    vitesse(i) = vitesse(i-1) + (h/2)*(Theta_t(i-1) + Theta_t(i)); 
end

% Tracé de la vitesse latérale
figure(3);
plot(t,vitesse,'b');
xlabel("Temps (s)");
ylabel("Vitesse latérale (m/s)");
title("Vitesse latérale du véhicule");
grid on
% Initialisation du vecteur position
position = zeros(size(t));

% Calcul de la position en utilisant la méthode des trapèzes
for i = 2:length(t)
    position(i) = position(i-1) + (h/2)*(vitesse(i-1) + vitesse(i)); 
end

% Tracé de la position latérale
figure(4);
plot(t,position/50,'m');
xlabel('Temps (s)');
ylabel('Position latérale (m)');
title("Position latérale du véhicule");
grid on
