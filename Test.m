% Définition des paramètres
L = 2.84; % longueur du véhicule en m
Vx = 25; % vitesse longitudinale constante en m/s
T = 10; % durée de la simulation en s
dt = 0.01; % pas de temps en s
t = 0:dt:T; % vecteur temps
N = length(t); % nombre de points

% Définition de l'angle de braquage
beta_v = 10 * sin(2 * pi * t/T); % angle de braquage en radian

% Initialisation des variables
psi = zeros(1,N); % angle de lacet en radian
y_G = zeros(1,N); % position latérale en m

% Boucle de simulation
for i = 2:N
    % Calcul de la vitesse de lacet
    psi_dot = Vx/L * tan(beta_v(i));
    
    % Calcul de l'angle de lacet
    psi(i) = psi(i-1) + psi_dot * dt;
    
    % Calcul de la dérivée de la position latérale
    y_G_dot = Vx * sin(psi(i));
    
    % Calcul de la position latérale
    y_G(i) = y_G(i-1) + y_G_dot * dt;
end

% Affichage des résultats
figure
subplot(2,1,1)
plot(t,psi*180/pi)
xlabel('Temps (s)')
ylabel('Angle de lacet (deg)')
title('Vitesse du lacet')
grid on

subplot(2,1,2)
plot(t,y_G)
xlabel('Temps (s)')
ylabel('Position latérale (m)')
title('Trajectoire du modèle cinématique')
grid on
