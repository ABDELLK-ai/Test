% Définir l'intervalle de temps et les conditions initiales pour le problème de Cauchy
tspan = [0  10];
y0 = [0 0 0 0];

% Résoudre le système d'équations différentielles ordinaires avec la fonction ode45
[t, y] = ode45(@Model_NL, tspan, y0);

% Calculer l'angle de volant (thetav) en fonction du temps
f = 1 / 8;
thetav = 2.2 * sin(2 * pi * f * (t - 0.25));
thetav(t < 0.25 | t > 8.25) = 0;


figure(1)
subplot(2,2,1)
plot(t, thetav ,'g')
title('Angle volant')
xlabel('Temps')
ylabel('Angle volant')

subplot(2,2,2)
plot(t, y(:, 4), 'magenta')
title('Yg (position latérale)')
xlabel('Temps')
ylabel('Position latérale')

subplot(2,2,3)
plot(t, rad2deg(y(:, 1)) ,'r')
title('Vitesse de lacet')
xlabel('Temps')
ylabel('Vitesse de lacet')

Fonction_nl;
% y(:, 2) = - y(:, 2);
Model_NL(t,y);


