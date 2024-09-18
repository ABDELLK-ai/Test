function dydt = Model_NL(t, y)
    % y(1) : phi. = vitesse de lacet
    % y(2) : vy = vitesse latérale
    % y(3) : psi = angle de lacet
    % y(4) : Yg = position latérale

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
    Fonction_nl;
    l = 1; % Largeur du véhicule
   

    % Convertir la vitesse longitudinale de km/h en m/s
    Vxt = 90 / 3.6;
    Vx0 = Vxt;

    % Calculer l'angle de volant (thetav) en fonction du temps
    f = 1 / 8;
    thetav = 2.2 * sin(2 * pi * f * (t - 0.25));
    thetav(t < 0.25 | t > 8.25) = 0;

    % Calculer l'angle de braquage des roues avant (Bvt) en fonction de thetav
    lambda = 16;
    Bvt = thetav / lambda;
    Bvt = deg2rad(Bvt);
    
    % y(1) : phi. = vitesse de lacet
    % y(2) : vy = vitesse latérale
    % y(3) : psi = angle de lacet
    % y(4) : Yg = position latérale    
    
%     y(1) = rad2deg(y(1));
    % Calculer les angles de dérive des pneumatiques
    delta11 = Bvt - atan((y(2) + Lf * y(1)) / (Vx0 - l * y(1)));
    delta12 = Bvt - atan((y(2) + Lf * y(1)) / (Vx0 + l * y(1)));
    delta21 = - atan((y(2) - Lr * y(1)) / (Vx0 - l * y(1)));
    delta22 = - atan((y(2) - Lr * y(1)) / (Vx0 + l * y(1)));

    
    figure(2)
    subplot(2,2,1)
    plot(t, delta11);
    subplot(2,2,2)
    plot(t, delta12);
    subplot(2,2,3)
    plot(t, delta21);
    subplot(2,2,4)
    plot(t, delta22);

    
    % Calculer les efforts latéraux des pneumatiques
    Fy11 = mu * D * sin(C * atan(B * (1 - E) * delta11 + (E / B) * atan(B * delta11)));
    Fy12 = mu * D * sin(C * atan(B * (1 - E) * delta12 + (E / B) * atan(B * delta12)));
    Fy21 = mu * D * sin(C * atan(B * (1 - E) * delta21 + (E / B) * atan(B * delta21)));
    Fy22 = mu * D * sin(C * atan(B * (1 - E) * delta22 + (E / B) * atan(B * delta22)));
%     Fy21 = 0.32 * Fy11 %mu * D * sin(C * atan(B * (1 - E) * delta21 + (E / B) * atan(B * delta21)));
%     Fy22 = 0.32 * Fy12 %mu * D * sin(C * atan(B * (1 - E) * delta22 + (E / B) * atan(B * delta22)));
    
%     Fy = Fy11 + Fy12 + Fy21 + Fy22 + f0y;
%     alat_n = Fy / Mt;

%     % Tracer alat_n
%     figure(1)
%     subplot(2,2,4)
%     plot(t, alat_n);
%     xlabel('Temps');
%     ylabel('Acceleration latérale');
%     title('Evolution de l''accélération latérale');
%     
%     % Tracer alat_n
%     figure(2)
%     subplot(2,2,1)
%     plot(t, Fy11);
%     subplot(2,2,2)
%     plot(t, Fy12);
%     subplot(2,2,3)
%     plot(t, Fy21);
%     subplot(2,2,4)
%     plot(t, Fy22);


    % Calculer la somme des couples extérieurs
    Cz = Lf * (Fy11 + Fy12) - Lr * (Fy21 + Fy22);
    
    % y(1) : phi. = vitesse de lacet
    % y(2) : vy = vitesse latérale
    % y(3) : psi = angle de lacet
    % y(4) : Yg = position latérale

    % Calculer les dérivées des variables d'état
    dydt1 = Cz /Iz;
    dydt2 = (Fy11 + Fy12 + Fy21 + Fy22 + f0y) / Mt - (Vx0 * y(1));
    dydt3 = y(1);
    dydt4 = Vx0 * sin(y(3)) + y(2) * cos(y(3));
    

    % Retourner les dérivées dans un vecteur colonne
    dydt = [dydt1; dydt2; dydt3; dydt4];

end
