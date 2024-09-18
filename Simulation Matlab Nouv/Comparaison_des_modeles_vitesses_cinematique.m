clear all,close all,clc
Fig=1;
Mt=1759; Iz= 2638.5; Lf=0.71; Lr= 2.13; c_yf=94446; c_yr=48699; lmbd= 16; Dx=200;L=Lr+Lf;
Vx0=[10 30 50 70 90 110 130]/3.6;
T_t=Dx./Vx0; 
Amplitude=[3.5/2.42 3.5/2.30 3.5/2.10 3.5/1.84 3.5/1.58 3.5/1.34 3.5/1.14];

y3ode=[0;0];
opts = odeset('RelTol',1e-10);
[t1ode,Y1ode]=ode23(@(t,y) ModeCin(t,y,Vx0(1), Dx,T_t(1),L,Amplitude(1)),[0 100],y3ode,opts);
[t2ode,Y2ode]=ode23(@(t,y) ModeCin(t,y,Vx0(2), Dx,T_t(2),L,Amplitude(2)),[0 100],y3ode,opts);
[t3ode,Y3ode]=ode23(@(t,y) ModeCin(t,y,Vx0(3), Dx,T_t(3),L,Amplitude(3)),[0 100],y3ode,opts);
[t4ode,Y4ode]=ode23(@(t,y) ModeCin(t,y,Vx0(4), Dx,T_t(4),L,Amplitude(4)),[0 100],y3ode,opts);
[t5ode,Y5ode]=ode23(@(t,y) ModeCin(t,y,Vx0(5), Dx,T_t(5),L,Amplitude(5)),[0 100],y3ode,opts);
[t6ode,Y6ode]=ode23(@(t,y) ModeCin(t,y,Vx0(6), Dx,T_t(6),L,Amplitude(6)),[0 100],y3ode,opts);
[t7ode,Y7ode]=ode23(@(t,y) ModeCin(t,y,Vx0(7), Dx,T_t(7),L,Amplitude(7)),[0 100],y3ode,opts);

for i=1:length(t1ode)
    dxi1(i)= (Vx0(1)/L)*tan(beta(t1ode(i),Vx0(1),Dx,T_t(1),Amplitude(1))/16);
end
for i=1:length(t2ode)
    dxi2(i)= (Vx0(2)/L)*tan(beta(t2ode(i),Vx0(2),Dx,T_t(2),Amplitude(2))/16);
end
for i=1:length(t3ode)
    dxi3(i)= (Vx0(3)/L)*tan(beta(t3ode(i),Vx0(3),Dx,T_t(3),Amplitude(3))/16);
end
for i=1:length(t4ode)
    dxi4(i)= (Vx0(4)/L)*tan(beta(t4ode(i),Vx0(4),Dx,T_t(4),Amplitude(4))/16);
end
for i=1:length(t5ode)
    dxi5(i)= (Vx0(5)/L)*tan(beta(t5ode(i),Vx0(5),Dx,T_t(5),Amplitude(5))/16);
end
for i=1:length(t6ode)
    dxi6(i)= (Vx0(6)/L)*tan(beta(t6ode(i),Vx0(6),Dx,T_t(6),Amplitude(6))/16);
end
for i=1:length(t7ode)
    dxi7(i)= (Vx0(7)/L)*tan(beta(t7ode(i),Vx0(7),Dx,T_t(7),Amplitude(7))/16);
end

% % Affichage à 90 km/h
% figure(Fig)
% Fig=Fig+1;
% plot(t5ode,(180/pi)*Y5ode(:,1))
% title('Angle de lacet')
% xlabel('t(s)')
% ylabel('Xi(°)')
% axis([0 10 0 3.5])
% grid on
% 
% figure(Fig)
% Fig=Fig+1;
% plot(t5ode,Y5ode(:,2))
% title('Position latérale')
% xlabel('t(s)')
% ylabel('Y_G(m)')
% axis([0 10 0 6])
% grid on
% 
% figure(Fig)
% Fig=Fig+1;
% plot(t5ode,(180/pi)*dxi5)
% title('Vitesse de lacet')
% xlabel('t(s)')
% ylabel('dXi(°/s)')
% axis([0 10 -1.5 1.5])
% grid on

%% Affichage pour différentes vitesses
figure(Fig)
Fig=Fig+1;
plot(t1ode*Vx0(1),(180/pi)*dxi1)
hold on
plot(t2ode*Vx0(2),(180/pi)*dxi2)
hold on
plot(t3ode*Vx0(3),(180/pi)*dxi3)
hold on
plot(t4ode*Vx0(4),(180/pi)*dxi4)
hold on
plot(t5ode*Vx0(5),(180/pi)*dxi5)
hold on
plot(t6ode*Vx0(6),(180/pi)*dxi6)
hold on
plot(t7ode*Vx0(7),(180/pi)*dxi7)
hold off
title('Modèle cinématique : Vitesse de lacet')
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
xlabel('X_G(m)')
ylabel('d\Psi(m)')
axis([0 210 -2.5 2.5])
grid on

figure(Fig)
Fig=Fig+1;
plot(t1ode*Vx0(1),Y1ode(:,2))
hold on
plot(t2ode*Vx0(2),Y2ode(:,2))
hold on
plot(t3ode*Vx0(3),Y3ode(:,2))
hold on
plot(t4ode*Vx0(4),Y4ode(:,2))
hold on
plot(t5ode*Vx0(5),Y5ode(:,2))
hold on
plot(t6ode*Vx0(6),Y6ode(:,2))
hold on
plot(t7ode*Vx0(7),Y7ode(:,2))
hold off
title('Modèle cinématique : Position latérale')
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
xlabel('X_G(m)')
ylabel('Y_G(m)')
axis([0 210 0 8])
grid on

%% Partie fréquentielle

S11=tf([Vx0(1)^2],[lmbd*L 0 0]);
S12=tf([Vx0(2)^2],[lmbd*L 0 0]);
S13=tf([Vx0(3)^2],[lmbd*L 0 0]);
S14=tf([Vx0(4)^2],[lmbd*L 0 0]);
S15=tf([Vx0(5)^2],[lmbd*L 0 0]);
S16=tf([Vx0(6)^2],[lmbd*L 0 0]);
S17=tf([Vx0(7)^2],[lmbd*L 0 0]);

S21=tf([Vx0(1)],[lmbd*L]);
S22=tf([Vx0(2)],[lmbd*L]);
S23=tf([Vx0(3)],[lmbd*L]);
S24=tf([Vx0(4)],[lmbd*L]);
S25=tf([Vx0(5)],[lmbd*L]);
S26=tf([Vx0(6)],[lmbd*L]);
S27=tf([Vx0(7)],[lmbd*L]);

figure(Fig)
Fig=Fig+1;
bode(S11,S12,S13,S14,S15,S16,S17,{10^-4,10^4})
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title('Modèle cinématique : Bode entre en la position latérale et l''angle du volant')
grid on

figure(Fig)
Fig=Fig+1;
bode(S21,S22,S23,S24,S25,S26,S27,{10^-4,10^4})
legend('10 km/h','30 km/h','50 km/h','70 km/h','90 km/h','110 km/h','130 km/h')
title('Modèle cinématique : Bode entre en la vitesse de lacet et l''angle du volant')
grid on

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






