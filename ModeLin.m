function dY = ModeLin(t,Y,Mt, Iz, Lf, Lr, c_yf, c_yr, lmbd, Vx0, Dx,T_t,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i)  -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    else
        dY(1,i)= Y(2,i);
        dY(2,i)= -2*((Lf^2*c_yf + Lr^2*c_yr)/(Vx0*Iz))*Y(2,i) + 2*((-Lf*c_yf  +Lr*c_yr)/(Vx0*Iz))*Y(3,i) +2*Lf*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lmbd*Iz);
        dY(3,i)= ((-2*(Lf*c_yf-Lr*c_yr)/(Mt*Vx0))-Vx0)*Y(2,i) -2*((c_yf+c_yr)/(Vx0*Mt))*Y(3,i) + 2*c_yf*beta(t(i),Vx0,Dx,T_t,A)/(lmbd*Mt);
        dY(4,i)= Vx0*Y(1,i) + Y(3,i);

    end
end
end
