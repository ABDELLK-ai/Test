function dY = ModeLinCir(t,Y,Mt,Lf,Lr,c_yf,c_yr,lmbd,Vx0,Dx,T_t,A)

L=Lr+Lf;
e=Mt*((Lr*c_yr-Lf*c_yf)/(2*L*c_yr*c_yf))+(L/(Vx0*Vx0));

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        dY(1,i)= 0;
        dY(2,i)= Y(1,i);
    else
        dY(1,i)= beta(t(i),Vx0,Dx,T_t,A)/(lmbd*e);
        dY(2,i)= Y(1,i);
    end
end
end
