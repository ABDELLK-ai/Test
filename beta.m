function Y= beta(t,Vx0,Dx,T_t,A)

for i=1:length(t)
    if t(i)< (5/Vx0) || t(i)> ((Dx+5)/Vx0)
        Y(i)=0;
    else
        Y(i)= (pi/180)*A*sin(2*pi*(t(i)-(5/Vx0))/T_t);
    end
end
end
