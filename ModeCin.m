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