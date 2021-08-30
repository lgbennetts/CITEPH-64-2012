function [ lambda ] = spectralradius( u )

    g=9.81;
    
    if u(1)<0
        u(1)=0;
    end
    
    if isnan(u(2)/u(1))==true || isinf(u(2)/u(1))==true
        lambda=[0,0];
    else
        lambda=[u(2)/u(1)+sqrt(g*u(1)),u(2)/u(1)-sqrt(g*u(1))];
    end


end