function [ result ] = f( uj )

    g=9.81;

    if uj(1)==0
        result=[0;0];
    else
        result=[uj(2);1/2*g*uj(1)^2+uj(2)^2/uj(1)];
    end

end

