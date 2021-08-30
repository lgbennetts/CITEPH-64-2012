function [ ux ] = minmod( u,dx,theta )

    ux=u*0;
    
    for i=2:length(u)-1
        
        term1=(u(i)-u(i-1))/1*theta;
        term2=(u(i+1)-u(i-1))/2/1;
        term3=(u(i+1)-u(i))/1*theta;
    
        if term1>0 && term2>0 && term3>0
            ux(i)=min(term1,min(term2,term3));
        elseif term1<0 && term2<0 && term3<0
            ux(i)=max(term1,max(term2,term3));
        else
            ux(i)=0;
        end
        
    end
    
    i=length(ux);
    ux(i)=(u(i)-u(i-1))/1*theta;
    
    ux(1)=(u(2)-u(1))/1*theta;
    

end

