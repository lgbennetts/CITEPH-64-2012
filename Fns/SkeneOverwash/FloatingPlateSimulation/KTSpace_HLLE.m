function [ spacediscrete,amax ] = KTSpace_HLLE( u,dx,x,theta,mindepth )

%     u(1,:)=u(1,:)*pi*sqrt(2/24)
%     u(2,:)=u(2,:)*sqrt(29/24)*sqrt(pi*sqrt(2/24))

   g=9.81;

    %theta=1;

    ux=u*0;
    
    for i=1:length(u(:,1))
        ux(i,:)=minmod(u(i,:),dx,theta);
    end
   
    uplus=u*0;
    uneg=u*0;
    
    for i=2:length(u(1,:))
        uplus(:,i)=u(:,i)-1/2*ux(:,i);
        uneg(:,i)=u(:,i-1)+1/2*ux(:,i-1);
    end
    
    uneg(:,i)=u(:,i-1)+1/2*ux(:,i-1);
    uplus(:,i)=u(:,i);
    i=1;
    uplus(:,i)=u(:,i)-1/2*ux(:,i);
    uneg(:,i)=u(:,i);
    
    aneg=u(1,:)*0;
    aplus=aneg;
    
    for i=1:length(uneg)
        a=spectralradius(uplus(:,i));
        b=spectralradius(uneg(:,i));
        aneg(i)=min([a,b,0]);
        aplus(i)=max([a,b,0]);
    end
            
    amax=max([abs(aneg),abs(aplus)]);
    
    
    spacediscrete=u*0;
    
    
    i=2;
    
    hplus=(aplus(i)*f(uneg(:,i))-aneg(i)*f(uplus(:,i)))/(aplus(i)-aneg(i))...
                +(aplus(i)*aneg(i))/(aplus(i)-aneg(i))*(uplus(:,i)-uneg(:,i));
    hneg=(aplus(i-1)*f(uneg(:,i-1))-aneg(i-1)*f(uplus(:,i-1)))/(aplus(i-1)-aneg(i-1))...
                +(aplus(i-1)*aneg(i-1))/(aplus(i-1)-aneg(i-1))*(uplus(:,i-1)-uneg(:,i-1));
            
    if aplus(i-1)-aneg(i-1)==0
        hplus=[0;0];
    end
    if aplus(i)-aneg(i)==0
        hneg=[0;0];
    end
    spacediscrete(:,i-1)=-(hplus-hneg)/dx;
            
    for i=3:length(u)
        
        hneg=hplus;
        if aplus(i) - aneg(i)==0
            hplus=[0;0];
        else
            hplus=(aplus(i)*f(uneg(:,i))-aneg(i)*f(uplus(:,i)))/(aplus(i)-aneg(i))...
                +(aplus(i)*aneg(i))/(aplus(i)-aneg(i))*(uplus(:,i)-uneg(:,i));
        end
            
        spacediscrete(:,i-1)=-(hplus-hneg)/dx;
        
    end
    
    %spacediscrete=spacediscrete*0;
    
    spacediscrete(:,1)=[0;0];
    
    fplus=[0;0];
    fneg=[0;0];
    for i=[1,2]%,length(x)-2,length(x)-1]
        fneg=fplus;   
        
        hl=u(1,i);
        cl=sqrt(abs(hl*g));
        hr=u(1,i+1);
        cr=sqrt(abs(hr*g));
        
        if hl>mindepth
           ul=u(2,i)/u(1,i); 
           ubarl=u(2,i)/sqrt(u(1,i));
        else
           ul=0;
           ubarl=0;
        end
        if hr>mindepth
            ur=u(2,i+1)/u(1,i+1);
            ubarr=u(2,i+1)/sqrt(u(1,i+1));
        else
            ur=0;
            ubarr=0;
        end
        
        cbar=sqrt(abs(g/2*(hl+hr)));
        
        if hl>0 && hr>0
            ubar=(ubarr+ubarl)/(sqrt(hl)+sqrt(hr));
        else
            ubar=0;
        end
        
        a1=ubar+cbar;
        a2=ubar-cbar;
        
        if cbar~=0
        alpha1=(u(2,i+1)-u(2,i)+(-ubar+cbar)*(u(1,i+1)-u(1,i)))/(2*cbar);
        alpha2=(u(2,i+1)-u(2,i)+(-ubar-cbar)*(u(1,i+1)-u(1,i)))/(-2*cbar);
        else
            alpha1=0;
            alpha2=0;
        end
        
        br=max(ur+cr,a1);
        bl=min(ul-cl,a2);
        
        bplus=max(br,0);
        bneg=min(bl,0);
        
        e1=[1;ubar+cbar];
        e2=[1;ubar-cbar];
        
        if bplus > bneg
            p1=((bplus+bneg)*a1-2*bplus*bneg)/(bplus-bneg);
            p2=((bplus+bneg)*a2-2*bplus*bneg)/(bplus-bneg);
        else
            p1=0;
            p2=0;
        end
        
        fplus=1/2*(f(u(:,i+1))+f(u(:,i))-p1*alpha1*e1-p2*alpha2*e2);
        
         if i==2 || i==length(x)-1
            spacediscrete(:,i)=-(fplus-fneg)/dx;
         end
         amax=max([amax,abs(br),abs(bl)]);
         if sum(sum(isnan(spacediscrete)))>1
            %keyboard 
         end  
    end  

end

