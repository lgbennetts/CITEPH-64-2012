close all;
clear all;

NLB = 200;
th_res=100;
SURGE=0;
terms_grn=100;
extra_pts=[];     
rigid = 4; 
Param = ParamDef_Oceanide(rigid); 
Param = ModParam_def(Param,NLB,NLB,extra_pts,terms_grn,th_res); 

PhysVars.omega = 2*pi/0.95;
PhysVars.g = Param.g;
PhysVars.h = Param.bed;
PhysVars.rho0= Param.rho_0;
PhysVars.rho = Param.rho;
PhysVars.thickness= Param.thickness;
PhysVars.E = Param.E;
PhysVars.nu = Param.nu;
PhysVars.L = Param.floe_diam/2;

ModelVars.RootFindThreshold = 10^-12; 
ModelVars.nBeamTerms = 200; 
ModelVars.nWaterTerms = 10^5; 
% ModelVars.nBeamLocations = 100;  

nS = 6:22;
nS = 6:22;
% nS = 20:22;
dxs =  zeros(size(nS));
rPBCWWA = zeros(size(nS));
rPBCWWN = zeros(size(nS));
rFE2WW = zeros(size(nS));
rFE3WW = zeros(size(nS));

rPBCLB = zeros(size(nS));
rFE2LB = zeros(size(nS));
rFE3LB = zeros(size(nS));

D = PhysVars.E*PhysVars.thickness^2/(12*PhysVars.rho*(1-PhysVars.nu^2));
alpha = PhysVars.omega^2/PhysVars.g;
beta = (PhysVars.rho*PhysVars.thickness*D) / (PhysVars.rho0*PhysVars.g);
gamma = PhysVars.rho*PhysVars.thickness/(PhysVars.rho0);

for i = 1: length(nS)
    
    N = round(2^(nS(i)/2));
    NLB = round(2^(nS(i)/3));
    Param = ModParam_def(Param,NLB,NLB,extra_pts,terms_grn,th_res); 
%     Sol_Vec = fn_Floating_Elastic_Plate_Sol_FixedExpansions(PhysVars,ModelVars);
%     [xi,displacement,potential,R,T,f]
    
    [xi,displacement,potential,R,T,f,ab,lambda,diff,D4Disp] = elastic_plate_modes(alpha,beta,gamma,PhysVars.h,PhysVars.L,ModelVars.nBeamTerms,N);
    x=-PhysVars.L:2*PhysVars.L/(length(displacement)-1):PhysVars.L;
    Disp = 1i/sqrt(alpha)*displacement;
    D4DispR = 1i/sqrt(alpha)*D4Disp;
    Pot = 1i*PhysVars.g/PhysVars.omega*potential;
%     Disp = 1i*displacement;
%     D4DispR = 1i*D4Disp;
%     Pot = 1i*sqrt(PhysVars.g)*potential;
%      Disp = displacement;
%      D4DispR = D4Disp;
%      Pot = potential;
    
    
%     out = fn_ElasticRaft2d('freq',PhysVars.omega/(2*pi),Param,'disk',SURGE,0,1,1,0,'r',x);
%     TN = (out(1).value);
%     RN = (out(2).value);
%     EtaOut = out(3).value;
%     PhiOut = (1i*PhysVars.g/PhysVars.omega)*out(4).value;
    % {x,xi,displacement,potential,R,T,f,ab,lambda,diff}
%     x = Sol_Vec{1};
%     Displacement = Sol_Vec{3};
%     Potential = Sol_Vec{4};
% 

    
    dx = x(2) - x(1);
    D =  PhysVars.E*PhysVars.thickness^3/(12*(1-PhysVars.nu^2));
    [outPBCAD4WW,outPBCND4WW,outFE2WW,outFE3WW] = Measure(1i*PhysVars.omega*Disp,1i*PhysVars.omega*D4DispR,Pot,D,PhysVars.rho,PhysVars.rho0,PhysVars.thickness,PhysVars.omega,PhysVars.g,dx,PhysVars.L,x);
%     [outPBCLB,outFE2LB,outFE3LB] = Measure(1i*PhysVars.omega*EtaOut, PhiOut,D,PhysVars.rho,PhysVars.rho0,PhysVars.thickness,PhysVars.omega,PhysVars.g,dx,PhysVars.L,x);

    rPBCWWA(i) = outPBCAD4WW;
    rPBCWWN(i) = outPBCND4WW;
    rFE2WW(i) = outFE2WW;
    rFE3WW(i) = outFE3WW;
    xs{i} = x;
%     rT1WW{i} = T1WW;
%     rT2WW{i} = T2WW;
%     D4Disps{i} = D4Disp;
%     rPBCLB(i) = outPBCLB;
%     rFE2LB(i) = outFE2LB;
%     rFE3LB(i) = outFE3LB;
    dxs(i) = dx;
    
end


% figure('DefaultAxesFontSize',18);
% % set(gca,'FontSize',24) 
% hold on;
% for i = 1:length(xs)
% plot(xs{i},real(rT1WW{i}),'--','MarkerSize',12,'DisplayName', ['T1 (d/dz pot) WW | dx = ',num2str(dxs(i),'%3.2e')]);
% plot(xs{i},real(rT2WW{i}),'-.','MarkerSize',12,'DisplayName', ['T2 (pot) WW | dx = ', num2str(dxs(i),'%3.2e') ]);
% % plot(xs{i},real(1i*PhysVars.omega*D4Disps{i}),':','MarkerSize',12,'DisplayName', ['D4 Disp']);
% end

% xlabel('$\Delta x$ - (Beam Locations)','Interpreter','latex');
% ylabel('d/dz Pot , Pot Terms in Free Surface BC');
% namestr = ['Root Threshold : ', num2str(ModelVars.RootFindThreshold,'%3.2e') , ' ,  Beam Modes ', num2str(ModelVars.nBeamTerms), ',  Water Modes ',num2str(ModelVars.nWaterTerms,'%3.2e')];
% paramstr = ['g: ', num2str(PhysVars.g,'%3.2f'),', $\rho$: ', num2str(PhysVars.rho,'%3.2e'),', $\rho_0$: ', num2str(PhysVars.rho0,'%3.2e'), ...
%     ', E: ', num2str(PhysVars.E,'%3.2e'), ', d: ', num2str(PhysVars.thickness,'%3.2e') , ', H: ', num2str(PhysVars.h,'%3.2f') ,...
%     ', $\nu$: ', num2str(PhysVars.nu,'%3.2f'), ', $\omega$: ', num2str(PhysVars.omega,'%3.2f')];
% title({'Plot DzPot | Pot ', namestr,paramstr},'Interpreter','latex');
% legend();


figure('DefaultAxesFontSize',18);
% set(gca,'FontSize',24) 
loglog(dxs,rPBCWWA,'.k','MarkerSize',12,'DisplayName', 'Plate Boundary WW (Analytic D4)');
hold on;
loglog(dxs,rPBCWWN,'xk','MarkerSize',12,'DisplayName', 'Plate Boundary WW (Numeric D4)');
loglog(dxs,rFE2WW,'.r','MarkerSize',12,'DisplayName', 'D2 Free Edge Condition WW');
loglog(dxs,rFE3WW,'.b','MarkerSize',12,'DisplayName', 'D3 Free Edge Condition WW');
% loglog(dxs,rPBCLB,'xk','MarkerSize',12,'DisplayName', 'Plate Boundary LB');
% loglog(dxs,rFE2LB,'xr','MarkerSize',12,'DisplayName', 'D2 Free Edge Condition LB');
% loglog(dxs,rFE3LB,'xb','MarkerSize',12,'DisplayName', 'D3 Free Edge Condition LB');
xlabel('$\Delta x$ - (Beam Locations)','Interpreter','latex');
ylabel('Error Measure');
namestr = ['Root Threshold : ', num2str(ModelVars.RootFindThreshold,'%3.2e') , ' ,  Beam Modes ', num2str(ModelVars.nBeamTerms), ',  Water Modes ',num2str(ModelVars.nWaterTerms,'%3.2e')];
paramstr = ['g: ', num2str(PhysVars.g,'%3.2f'),', $\rho$: ', num2str(PhysVars.rho,'%3.2e'),', $\rho_0$: ', num2str(PhysVars.rho0,'%3.2e'), ...
    ', E: ', num2str(PhysVars.E,'%3.2e'), ', d: ', num2str(PhysVars.thickness,'%3.2e') , ', H: ', num2str(PhysVars.h,'%3.2f') ,...
    ', $\nu$: ', num2str(PhysVars.nu,'%3.2f'), ', $\omega$: ', num2str(PhysVars.omega,'%3.2f')];
title({'Convergence Plot ', namestr,paramstr},'Interpreter','latex');
legend();





function D4 = FourthDeriv(fx,dx)
D4 = zeros(size(fx));
for i = 1:2
D4(i) = (3*fx(i)-14*fx(i+1)+26*fx(i+2)-24*fx(i+3)+11*fx(i+4) -2*fx(i+5)) / dx^4;
D4(end + 1 -i) = (3*fx(end + 1 -i)-14*fx(end + 1 -(i+1))+26*fx(end + 1 -(i+2))-24*fx(end + 1 - (i+3))+11*fx(end + 1- (i+4)) - 2*fx(end + 1- (i+5))) / dx^4;
end
%
for i = 3:length(fx)-2
    D4(i) = (fx(i-2)-4*fx(i-1)+6*fx(i)-4*fx(i+1)+fx(i+2)) / dx^4;
end
end

function D2 = SecondDeriv(fx,dx)
D2 = zeros(size(fx));
for i = 1:1
%     2 	−5 	4 	−1
D2(i) = (2*fx(i)-5*fx(i+1)+4*fx(i+2)-1*fx(i+3)) / dx^2;
D2(end + 1 - i) = (2*fx(end + 1 -i)-5*fx(end + 1 -(i+1))+4*fx(end + 1 -(i+2))-1*fx(end + 1 -(i+3))) / dx^2;
end
%
for i = 2:length(fx)-1
    D2(i) = (fx(i-1)-2*fx(i)+ fx(i+1)) / dx^2;
end
end

function D3 = ThirdDeriv(fx,dx)
D3 = zeros(size(fx));
for i = 1:2
% −5/2 	9 	−12 	7 	−3/2 
D3(i) = (-5/2*fx(i)+9*fx(i+1)-12*fx(i+2)+7*fx(i+3)-3/2*fx(i+4)) / dx^3;
D3(end + 1 -i) = (5/2*fx(end + 1 -i)-9*fx(end + 1 - (i+1))+12*fx(end + 1 - (i+2))-7*fx(end + 1 - (i+3))+3/2*fx(end + 1 -(i+4))) / dx^3;
end
%−1/2 	1 	0 	−1 	1/2
for i = 3:length(fx)-2
    D3(i) = (-1/2*fx(i-2)+1*fx(i-1)-1*fx(i+1)+1/2*fx(i+2)) / dx^3;
end
end




function [LHSA,LHSN] =  PlateBoundary(DzPotential,Dx4DzPotential,Potential,D,rho,rho0,thickness,omega,g,dx)
 
D4Num = FourthDeriv(DzPotential,dx);
D4Expr = Dx4DzPotential;
% figure();
% plot(real(D4Num));
% hold on;
% plot(real(D4Expr));
% DzTerms = (D*FourthDeriv(DzPotential,dx)); %+ (rho0*g- rho*thickness*omega^2)*DzPotential;

alpha = omega^2/g;
%  D =  PhysVars.E*PhysVars.thickness^3/(12*(1-PhysVars.nu^2));
beta = (D) / (rho0*g);
gamma = rho*thickness/(rho0);

% DzTermsA = D*D4Expr + (rho0*g- rho*thickness*omega^2)*DzPotential;
% DzTermsN = D*D4Num  + (rho0*g- rho*thickness*omega^2)*DzPotential;
% Rest = rho0*omega^2*Potential;

DzTermsA = beta*D4Expr - (alpha*gamma -1)*DzPotential;
DzTermsN = beta*D4Num  - (alpha*gamma -1)*DzPotential;
Rest = alpha*Potential;


figure()
plot(DzTermsA,'-b','DisplayName','analytic d/dz pot')
hold on;
plot(DzTermsN,'-r','DisplayName','Numerical d/dz pot')
plot(Rest,'--k','DisplayName','pot')
xlabel('Index')
ylabel('Function value at index')
paramstr = ['D: ', num2str(D,'%3.2e'),', $rho$: ', num2str(rho,'%3.2e'),', rho0: ', num2str(rho0,'%3.2e'), ...
    ', thickness: ', num2str(thickness,'%3.2e'), ', omega: ', num2str(omega,'%3.2e') , ', g: ', num2str(g,'%3.2f') ,...
    ', dx: ', num2str(dx,'%3.2e')];
title(paramstr)
legend();

LHSA =  (DzTermsA - Rest);
LHSN =  (DzTermsN - Rest);
end

function [outPBCAD4,outPBCND4,outFE2,outFE3] = Measure(DzPotential,Dx4DzPotential,Potential,D,rho,rho0,thickness,omega,g,dx,L,x)
 [LHSA,LHSN] =  PlateBoundary(real(DzPotential),real(Dx4DzPotential),real(Potential),D,rho,rho0,thickness,omega,g,dx);
 outPBCAD4 = trapz(x,abs(LHSA))./(2*L);
 outPBCND4= trapz(x,abs(LHSN))./(2*L);
%  outPBC = outPBC / (rho0*omega^2);
 D2 = SecondDeriv(DzPotential,dx);
 outFE2 =  abs(real(D2(1) + D2(end))); 
 D3 = ThirdDeriv(DzPotential,dx);
 outFE3 =  abs(real(D3(1) + D3(end))); 
%  Mean = mean(LHS);
end