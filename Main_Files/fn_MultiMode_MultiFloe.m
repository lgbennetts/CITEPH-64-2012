function [Rm,Tm,Rp,Tp,v_vec,u_vec,k0,wt_0,xNm] = ...  %displ_fs, floe_displ, 
    fn_MultiMode_MultiFloe(...
    PVec, Vert_Dim, evs, Geom_Vec, kappa, beta_vec, thick_vec, ...
    draft_vec, R_vec, posits, rad_vec, th_vec, x_vec, y_vec, FS_mesh, ...
    res, extra_pts)

disp(['Vertical modes = ' int2str(Vert_Dim)])

% - Geom_Vec = [scaling No.plates depth thickness l w]
% - Param_Vec = [g rho_water draft beta gamma rho_ice nu E D];

% - 27.09.10 created by LGB
% - adapted from fn_SingleMode_MultiFloe_v2

Tol_vec(1) = 1e-16; % - Real root error - %
Tol_vec(2) = 1e-16; % - Imag root error (init) - %
Tol_vec(3) = 1e-1; % - Az_Dim tol - %

scat = 1; % - scat=0 (scat only) scat=1 (scat+inc)

%k0 = 2*pi/lam0;

%%% Number of floes %%%

Np = length(thick_vec); disp(['Plates = ' int2str(Np)])

%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tank Dimensions %%%

width = Geom_Vec(2); bed = Geom_Vec(3); lth = Geom_Vec(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter_vector = [PVec, kappa, bed];

%%% Roots in the free surface %%%

k0 = zeros(1,Vert_Dim); wt_0 = k0;

for loop_Dim=1:Vert_Dim
    k0(loop_Dim) = GetRootsMMA_FS_PWC(parameter_vector, loop_Dim, Tol_vec);
    wt_0(loop_Dim) = weight_0_PWC(parameter_vector(7), k0(loop_Dim));   
end

%%% Roots in the plates %%%

kk = zeros(Vert_Dim,Np); wt = kk; %wt_0 = wt;

mu_0 = zeros(1,Np); mu_1 = mu_0;

for loop_P=1:Np
    pv_plate = [parameter_vector,draft_vec(loop_P),thick_vec(loop_P),...
        bed-draft_vec(loop_P),draft_vec(loop_P),beta_vec(loop_P)];
 for loop_Dim = 1:Vert_Dim
    kk(loop_Dim,loop_P) = GetRootsMMA_PWC(pv_plate, loop_Dim, Tol_vec);
    wt(loop_Dim,loop_P) = weight_PWC(pv_plate, kk(loop_Dim,loop_P));
 end
    [mu_0(loop_P), mu_1(loop_P)] = mu_new_PWC(pv_plate, Vert_Dim, ...
        kk(:,loop_P), wt(:,loop_P));
    clear pv_plate
end

%%%% - Amplitude of the x-derivative !!!!! - %%%%%
Amp0 = 1;

%%%%%%%%%%%%%%%%%%%%%%%

%%% Angular dimension %%%

Az_Dim_vec = 2*Tol_vec(3)*ones(3,1);
Az_Dim = 0;
count = 0;
while max(Az_Dim_vec)>Tol_vec(3)    
 Az_Dim_vec(count+1) = abs(besselj(Az_Dim,k0(1)*max(R_vec)));
 Az_Dim=Az_Dim+1; count = count+1; count = mod(count,length(Az_Dim_vec));
end

clear Az_Dim_vec count

Az_Dim=Az_Dim-2;

disp(['Azimuthal modes = ' int2str(Az_Dim) ' (x2 + 1)'])
disp(['Extra points = ' int2str(extra_pts)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Greens functions evauated around the plates %%%
%%% - PhiInc = sum_{m=Az_Dims} J_m*PhiInc_{m}*exp{i*m*theta}
%%% - PhiInc_dr = sum_{m=Az_Dims} PhiInc_dr_{m}*exp{i*m*theta}

[PhiInc, G, Gdr_add, v_vec, u_vec, theta_inc, ExtraOut] = ... 
    fn_IntEqn_Reson(Vert_Dim, Az_Dim, k0, R_vec, width, lth, posits, ...
    res, Amp0, evs, extra_pts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set Y-DIM %%%
Y_Dim=length(u_vec(1,:)); Y0_Dim=Y_Dim-evs; new_res = ExtraOut{2};
%%%%%%%%%%%%%%%%%

disp(['Tank modes = ',num2str(Y0_Dim)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dirichlet 2 Neumann on plate boundaries %%%
%%% - Matrices are block tri-diag -

PsiMat = zeros(Np*Vert_Dim*(2*Az_Dim+1)); PsidrMat = PsiMat;
AW0 = zeros(Np*Vert_Dim*(2*Az_Dim+1)); AW = AW0; VT0 = AW; VT = AW;
Disp_mat = zeros(Vert_Dim*(2*Az_Dim+1),length(rad_vec),Np);
Amp_mu_Mat = zeros(2,Vert_Dim,2*Az_Dim+1,Np);

v1 = 1:Vert_Dim*(2*Az_Dim+1);

for loop_P=1:Np
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% NO SCATTERER TEST !!!!
%  for loop_N=1:Vert_Dim
%   v1a = v1(loop_N:Vert_Dim:Vert_Dim*(2*Az_Dim+1)); 
%   PsiMat(v1a,v1a) = diag(besselj(-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P))); 
%   PsidrMat(v1a,v1a) = k0(loop_N)*diag(Bessel_dz(@besselj,-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P)));
%  end
%  AW0(v1,v1) = eye(Vert_Dim*(2*Az_Dim+1)); AW(v1,v1) = AW0(v1,v1); 
%  VT0(v1,v1) = AW0(v1,v1); VT(v1,v1) = AW0(v1,v1);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% SOUND SOFT CYLINDER
%  for loop_N=1:Vert_Dim
%   v1a = v1(loop_N:Vert_Dim:Vert_Dim*(2*Az_Dim+1)); 
%   PsiMat(v1a,v1a) = diag(besselj(-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P))); 
%   PsidrMat(v1a,v1a) = k0(loop_N)*diag(Bessel_dz(@besselj,-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P)));
%  end
%  AW0(v1,v1) = eye(Vert_Dim*(2*Az_Dim+1)); AW(v1,v1) = AW0(v1,v1); 
%  VT0(v1,v1) = eye(Vert_Dim*(2*Az_Dim+1)); VT(v1,v1) = zeros(Vert_Dim*(2*Az_Dim+1));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% SURFACE PIERCING CYLINDER
 for loop_N=1:Vert_Dim
    v1a = v1(loop_N:Vert_Dim:Vert_Dim*(2*Az_Dim+1)); 
    PsiMat(v1a,v1a) = diag(besselj(-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P))); 
    PsidrMat(v1a,v1a) = k0(loop_N)*diag(Bessel_dz(@besselj,-Az_Dim:Az_Dim,k0(loop_N)*R_vec(loop_P)));
 end
 AW(v1,v1) = eye(Vert_Dim*(2*Az_Dim+1)); AW0(v1,v1) = zeros(Vert_Dim*(2*Az_Dim+1));
 VT0(v1,v1) = eye(Vert_Dim*(2*Az_Dim+1)); VT(v1,v1) = VT0(v1,v1);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% ELASTIC PLATE
%  pv_plate = [parameter_vector,PVec(loop_P,3),thick_vec(loop_P),...
%         bed-PVec(loop_P,3),PVec(loop_P,5),PVec(loop_P,4)];
%  [PsiMat(v1,v1),PsidrMat(v1,v1),Disp_mat(:,:,loop_P),...
%      Amp_mu_Mat(:,:,:,loop_P)] = ...
%     fn_SingleFloe(pv_plate, Vert_Dim, Az_Dim, kk(:,loop_P), ...
%      wt(:,loop_P), mu_0(loop_P), mu_1(loop_P), R_vec(loop_P), rad_vec(loop_P,:));
%  [AW0(v1,v1),AW(v1,v1),VT0(v1,v1),VT(v1,v1)] = fn_JumpMats(pv_plate,...
%     Vert_Dim, Vert_Dim, Az_Dim, k0, kk(:,loop_P), wt_0, wt(:,loop_P));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 v1=v1+Vert_Dim*(2*Az_Dim+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(ExtraOut)>2
 ExtraIn = {ExtraOut{3:end}};
else
 ExtraIn = [];
end
clear ExtraOut

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOLVE THE INTEGRAL SYSTEM %%%

%%% - Phi=sum Phi_m*exp(im*theta)
%%% - Psi=sum hat{J}_(k0*R)*Psi_m*exp(im*theta)
%%% - (d/dr)Phi=sum Phi_dr_m*exp(im*theta)
%%% - (d/dr)Psi=sum Psi_dr_m*exp(im*theta)

[Phi, Psi, Phi_dr, Psi_dr, ExtraOut] = ...
    fn_Solve(Vert_Dim, Az_Dim, k0, R_vec, width, PhiInc, G, Gdr_add, ...
    PsiMat, PsidrMat, AW0, AW, VT0, VT, extra_pts, ExtraIn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(1); plot(k0,abs(Phi-PhiInc),'.')
%figure(2); plot(k0,abs(Phi_dr-PhiInc_dr),'.')

%[abs(Phi-PhiInc), abs(Psi-PhiInc), abs(Phi_dr-PhiInc_dr), abs(Psi_dr-PhiInc_dr)]

%% -- Scattered amplitudes -- %%

if isempty(ExtraOut)==1
 ExtraIn = [];
else
 ExtraIn = [ExtraIn{1}{1}(1), ExtraOut];
 % - ExtraIn = [skip, const]
end
clear ExtraOut

%% - Plot - %%

%%% Floe %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% A MORE STABLE (BUT LONG) WAY OF CALC AMPS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% (IT USES PHI AND PHI_DR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
floe_displ = zeros(length(th_vec),length(rad_vec),Np);

v1 = 1:Vert_Dim*(2*Az_Dim+1);

for loop_Y=1:Y_Dim
 for loop_P=1:Np

  Amps = reshape(Psi(v1,loop_Y), Vert_Dim, 2*Az_Dim+1);
 
  Amps_dr = PsidrMat(v1,v1)\Psi_dr(v1,loop_Y); 
  Amps_dr = reshape(Amps_dr, Vert_Dim, 2*Az_Dim+1);

  tilAmps = zeros(Vert_Dim+2, 2*Az_Dim+1);
  tilAmps_dr = zeros(Vert_Dim+2, 2*Az_Dim+1);

  for loop_Az=1:Az_Dim
   tilAmps(:,Az_Dim+1-loop_Az) = [Amps(:,Az_Dim+1-loop_Az); ...
     Amp_mu_Mat(:,:,loop_Az+1)*Amps(:,Az_Dim+1-loop_Az)]; 
   tilAmps_dr(:,Az_Dim+1-loop_Az) = [Amps_dr(:,Az_Dim+1-loop_Az); ...
     Amp_mu_Mat(:,:,loop_Az+1)*Amps_dr(:,Az_Dim+1-loop_Az)]; 
  end
 
  for loop_Az=0:Az_Dim
   tilAmps(:,Az_Dim+1+loop_Az) = [Amps(:,Az_Dim+1+loop_Az); ...
     Amp_mu_Mat(:,:,loop_Az+1)*Amps(:,Az_Dim+1+loop_Az)];
   tilAmps_dr(:,Az_Dim+1+loop_Az) = [Amps_dr(:,Az_Dim+1+loop_Az); ...
     Amp_mu_Mat(:,:,loop_Az+1)*Amps_dr(:,Az_Dim+1+loop_Az)];
  end

  tilH = besselh(-Az_Dim:Az_Dim, [kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)]*R_vec(loop_P));
  tilHr = diag([kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)])*...
     Bessel_dz(@besselh, -Az_Dim:Az_Dim, [kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)]*R_vec(loop_P));
  tilJ = besselj(-Az_Dim:Az_Dim, [kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)]*R_vec(loop_P));
  tilJr = diag([kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)])*...
     Bessel_dz(@besselj, -Az_Dim:Az_Dim, [kk(:,loop_P);mu_0(loop_P);mu_1(loop_P)]*R_vec(loop_P));

  Amps = pi*R_vec(loop_P)*(tilHr.*tilJ.*tilAmps - tilH.*tilJr.*tilAmps_dr)/2i;
  Amps = Amps(1:Vert_Dim,:);
 
  Amps = reshape(Amps,1,Vert_Dim*(2*Az_Dim+1)); 
 
%  floe_displ(:,:,loop_P)=...
%      fn_PlotFloe(Az_Dim, Amps, th_vec, Disp_mat(:,:,loop_P));
 
  v1 = v1 + Vert_Dim*(2*Az_Dim+1);

 end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% - Free surface - %%%

xNm = [min(posits(1,:)-R_vec),max(posits(1,:)+R_vec)];

%displ_fs = zeros(size(FS_mesh));

%%% AMPLITUDES: 1st DIM = VERT DIM; 2nd DIM = Y DIM; 3rd DIM = INC WAVE
Bm = zeros(Vert_Dim,Y_Dim,2*Vert_Dim*Y_Dim); Ap = Bm; Am = Bm; Bp = Bm;

u1 = 1:Vert_Dim:Vert_Dim*Y_Dim;

for loop_Dim=1:Vert_Dim
 for loop_Y=1:Y_Dim
  Am(loop_Dim,loop_Y,u1(loop_Y)) = ...
      exp(1i*k0(loop_Dim)*v_vec(loop_Dim,loop_Y)*xNm(2));
  Bp(loop_Dim,loop_Y,Vert_Dim*Y_Dim+u1(loop_Y)) = ...
      exp(1i*k0(loop_Dim)*v_vec(loop_Dim,loop_Y)*lth)*...
      exp(-1i*k0(loop_Dim)*v_vec(loop_Dim,loop_Y)*xNm(1));
 end
 u1=u1+1;
end

Am = Am*Amp0/1i/k0(1); Bp = Bp*Amp0/1i/k0(1);

for loop_Y=1:2*Y_Dim*Vert_Dim
  
 v1 = 1:Vert_Dim*(2*Az_Dim+1);
 for loop_P=1:Np
    
  [dum_Bm, dum_Ap] = ...
    find_Amps(Vert_Dim, Az_Dim, Phi(v1,loop_Y), Phi_dr(v1,loop_Y),  ...
    k0, Y_Dim, [R_vec(loop_P),width], posits(:,loop_P), ...
    diag(k0)*u_vec, diag(k0)*v_vec, theta_inc, xNm, ExtraIn);

%  displ_fs = displ_fs + fn_PlotFS(parameter_vector, Vert_Dim, Az_Dim, Y_Dim, k0,... 
%     wt_0, width, R_vec(loop_P), dum_Am, dum_Bm, dum_Ap,...
%     posits(:,loop_P), u_vec, v_vec, x0, x_vec, y_vec, FS_mesh, ...
%     Phi(v1), Phi_dr(v1), res, ExtraIn);

  Bm(:,:,loop_Y) = Bm(:,:,loop_Y) + dum_Bm; 
  Ap(:,:,loop_Y) = Ap(:,:,loop_Y) + dum_Ap;
 
  v1 = v1 + Vert_Dim*(2*Az_Dim+1);
 end
end

%% - ENERGY CHECK - %%

real_inds = find(real(v_vec(1,:))~=0); 

epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;

for loop_Y=1:Vert_Dim:Y0_Dim*Vert_Dim

 EngErr = sum(epsm.*v_vec(1,real_inds).*...
     (abs(Am(1,real_inds,loop_Y)).^2 + abs(Bp(1,real_inds,loop_Y)).^2 ...
    - abs(Am(1,real_inds,loop_Y)+Ap(1,real_inds,loop_Y)).^2 ...
    - abs(Bp(1,real_inds,loop_Y)+Bm(1,real_inds,loop_Y)).^2));

 if abs(EngErr)>1e-3
  disp(['Energy error = ',num2str(abs(EngErr)),'!!!!!!!!!!!!!!!!'])
  disp(['Mode = ' int2str(loop_Y)])
 end   
 
 EngErr = sum(epsm.*v_vec(1,real_inds).*...
     (abs(Am(1,real_inds,Vert_Dim*Y_Dim+loop_Y)).^2 + abs(Bp(1,real_inds,Vert_Dim*Y_Dim+loop_Y)).^2 ...
    - abs(Am(1,real_inds,Vert_Dim*Y_Dim+loop_Y)+Ap(1,real_inds,Vert_Dim*Y_Dim+loop_Y)).^2 ...
    - abs(Bp(1,real_inds,Vert_Dim*Y_Dim+loop_Y)+Bm(1,real_inds,Vert_Dim*Y_Dim+loop_Y)).^2));

 if abs(EngErr)>1e-3
  disp(['Energy error = ',num2str(abs(EngErr)),'!!!!!!!!!!!!!!!!'])
  disp(['Mode = ' int2str(Y_Dim+loop_Y)])
 end   
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - For output - %%

%%% TAKE OUT AMPLITUDE %%%

Am = (1i*k0(1)/Amp0)*Am; Bm = (1i*k0(1)/Amp0)*Bm;
Ap = (1i*k0(1)/Amp0)*Ap; Bp = (1i*k0(1)/Amp0)*Bp;

Rm = zeros(Vert_Dim*Y_Dim); Rp=Rm; Tm=Rm; Tp=Rm;

v1 = 1:Vert_Dim;
for loop=1:Y_Dim
 Rm(v1,:) = squeeze(Bp(:,loop,1:Vert_Dim*Y_Dim)+Bm(:,loop,1:Vert_Dim*Y_Dim)); 
 Tp(v1,:) = squeeze(Bp(:,loop,Vert_Dim*Y_Dim+1:2*Vert_Dim*Y_Dim)+Bm(:,loop,Vert_Dim*Y_Dim+1:2*Vert_Dim*Y_Dim)); 
 Tm(v1,:) = squeeze(Am(:,loop,1:Vert_Dim*Y_Dim)+Ap(:,loop,1:Vert_Dim*Y_Dim)); 
 Rp(v1,:) = squeeze(Am(:,loop,Vert_Dim*Y_Dim+1:2*Vert_Dim*Y_Dim)+Ap(:,loop,Vert_Dim*Y_Dim+1:2*Vert_Dim*Y_Dim)); 
 v1 = v1+Vert_Dim;
end

if scat==1
 
 if 0
  disp('NEEDS REWRITING !!!!!!')   
  displ_fs = displ_fs + fn_PlotIncWv(parameter_vector, k0, wt_0, ...
     Amp0/1i/k0(1), x_vec, y_vec, FS_mesh);
 end
 
 Ap = Ap + Am; Bm = Bm + Bp;
end

u1 = 1:Vert_Dim:2*Vert_Dim*Y_Dim;

%Am = abs(Am(1,real_inds,u1)); Bp = abs(Bp(1,real_inds,u1));
Ap = abs(Ap(1,real_inds,u1)); Bm = abs(Bm(1,real_inds,u1)); 

En = zeros(1,2*length(real_inds));

for loop_Y=1:length(real_inds)
 En(loop_Y) = sum(epsm.*v_vec(1,real_inds).*(squeeze(Ap(1,:,loop_Y)).^2)/k0(1));
 %En = En + sum(epsm.*v_vec(1,real_inds).*(Bm.^2)/k0(1));
 En(Y0_Dim+loop_Y) = sum(epsm.*v_vec(1,real_inds).*(squeeze(Bm(1,:,Y_Dim+loop_Y)).^2)/k0(1));
end

display(['Transmitted energy    : ' num2str(abs(En))])

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Subfns - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% - THE INTEGRAL EQN STUFF - %%

function [PhiInc, G, Gdr_add, al_vec, be_vec, theta_inc, ExtraOut] = ...
    fn_IntEqn_Reson(Vert_Dim, Az_Dim, k0, Rad_vec, width, lth, ...
    posits, res, Amp0, evs, extra_pts)

%%% INPUTS
%
% evs = no. of evanescent waves included in inc wave (evs=0,1,2,...)

%%% FIND alpha_m beta_m & RESONANCES %%%

tol = 7.5e-4; %5e-2; % - tol on resonances

skip=[]; skipinfo=[];
count=1; al_vec(:,1)=ones(Vert_Dim,1); be_vec(:,1)=zeros(Vert_Dim, 1); Y_Dim=0; 
while be_vec(1,count)-1<tol %imag(vm)<Tols(5)
 be_vec(:,count+1) = (count*pi/width./k0).';
 al_vec(:,count+1) = sqrt(1-be_vec(:,count+1).^2);
    
 if and(Y_Dim==0,imag(al_vec(1,count+1))~=0)
  Y_Dim=count;
 end
    
 if abs(al_vec(1,count+1)-1)<tol %abs(vm(count))<Tols(5)
  skipinfo(1) = count;
  skipinfo(2) = be_vec(1,count+1);
  skipinfo(3) = al_vec(1,count+1);
  skip=count;
 end
 count=count+1;
end

al_vec(:,end) = []; be_vec(:,end) = [];
count = count-1; Y_Dim = count;

for loop_ev=1:evs
 be_vec(:,count+1) = count*pi/width./k0;
 al_vec(:,count+1) = sqrt(1-be_vec(:,count+1).^2);
 count=count+1;
end

clear count tol

theta_inc = zeros(size(al_vec)); % Nb 1st entry should be 0

for loop=1:length(al_vec(1,:))
 for loop_Dim=1:Vert_Dim   
 theta_inc(loop_Dim,loop) = GrafThetaVal(be_vec(loop_Dim,loop), ...
     al_vec(loop_Dim,loop), [1, 1]);
 end
end
  
%%% GREENS FNS & INC WAVE %%%

Np = length(Rad_vec);

G = zeros(Np*Vert_Dim*(2*Az_Dim+1)+Np*extra_pts, Np*Vert_Dim*(2*Az_Dim+1));
PhiInc = zeros(Np*Vert_Dim*(2*Az_Dim+1)+Np*extra_pts,2*Vert_Dim*(Y_Dim+evs));
Gdr_add = zeros(Np*extra_pts, Np*Vert_Dim*(2*Az_Dim+1));

v1=1:Vert_Dim*(2*Az_Dim+1)+extra_pts; 
v1i=[1:Vert_Dim:Vert_Dim*(2*Az_Dim+1),Vert_Dim*(2*Az_Dim+1)+ [1:extra_pts]];  
w1=1:extra_pts; 

for loop_p1=1:Np
 % - ORIGIN OF FIELD PT:   
 cc = posits(:,loop_p1); Rad=Rad_vec(loop_p1);  
 irreg_pts = fn_GenRandPts(cc, Rad, extra_pts);
 
 v2=1:Vert_Dim*(2*Az_Dim+1);     
 for loop_p2=1:Np
  % - ORIGIN OF SOURCE PT:   
  c0 = posits(:,loop_p2); Rad0=Rad_vec(loop_p2);
  if loop_p1==loop_p2
   [G(v1,v2), Gdr_add(w1,v2), ExtraOut] = fn_getG(...
       Vert_Dim, Az_Dim, k0, Rad, width, c0, res, irreg_pts);
  else
   [G(v1,v2), Gdr_add(w1,v2), ExtraOut] = fn_getGmix(...
       Vert_Dim, Az_Dim, k0, Rad, Rad0, width, cc, c0, res, irreg_pts);
  end
  v2=v2+Vert_Dim*(2*Az_Dim+1);
 end
 PhiInc(v1,:) = get_IncWv(Vert_Dim, Az_Dim, Y_Dim+evs, Amp0, k0, ...
     posits(:,loop_p1), al_vec, be_vec, theta_inc, lth, irreg_pts);  
 
 v1 =v1 +Vert_Dim*(2*Az_Dim+1)+extra_pts; 
 v1i=v1i+Vert_Dim*(2*Az_Dim+1)+extra_pts;
 w1 =w1 +extra_pts;
 
end

% if length(ExtraOut)>2
%  if isempty(ExtraOut{3})==3
%   disp('resonant case !!!!!!!!!')
%  end
%  if length(ExtraOut)>3
%   disp('irreg freq !!!!!!!!!')
%  end
% end

return

%

function PhiInc = get_IncWv(Vert_Dim, Az_Dim, Y_Dim, Amp0, kk, c0, ... 
    al_vec, be_vec, theta_inc, lth, irreg_pts) 

%%% DESCRIPTION
%
% Inc Wave: PhiInc = (Amp0/1i/kk)*Sum_{m=0}^{Phat} exp(i*alpha_m*x)cos(\beta_m*y)
% Inner prods taken wrt -exp(i*m*theta)/2*pi 
% on the circular bdy c0 + R*(cos(theta),sin(theta))

%%% INPUTS
%
% 

%%% COMMENTS
%
% - 27.05.10 - after noting mistake in Green's fn (but no y-dep so no change here)
% - 27.04.10 - take out the bessel fn weighting (for irreg freqs)

PhiInc = zeros(2*Az_Dim+1,2*Vert_Dim*Y_Dim);

v1 = 1:2*Az_Dim+1; u1 = 1:Vert_Dim;

%%% POSITIVE x DIRECTION

Sc = exp(1i*kk.*transpose(al_vec(:,1))*c0(1));

for loop_Dim=1:Vert_Dim
 PhiInc(v1,u1(loop_Dim)) = (1i.^[-Az_Dim:Az_Dim]).';
end
PhiInc(v1,u1) = PhiInc(v1,u1)*diag(Sc);
 
for loop=2:Y_Dim

 u1 = u1 + Vert_Dim;
 
 Sc = exp(1i*kk.*transpose(al_vec(:,loop))*c0(1));

 for loop_Dim=1:Vert_Dim
  PhiInc(v1,u1(loop_Dim)) = (cos([-Az_Dim:Az_Dim]*theta_inc(loop_Dim,loop)...
     -kk(loop_Dim)*be_vec(loop_Dim,loop)*c0(2)).*(1i.^[-Az_Dim:Az_Dim])).';
 end
 
 PhiInc(v1,u1) = PhiInc(v1,u1)*diag(Sc);

end

%%% NEGATIVE x DIRECTION

u1 = u1 + Vert_Dim;

theta_inc = theta_inc + pi;

Sc = exp(-1i*kk.*transpose(al_vec(:,1))*c0(1));
Sc = exp(1i*kk.*transpose(al_vec(:,1))*lth).*Sc;
for loop_Dim=1:Vert_Dim
 PhiInc(v1,u1(loop_Dim)) = ((-1i).^[-Az_Dim:Az_Dim]).';
end
PhiInc(v1,u1) = PhiInc(v1,u1)*diag(Sc);

for loop=2:Y_Dim

 u1 = u1 + Vert_Dim;
 
 Sc = exp(-1i*kk.*transpose(al_vec(:,loop))*c0(1));
 Sc = exp(1i*kk.*transpose(al_vec(:,loop))*lth).*Sc;
 
 for loop_Dim=1:Vert_Dim
  PhiInc(v1,u1(loop_Dim)) = (cos([-Az_Dim:Az_Dim]*theta_inc(loop_Dim,loop)...
     +kk(loop_Dim)*be_vec(loop_Dim,loop)*c0(2)).*(1i.^[-Az_Dim:Az_Dim])).';
 end
 PhiInc(v1,u1) = PhiInc(v1,u1)*diag(Sc);

end

% PhiInc_dr = kk*(1i.^[-Az_Dim:Az_Dim]).*...
%         Bessel_dz(@besselj,-Az_Dim:Az_Dim,kk*Rad);
%  
% PhiInc_dr = Sc*PhiInc_dr;
% 
% PhiInc_dr = reshape(PhiInc_dr, 2*Az_Dim+1,1);

%%% THE EXTRA PTS (FOR IRREGULAR FREQS) %%%

if isempty(irreg_pts)~=1

v1 = 2*Az_Dim+2;

for loop=1:length(irreg_pts(1,:))
 
 for loop_Y=1:length(al_vec)   
  PhiInc(v1,loop_Y) = exp(1i*al_vec(loop_Y)*irreg_pts(1,loop))*...
      cos(be_vec(loop_Y)*irreg_pts(2,loop));
 
  %PhiInc_dr(v1) = (Amp0)*exp(1i*kk(1)*irreg_pts(1,loop));
 end
 v1 = v1+1;
end

end

dum_PhiInc = (Amp0/1i/kk(1))*PhiInc;

PhiInc = zeros(Vert_Dim*(2*Az_Dim+1),2*Vert_Dim*Y_Dim);
u1i = 1:Vert_Dim:2*Vert_Dim*Y_Dim;
v1i = 1:Vert_Dim:Vert_Dim*(2*Az_Dim+1);
 
for loop=1:Vert_Dim
 
 PhiInc(v1i,u1i) = dum_PhiInc(:,u1i);
 
 u1i=u1i+1; v1i=v1i+1;
end
 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Solve the int eqns (noting possible resonance) - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi, Psi, Phi_dr, Psi_dr, ExtraOut] = fn_Solve(Vert_Dim, ...
    Az_Dim, k0, R_vec, width, PhiInc, G, Gdr_add, ...
    PsiMat, PsidrMat, AW0, AW, VT0, VT, extra_pts, ExtraIn)
% ------------------------------------------- %
Np=length(R_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(ExtraIn)==1 % - non-reson case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Gr = 0*G;

    BesJ = zeros(1,Vert_Dim*(2*Az_Dim+1)); BesJr=BesJ;
    
    w1 = 1:Vert_Dim*(2*Az_Dim+1);  
    w2 = 1:Vert_Dim*(2*Az_Dim+1);
   
    for loop_P=1:Np
     fac = ones(1,Vert_Dim*(2*Az_Dim+1))*(R_vec(loop_P));
     
     u1 = 1:Vert_Dim:Vert_Dim*(2*Az_Dim+1);
     for loop_Dim=1:Vert_Dim
      BesJ(u1) = besselj(-Az_Dim:Az_Dim, k0(loop_Dim)*R_vec(loop_P));
      BesJr(u1) = k0(loop_Dim)*Bessel_dz(@besselj, -Az_Dim:Az_Dim, k0(loop_Dim)*R_vec(loop_P));
      u1 = u1+1;
     end
     
     PhiInc(w1,:) = diag(BesJ)*PhiInc(w1,:);
     
     G(:,w2) = G(:,w2)*diag(fac);
     
     v1 = 1:Vert_Dim*(2*Az_Dim+1); 
     v2 = Vert_Dim*(2*Az_Dim+1)+[1:extra_pts];
     v3 = 1:extra_pts;
     for loop_Pr=1:Np
      Gr(v1,w2) = G(v1,w2)*diag(BesJr./BesJ);
     
      if extra_pts      
       Gr(v2,w2) = Gdr_add(v3,w2)*diag(fac);
      end
      v1 = v1 + Vert_Dim*(2*Az_Dim+1)+extra_pts;
      v2 = v2 + Vert_Dim*(2*Az_Dim+1)+extra_pts;
      v3 = v3 + extra_pts;
     end
      
     w1=w1+Vert_Dim*(2*Az_Dim+1)+extra_pts;
     w2=w2+Vert_Dim*(2*Az_Dim+1);
     
    end
    
    clear fac BesJ BesJr BesH BesHr u1 v1 v2 v3 w1 w2
    
    Ide = eye(Np*Vert_Dim*(2*Az_Dim+1)+Np*extra_pts);
    
    Idtil = zeros(Np*Vert_Dim*(2*Az_Dim+1)+Np*extra_pts,Np*Vert_Dim*(2*Az_Dim+1));
    
    v1 = 1:Vert_Dim*(2*Az_Dim+1)+extra_pts; v2 = 1:Vert_Dim*(2*Az_Dim+1);
    
    for loop=1:Np
    
     Idtil(v1,v2) = [eye((2*Az_Dim+1)*Vert_Dim);...
          zeros(extra_pts,(2*Az_Dim+1)*Vert_Dim)];
      
     v1 = v1 + Vert_Dim*(2*Az_Dim+1)+extra_pts;
     v2 = v2 + Vert_Dim*(2*Az_Dim+1);
     
    end
    
    %% - solve for auxiliary u (in terms of inc wave) - %%
    
    mat_Phi_Inc = ( Idtil + Gr )\Ide;
    
    mat_Phi_u = mat_Phi_Inc*(G*AW0);
    
    clear Jmat fac
    
    mat_Psi_u = (PsiMat/PsidrMat)*AW;
    
    %% - solve for auxiliary u (in terms of inc wave) - %%
    
    if max(max([AW0,AW]))==0
     u_aux = 0*PhiInc;
    else
     u_aux = (VT*mat_Psi_u-VT0*mat_Phi_u)\...
        VT0*mat_Phi_Inc*PhiInc;
    end
    
    %% - then solve for Phi_dr and Psi_dr - %%
    
    Phi_dr = AW0*u_aux; Psi_dr = AW*u_aux;
    
    %% - finally solve for Phi and Psi - %%
    
    Phi = mat_Phi_u*u_aux + mat_Phi_Inc*PhiInc;
    Psi = mat_Psi_u*u_aux;
    
    ExtraOut = [];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(ExtraIn)==1 % - so only reson freqs (ie vm=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% - subtract vp*const=(1/2*pi)*int(Gp'phi-Gp*phi')dtheta
% where Gp = exp(i*up*y)+exp(-i*up*y) (up=k for resonance)

    disp('needs to be rewritten')
    
    return

%     skipinfo=ExtraIn{1}{1};
%     
%     BesJ = besselj(-Az_Dim:Az_Dim, k0*Rad, 1);
%    
%     Gtil = diag(1./BesJ)*G;
%     Grtil = diag(1./BesJ)*G_dr*diag(BesJ);
%     
%     Id = eye((2*Az_Dim+1)*Vert_Dim);
%     
%     mat_Phi_u = (0.5*Id + (Rad/2/pi)*Grtil )\((Rad/2/pi)*Gtil*AW0);
%     
%     mat_Phi_Inc = Id/( 0.5*Id + (Rad/2/pi)*Grtil );
%     
%     % ----------------------------------------------------- %
%     % - The resonant part - %
%     % - Phi in terms of the constant in the expansion - %
%     %   const is lim = (1/2\pi*vp)*int{Gp'*phi-Gp*phi'}dtheta
%     %   -> Gp = exp{i*up*y}+exp{-i*up*y}
% 
%     Ipm = 1+(-1).^[-Az_Dim:Az_Dim];
%     Gp_mat(:,1) = Ipm.*besselj(-Az_Dim:Az_Dim,k0*Rad);
%     fac = Rad*(pi^2)/1i/width;
%     mat_Phi_const = -(0.5*Id + (Rad/2/pi)*Grtil )\(fac*Gp_mat);
%     
%     % ----------------------------------------------------- %
%     % - from the inner plate - %
%     
%     mat_Psi_u = (Id/PsidrMat)*AW;
%     
%     % ----------------------------------------------------- %
%     
% %     Gp_mat(1,:) = Ipm.*besselj(-Az_Dim:Az_Dim,k0*Rad);
%     Gp_mat = Gp_mat.';
%     Gpr_mat(1,:) = k0*Ipm.*Bessel_dz(@besselj,-Az_Dim:Az_Dim,k0*Rad);
%     
%     mat_const_const = Gpr_mat*diag(BesJ)*mat_Phi_const - skipinfo(3);
%     mat_const_Inc = -mat_const_const\Gpr_mat*diag(BesJ)*mat_Phi_Inc;
%     mat_const_u = -mat_const_const\(Gpr_mat*diag(BesJ)*mat_Phi_u ...
%                                    - Gp_mat*diag(BesJ)*AW0);
%     clear mat_const_const
%     
%     mat_Phi_Inc = mat_Phi_Inc + mat_Phi_const*mat_const_Inc;
%     mat_Phi_u = mat_Phi_u + mat_Phi_const*mat_const_u;
%     
%      %% - solve for auxiliary u (in terms of inc wave) - %%
%     
%     u_aux = (VT*PsiMat*mat_Psi_u-VT0*diag(BesJ)*mat_Phi_u)\...
%         VT0*diag(BesJ)*mat_Phi_Inc*PhiInc;
%     
%     %% - find the const limit - %%
%     
%     const = mat_const_Inc*PhiInc + mat_const_u*u_aux;
%     
%     %% - then solve for Phi_dr and Psi_dr - %%
%     
%     Phi_dr = AW0*u_aux; Psi_dr = AW*u_aux;
%     
%     %% - finally solve for Phi and Psi - %%
%     
%     Phi = mat_Phi_u*u_aux + mat_Phi_Inc*PhiInc;
%     Psi = mat_Psi_u*u_aux;
%     
%     % - Test - %
%     %(Gpr_mat*diag(BesJ)*Phi-Gp_mat*diag(BesJ)*Phi_dr)/2/pi;
%     % --------------------------------------------- %
%     
%     ExtraOut = const;
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elseif length(ExtraIn)==2 % - so only irreg freqs (ie. Jm=0)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     disp('needs to be rewritten')
%     
%     return
%     
%     irregs=ExtraIn{2}{1};
%     irregs_long = unique([Az_Dim+1+irregs,Az_Dim+1-irregs]);
%     
%     BesJ = besselj(-Az_Dim:Az_Dim, k0*Rad, 1);
%     
%     % - 'blank out' 'irregular' rows  - %
%     invBesJ = 1./BesJ; invBesJ(irregs_long)=0; 
%     % --------------------------------- %
%     Gtil = diag(invBesJ)*G; 
%     Grtil = diag(invBesJ)*G_dr;
%     % -------------------------------------------------- %
%     
%     Id = eye((2*Az_Dim+1)*Vert_Dim);
%     
%     mat_Phi_u = (0.5*Id + (Rad/2/pi)*Grtil )\((Rad/2/pi)*Gtil*AW0);
%     
%     mat_Phi_Inc = Id/( 0.5*Id + (Rad/2/pi)*Grtil );
%     
%     % ----------------------------------------------------- %
%     Id_ir = Id(:,irregs_long);
%     mat_Phi_ir = -(Rad/2/pi)*mat_Phi_Inc*Id_ir; % - size (2*Az_Dim+1)x|irregs|
%     % ----------------------------------------------------- %
%    
%     mat_Psi_u = (Id/PsidrMat)*AW;
%     
%     %% - use in extra eqns int{Gp'*phi - Gp*phi'}exp(-im\theta)d\theta = c1*Jp + ...
%     
%     vecGp = G(irregs_long,:); 
%     vecGrp = G_dr(irregs_long,:);
%     vec_rhs = BesJ(irregs_long);
%     
%     mat_ir_ir = vecGrp*mat_Phi_ir - diag(vec_rhs);
%     mat_ir_u = -mat_ir_ir\(vecGrp*mat_Phi_u - vecGp*AW0);
%     mat_ir_Inc = -mat_ir_ir\vecGrp*mat_Phi_Inc;
%     
%     clear mat_ir_ir
%     
%     mat_Phi_Inc = mat_Phi_Inc + mat_Phi_ir*mat_ir_Inc;
%     mat_Phi_u = mat_Phi_u + mat_Phi_ir*mat_ir_u;
%     
%     %% - solve for auxiliary u (in terms of inc wave) - %%
%     
%     u_aux = (VT*PsiMat*mat_Psi_u-VT0*diag(BesJ)*mat_Phi_u)\...
%         VT0*diag(BesJ)*mat_Phi_Inc*PhiInc;
%     
%     %% - find the const limit - %%
%     
%     ir = mat_ir_Inc*PhiInc + mat_ir_u*u_aux;
%     
%     %% - then solve for Phi_dr and Psi_dr - %%
%     
%     Phi_dr = AW0*u_aux; Psi_dr = AW*u_aux;
%     
%     %% - finally solve for Phi and Psi - %%
%     
%     Phi = mat_Phi_u*u_aux + mat_Phi_Inc*PhiInc;
%     Psi = mat_Psi_u*u_aux;
%     
%     ExtraOut =[];
%     
% 
% else % - both reson freqs & irreg freqs
%     
%     disp('need to write this part!!!!')
    
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Calculate amps in the rectangular regions either side of a plate - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Bm, Ap] = ...
    find_Amps(Vert_Dim, Az_Dim, Phi, Phi_dn, k0, ...
    Y_Dim, vars, c0, u_vec, v_vec, theta_inc, xNm, ExtraIn)

% [Am, Bm, Ap] = find_Amps(Vert_Dim, Phi_I, Phi, Phi_dn, k0, v_vec, u_vec, vars, c0)
%
% - vars = [Rad, width]
% - Bm and Ap are vectors of size Vert_Dim x Y_Dim
% - ExtraIn = [skip, constant limit] 
%   where const lim = (1/2\pi*vp)*int{Gp'phi-Gpphi'}dtheta: Gp = exp{i*up*y}+exp{-i*up*y}
%
% - Phi = sum_{0}^{Y_Dim-1}Ap/m exp(+/-i*v_m x)cos(u_m y)

% - 28.05.10 - changed after noting mistake in Green's fn
% - 22.09.10 - no incident wave included

Rad = vars(1);
w = vars(2); % width

Phi = reshape(Phi, Vert_Dim, 2*Az_Dim+1);
Phi_dn = reshape(Phi_dn, Vert_Dim, 2*Az_Dim+1);

if isempty(ExtraIn)==1
 skip=0;
else
 skip=ExtraIn(1)+1;
end

% - scaling - %
%BesJ = besselj(-Az_Dim:Az_Dim, k0*Rad);
%BesJr = k0*Bessel_dz(@besselj,-Az_Dim:Az_Dim, k0*Rad);

%Phi = BesJ.*Phi; %Phi_dn = BesJ.*Phi_dn;
% ------------------------------------- %

% theta_pp = theta_vecs{1}; theta_mm = theta_vecs{2};
% theta_mp = theta_vecs{3}; theta_pm = theta_vecs{4}; clear theta_vecs

Bm = zeros(Vert_Dim, Y_Dim); Ap = Bm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% - for xi<x   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% - v2 -> scale v1 so that Bm normalised at c0(1)-R & Ap at c0(1)+R

x0 = c0(1)+[-Rad,Rad]; 

for loop_Dim = 1:Vert_Dim
     
 Az = 1;   
 for loop_Az=-Az_Dim:Az_Dim
     
  %%% Nb. theta_inc(1)=0
     
  Bm(loop_Dim,1) = Bm(loop_Dim,1) - exp(-1i*(pi/2-theta_inc(loop_Dim,1))*loop_Az)*(...
      k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
      besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
  
%   Bm(loop_Dim,1) = Bm(loop_Dim,1) - exp(-1i*theta_pp(loop_Dim,1)*loop_Az)*(...
%       k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
%       besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );

  Az = Az+1; 
  
 end
 
 %%% SCALE FIELD TO LOWER LIMIT x0(1) & SOURCE TO DISK CENTRE c0
 Bm(loop_Dim,1) = exp(1i*v_vec(loop_Dim,1)*Rad+1i*u_vec(loop_Dim,1)*c0(2))...
     *Bm(loop_Dim,1)/v_vec(loop_Dim,1);
    
 for loop_Y = 2:Y_Dim   
     
  if loop_Y~=skip
        
   Az = 1;   
   for loop_Az=-Az_Dim:Az_Dim
       
    Bm(loop_Dim,loop_Y) = Bm(loop_Dim,loop_Y) - ...
      exp(-1i*loop_Az*pi/2)*cos(u_vec(loop_Dim,loop_Y)*c0(2)...
          +loop_Az*theta_inc(loop_Dim,loop_Y))*(...
      k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
      besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
   
    if 0
    Bm(loop_Dim,loop_Y) = Bm(loop_Dim,loop_Y) - ( ...
       exp( 1i*u_vec(loop_Dim,loop_Y)*c0(2))*exp(-1i*theta_pp(loop_Dim,loop_Y)*loop_Az) ...
     + exp(-1i*u_vec(loop_Dim,loop_Y)*c0(2))*exp(-1i*theta_pm(loop_Dim,loop_Y)*loop_Az) )*(...
      k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
      besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
    end 
    Az = Az+1; 
  
   end
  
   Bm(loop_Dim,loop_Y) = 2*exp(1i*v_vec(loop_Dim,loop_Y)*Rad)*...
       Bm(loop_Dim,loop_Y)/v_vec(loop_Dim,loop_Y);
  
  else
      
   Bm(loop_Dim,loop_Y) = -ExtraIn(2);
   
  end
  
 end
 
end

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - for xi>x   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for loop_Dim = 1:Vert_Dim
    
  Az = 1;   
  for loop_Az=-Az_Dim:Az_Dim
      
     Ap(loop_Dim,1) = Ap(loop_Dim,1) - exp(-1i*(3*pi/2-theta_inc(loop_Dim,1))*loop_Az)*(...
      k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
      besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
   
%    Ap(loop_Dim,1) = Ap(loop_Dim,1) - exp(-1i*theta_mp(loop_Dim,1)*loop_Az)*(...
%       k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
%       besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
  
   Az = Az+1; 
  
  end
  
  Ap(loop_Dim,1) = exp(1i*v_vec(loop_Dim,1)*Rad+1i*u_vec(loop_Dim,1)*c0(2))*...
      Ap(loop_Dim,1)/v_vec(loop_Dim,1);   
  
  for loop_Y = 2:Y_Dim  
     
   if loop_Y~=skip
     
    Az = 1;   
    for loop_Az=-Az_Dim:Az_Dim
        
       Ap(loop_Dim,loop_Y) = Ap(loop_Dim,loop_Y) - ...
        exp(-3i*pi*loop_Az/2)*cos(u_vec(loop_Dim,loop_Y)*c0(2)...
            - loop_Az*theta_inc(loop_Dim,loop_Y))*(...
        k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
        besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) ); 
     
      %%% TEST 
      if 0     
       Ap(loop_Dim,loop_Y) = Ap(loop_Dim,loop_Y) - ( ...
        exp( 1i*u_vec(loop_Dim,loop_Y)*c0(2))*exp(-1i*theta_mp(loop_Dim,loop_Y)*loop_Az) ...
        + exp(-1i*u_vec(loop_Dim,loop_Y)*c0(2))*exp(-1i*theta_mm(loop_Dim,loop_Y)*loop_Az) )*(...
        k0(loop_Dim)*Bessel_dz(@besselj, -loop_Az, Rad*k0(loop_Dim))*Phi(loop_Dim,Az) - ...
        besselj(-loop_Az, Rad*k0(loop_Dim))*Phi_dn(loop_Dim,Az) );
      end   
     Az = Az+1; 
  
    end
  
    Ap(loop_Dim,loop_Y) = 2*exp(1i*v_vec(loop_Dim,loop_Y)*Rad)*...
        Ap(loop_Dim,loop_Y)/v_vec(loop_Dim,loop_Y);
  
   else
      
    Ap(:,loop_Y) = -2*ExtraIn(2);
   
   end
  
 end
    
end

%Am = (pi*Rad/2i/w)*Am; 
Bm = (pi*Rad/1i/w)*Bm; Ap = (pi*Rad/1i/w)*Ap;

%% - Test Amps - %%

% vr = find(real(v_vec(1,:))~=0); %vi = find(real(v_vec)==0);
% 
% epsm = ones(1,length(vr)); epsm(1) = 2;
% 
% Am = [Amp0/1i/k0, zeros(1,length(vr)-1)];
% 
% ApI = [exp(1i*k0*x0(2))*Amp0/1i/k0, zeros(1,length(vr)-1)];
% 
% EngErr = sum(epsm.*v_vec(1,vr).*(abs(Am(1,vr)).^2 - abs(ApI+Ap(1,vr)).^2 - abs(Bm(1,vr)).^2));
% 
% if abs(EngErr)>1e-3
%  disp(['Energy error = ',num2str(abs(EngErr)),'!!!!!!!!!!!!!!!!'])
% end   

%%% NORMALIZE TO x=0

Bm = exp(1i*v_vec*(x0(1)-xNm(1))).*Bm;
Ap = exp(-1i*v_vec*(x0(2)-xNm(2))).*Ap;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ab = fn_AmpsBetween_vReson(Vert_Dim, Az_Dim, Y_Dim, x0, c0, res, vars,...
    u_vec, v_vec, Phi, Phi_dr, ExtraIn)

% - vars = [Rad,width];
%   ExtraIn = [skip, const]
%
% - function to calc the amplitudes around the floe pertaining to the part
% of the Green's fn that cannot have x-vars separated
% - if phi=Sum_{n=1,Y_Dim}phihat_n(x0)*cos(un*y0)
%   then phihat_n = phiInchat_n - phihat1_n - phihat2_n
% where phihat1_n = fac*exp(i*vn*x0)*int{G1n'phi-G1nphi'}dtheta
%       phihat2_n = fac*int{G2n'phi-G2nphi'}dtheta
% where fac = Rad/4i (n=1) and Rad/2i (n~=1)
%       G1n=G1n(x,y)=exp(i*vn*x)*cos(un*y) -> done analytically (elsewhere!)
%       G2n=G2n(x,y)=exp(i*vn*|x-x0|)*cos(un*y)  -> done numerically (here!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LIKE V2 BUT WITH PHI UNSCALED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% - 27.05.10 changed by LGB

if isempty(ExtraIn)==1
 skip=0;
else
 skip=ExtraIn(1)+1;
end

Ab = zeros(Vert_Dim,Y_Dim);
k0 = v_vec(1); fac = vars(1)/2i/vars(2);

Phi = reshape(Phi, Vert_Dim, 2*Az_Dim+1);
Phi_dr = reshape(Phi_dr, Vert_Dim, 2*Az_Dim+1);

% - scaling - %
BesJ = besselj(-Az_Dim:Az_Dim, k0*vars(1));
BesJr = k0*Bessel_dz(@besselj,-Az_Dim:Az_Dim, k0*vars(1));

%Phi = BesJ.*Phi; %Phi_dr = BesJ.*Phi_dr;
% ------------------------------------- %

% - find the angle of x0 (on the bdy) - %
% - nb. this will never be 0 or pi - %

if abs(x0-c0(1))<vars(1)
 ang0 = acos((x0-c0(1))/vars(1));

 theta_vec0=linspace(-ang0,ang0,res+1); % - x>x0
 theta_vec1=linspace(ang0,2*pi-ang0,res+1); % - x<x0
 
else
    
 theta_vec=linspace(-pi,pi,res+1); 
 
end
         
if exist('theta_vec0','var')==1
 for loop_Y=1:Y_Dim
     
  if loop_Y~=skip   
  
   Gqp_vals0 = zeros(1,length(theta_vec0)); Gqp_vals1 = Gqp_vals0; 

   for loop=1:length(theta_vec0)
    x = c0(1)+vars(1)*cos(theta_vec0(loop)); 
    y = c0(2)+vars(1)*sin(theta_vec0(loop));   
    Gqp_vals0(loop) = fn_Gqp0(x-x0,y,u_vec(loop_Y),v_vec(loop_Y));
    x = c0(1)+vars(1)*cos(theta_vec1(loop)); 
    y = c0(2)+vars(1)*sin(theta_vec1(loop));   
     Gqp_vals1(loop) = fn_Gqp0(x-x0,y,u_vec(loop_Y),v_vec(loop_Y));
   end
 
   for loop_Az=-Az_Dim:Az_Dim
     Evec0 = exp(1i*loop_Az*theta_vec0);
     Evec1 = exp(1i*loop_Az*theta_vec1);
     G = NumInt(Gqp_vals0.*Evec0,res,2*ang0) + ...
          NumInt(Gqp_vals1.*Evec1,res,2*(pi-ang0)); 
     Gr = (BesJr(Az_Dim+1+loop_Az)/BesJ(Az_Dim+1+loop_Az))*G;
     Ab(:,loop_Y) = Ab(:,loop_Y) - (Gr*Phi(:,Az_Dim+1+loop_Az) - ...
         G*Phi_dr(:,Az_Dim+1+loop_Az) );
   end
 
   clear Gqp_vals0 Gqp_vals1 Gqp_vals_dr0 Gqp_vals_dr1
 
   Ab(:,loop_Y) = Ab(:,loop_Y)/v_vec(loop_Y);
   
  else
      
   Ab(:,loop_Y) = -pi*ExtraIn(2);   
 
  end
 end
 
 Ab = fac*Ab; Ab(:,1) = Ab(:,1)/2;
  
else
    
 for loop_Y=1:Y_Dim
  
  if loop_Y~=skip   
     
   Gqp_vals = zeros(1,length(theta_vec)); 
  
   for loop=1:length(theta_vec)
    x = c0(1)+vars(1)*cos(theta_vec(loop)); 
    y = c0(2)+vars(1)*sin(theta_vec(loop));   
    Gqp_vals(loop) = fn_Gqp0(x-x0,y,u_vec(loop_Y),v_vec(loop_Y));
   end
 
   for loop_Az=-Az_Dim:Az_Dim
     Evec = exp(1i*loop_Az*theta_vec);
     G = NumInt(Gqp_vals.*Evec,res,2*pi); 
     Gr = (BesJr(Az_Dim+1+loop_Az)/BesJ(Az_Dim+1+loop_Az))*G;
     Ab(:,loop_Y) = Ab(:,loop_Y) - (Gr*Phi(:,Az_Dim+1+loop_Az) - ...
         G*Phi_dr(:,Az_Dim+1+loop_Az) );
   end
 
   clear Gqp_vals Gqp_vals_dr
 
   Ab(:,loop_Y) = Ab(:,loop_Y)/v_vec(loop_Y);
 
  else
      
   Ab(:,loop_Y) = -pi*ExtraIn(2);
   
  end
   
 end
 
 Ab = fac*Ab; Ab(:,1) = Ab(:,1)/2;
  
end
  
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = fn_Gqp0(xm,y,um,vm)

%xm = xx(1)-xx(2);

G = exp(1i*vm*abs(xm))*cos(um*y);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% - NumInt - %%

function I = NumInt(vec, res, l)

dt = l/res;

I = vec(1) + vec(res+1);

I = I + 2*sum(vec(2:res));

I = I*dt/2;

return 

%

function irreg_pts = fn_GenRandPts(cc, Rad, extra_pts)

if extra_pts
    
irreg_pts = zeros(2,extra_pts);

for loop=1:extra_pts
 rr = 0.8*Rad*rand;   
 th = 2*pi*rand;
 
 irreg_pts(1,loop) = cc(1) + rr*cos(th);
 irreg_pts(2,loop) = cc(2) + rr*sin(th);
 
end

else
    
irreg_pts = [];

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - THE FLOE STUFF - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Big_Psi, Big_Psi_dr, disp, Amp_mu_Mat] = ...
    fn_SingleFloe(pv, Vert_Dim, Az_Dim, kk, wt, mu_0, mu_1, Rad, r_vec)

%% - The ice-covered disc (uniform case) - %
%% -- Coeffs - %

C_mat = mat_C_dr_PWC(pv, Vert_Dim, kk, mu_0, mu_1, wt, pv(12));

Psi_0 = zeros(Vert_Dim, Vert_Dim, Az_Dim+1); 
Psi_dr_0 = zeros(Vert_Dim, Vert_Dim, Az_Dim+1);

Amp_mu_Mat = zeros(2, Vert_Dim, 2*Az_Dim+1);
%ilPsi_Mat = zeros(Vert_Dim+2, 2*Az_Dim+1); tilPsidr_Mat = tilPsi_Mat;

for loop_Az=0:Az_Dim % - calc for each azi mode (>=0)

    % - Bend Moment & Shear Stress - %

    [bend_mu0, shear_mu0] = ...
        BendShear3dAxi_fn( pv, [1, loop_Az, 1], mu_0*Rad );
    [bend_mu1, shear_mu1] = ...
        BendShear3dAxi_fn( pv, [1, loop_Az, 1], mu_1*Rad );
    [bend_K, shear_K] = ...
        BendShear3dAxi_fn( pv, [1, loop_Az, 1], diag([kk*Rad]) );

    % - Bessel fns - %

    J_R = diag([besselj(loop_Az, kk*Rad, 1 )]);
    Jmu_R = diag( besselj(loop_Az, [ mu_0; mu_1 ].*Rad, 1) );

    Jz_R = diag([Bessel_dz(@besselj, loop_Az, kk*Rad, 1 )]);
    Jzmu_R = diag( Bessel_dz(@besselj, loop_Az, [ mu_0; mu_1 ].*Rad, 1 ) );

    Amp_mu(1,1:Vert_Dim) = -C_mat(Vert_Dim+1, 1:Vert_Dim)*...
        (bend_mu1*shear_K - shear_mu1*bend_K);
    Amp_mu(2,1:Vert_Dim) = C_mat(Vert_Dim+1, 1:Vert_Dim)*...
        (bend_mu0*shear_K - shear_mu0*bend_K);

    Amp_mu = -Amp_mu./(bend_mu0*shear_mu1 - bend_mu1*shear_mu0);
    % Amp_mu = -Amp_mu./pv(6);

    Amp_mu_Mat(:,:,Az_Dim+loop_Az+1) = Amp_mu;
    if loop_Az~=0
     Amp_mu_Mat(:,:,Az_Dim-loop_Az+1) = Amp_mu;
    end
    
    % - fn psi = %

    Psi_0(:,:,1+loop_Az) = J_R + ...
        C_mat(1:Vert_Dim, Vert_Dim+1:Vert_Dim+2)*Jmu_R*Amp_mu;
    Psi_dr_0(:,:,1+loop_Az) = diag(kk)*Jz_R + ...
        C_mat(1:Vert_Dim, Vert_Dim+1:Vert_Dim+2)*diag([mu_0, mu_1])*Jzmu_R*Amp_mu;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amalgamate azi dims (using fact that they are decoupled & symmetry) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Big_Psi = zeros((2*Az_Dim+1)*Vert_Dim); Big_Psi_dr = Big_Psi;
v1 = 1:Vert_Dim;

% - the negative ones (m<0)

for loop_Az=1:Az_Dim
    Big_Psi(v1,v1) = ((-1)^(Az_Dim+1-loop_Az))*Psi_0(:,:,Az_Dim+2-loop_Az);
    Big_Psi_dr(v1,v1) = ((-1)^(Az_Dim+1-loop_Az))*Psi_dr_0(:,:,Az_Dim+2-loop_Az);
    
    v1=v1+Vert_Dim;
end

% - the rest (m>=0)

for loop_Az=0:Az_Dim
    Big_Psi(v1,v1) = Psi_0(:,:,1+loop_Az);
    Big_Psi_dr(v1,v1) = Psi_dr_0(:,:,1+loop_Az);
    
    v1=v1+Vert_Dim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% - Disp in floe - %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_r = length(r_vec);

xi_vec = C_mat(Vert_Dim+1, 1:Vert_Dim);
chi_vec = C_mat(Vert_Dim+1, Vert_Dim+1:Vert_Dim+2);

disp = zeros(Vert_Dim*(2*Az_Dim+1), res_r);

for loop_r=1:res_r   
 for loop_Az=-Az_Dim:Az_Dim
     
  Amp_mu = diag(exp(-imag([mu_0 mu_1])*Rad))*Amp_mu_Mat(:,:,Az_Dim+loop_Az+1);
    
  J_mat = diag(besselj(loop_Az, kk*r_vec(loop_r) ));
  Jmu = diag(besselj(loop_Az, [ mu_0; mu_1 ].*r_vec(loop_r)) );
   
  disp(Az_Dim+loop_Az+[1:Vert_Dim],loop_r) = ( xi_vec*J_mat + chi_vec*Jmu*Amp_mu );
  
 end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - THE JCs - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Big_AW0,Big_AW,Big_VT0,Big_VT] = ...
    fn_JumpMats(pv,Vert_Dim,DimG,Az_Dim,k0,kk,wt_0,wt)

%% -- coeffs of the JCs at r=R -- %

mat_A = matrix_A_PWC(pv, Vert_Dim, kk, wt, pv(10));
mat_A0 = matrix_A_PWC(pv, Vert_Dim, k0, wt_0, pv(7));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% - Old way (vertical modes) - %%%
%[mat_W0, mat_W] = jump_W_PWC(pv, Vert_Dim, DimG, kk, k0, kk, wt_0, wt);
[mat_W0, mat_W] = jump_W_PWC_wtd(pv, Vert_Dim, kk, k0, kk, wt, wt_0, wt);
%%% - New way (Gegenbauer) - %%%
% if pv(9)==0
%     [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 0.5, DimG);
% else
%     [mat_W0, mat_W] = Jump_Gegen(k0, kk, wt_0, wt, pv(7), pv(8), 1/6, DimG);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -- Put jumps together in big matrices - %

Big_AW0 = zeros((2*Az_Dim+1)*Vert_Dim, (2*Az_Dim+1)*DimG); Big_AW = Big_AW0;
Big_VT0 = Big_AW0.'; Big_VT = Big_VT0;

v1 = 1:Vert_Dim;
v2 = 1:DimG;

for loop_Az = 1:2*Az_Dim+1
    Big_AW0(v1, v2) = mat_A0\mat_W0;
    Big_AW(v1, v2) = mat_A\mat_W;
    
    Big_VT0(v2,v1) = mat_W0.';
    Big_VT(v2,v1) = mat_W.';
    
    v1 = v1+Vert_Dim; v2 = v2+DimG;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Plot the disp in floe - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displ = fn_PlotFloe(Az_Dim, Amps, th_vec, dum_disp)

res_th = length(th_vec); 

dum_disp = diag(Amps)*dum_disp;

displ = ones(res_th,1)*dum_disp(Az_Dim+1,:);

for loop_Az=1:Az_Dim;    
    displ = displ + exp(1i*loop_Az*th_vec)*dum_disp(Az_Dim+1+loop_Az,:)+...
        exp(-1i*loop_Az*th_vec)*dum_disp(Az_Dim+1-loop_Az,:);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Plot the surrounding wavefield - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displ = fn_PlotFS(pv, Vert_Dim, Az_Dim, Y_Dim, k0, wt_0, width, Rad,... 
    Am, Bm, Ap, c0, u_vec, v_vec, x0, x_vec, y_vec, FS_mesh, ...
    Phi_bdy, Phi_dr_bdy, res, ExtraIn)

x_res = length(x_vec); y_res = length(y_vec);

displ = zeros(y_res,x_res);

% - need to mesh the area surrounding the floe - %

% x_vec = linspace(0,c0(1)+2*Rad,x_res);
% y_vec = linspace(-width,width,y_res);

for loop_x=1:x_res
 xx = x_vec(loop_x);
 if xx<=c0(1)-Rad % - between wave maker & near edge of floe
  for loop_y=1:y_res
   if FS_mesh(loop_y,loop_x)==0   
   yy = y_vec(loop_y);
   Phi=diag(exp(1i*v_vec(:,1)*xx))*Am(:,1) ...
        + diag(exp(-1i*v_vec(:,1)*(xx-x0(1))))*Bm(:,1); 
    for loop_Y=2:Y_Dim 
     Phi=Phi + cos(u_vec(loop_Y)*yy)*( ...
         diag(exp(1i*v_vec(:,loop_Y)*xx))*Am(:,loop_Y) ...
         + diag(exp(-1i*v_vec(:,loop_Y)*(xx-x0(1))))*Bm(:,loop_Y)); 
    end
    displ(loop_y,loop_x) = sum(wt_0.*cosh(k0*pv(7)).*transpose(Phi));
    else
     displ(loop_y,loop_x) = nan+1i*nan;
    end
  end
 elseif xx>=c0(1)+Rad % - far edge of floe & far-field
  for loop_y=1:y_res
   yy = y_vec(loop_y);
    Phi=diag(exp(1i*v_vec(:,1)*(xx-x0(2))))*Ap(:,1); 
    for loop_Y=2:Y_Dim 
     Phi=Phi + cos(u_vec(loop_Y)*yy)*...
         diag(exp(1i*v_vec(:,loop_Y)*(xx-x0(2))))*Ap(:,loop_Y); 
    end
    displ(loop_y,loop_x) = sum(wt_0.*cosh(k0*pv(7)).*transpose(Phi));
  end
 else % - around floe
  % - version 2: Integral expression method - %
  Ab = fn_AmpsBetween_vReson(Vert_Dim, Az_Dim, Y_Dim, xx, c0, res,... 
      [Rad,width], u_vec, v_vec, Phi_bdy, Phi_dr_bdy, ExtraIn);
   for loop_y=1:y_res
    if FS_mesh(loop_y,loop_x)==0
     yy = y_vec(loop_y);
     Phi=diag(exp(1i*v_vec(:,1)*xx))*Am(:,1) + Ab(:,1); 
      for loop_Y=2:Y_Dim 
       Phi=Phi + cos(u_vec(loop_Y)*yy)*( ...
         diag(exp(1i*v_vec(:,loop_Y)*xx))*Am(:,loop_Y) + Ab(:,loop_Y)); 
      end
      displ(loop_y,loop_x) = sum(wt_0.*cosh(k0*pv(7)).*transpose(Phi)); 
    else
      displ(loop_y,loop_x) = nan+1i*nan;
    end
   end  
 end
end

% figure(fig)
% 
% surf(x_vec,y_vec,real(displ))
% 
% shading interp
% 
% xlabel('x'); ylabel('y')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Plot the incident wave in the free-surf region - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displ = fn_PlotIncWv(pv, k0, wt_0, IncAmp, x_vec, y_vec, FS_mesh)

displ = FS_mesh;

for loop_x=1:length(x_vec)
 xx = x_vec(loop_x);
 Phi = exp(1i*k0*xx)*IncAmp;
 
 for loop_y=1:length(y_vec)
  
  if FS_mesh(loop_y,loop_x)==0
   displ(loop_y,loop_x) = sum(wt_0*cosh(k0*pv(7))*Phi);
  end
  
 end
 
end

return
