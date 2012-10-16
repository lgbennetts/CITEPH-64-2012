function [G, Gdr_add, ExtraOut] = fn_getGmix(...
    Dim, Az_Dim, kk, Rad, Rad0, width, cc, c0, res, irreg_vals)

skip = []; % INVESTIGATE LATER (12.10.12)

% - 27.09.12 - modified from find_getGmix for irregular values
% - 27.05.10 - after noting mistake in Green's fn
% - 27.04.10 - take care of reson freqs (vm=0) and irreg freqs (Jm=0) with
%              perturbation methods
% - 26.04.10 - converge with respect to the numerical integration
% - nb calc amplitudes free of Bessel fns to avoid problems at zeros

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - nb - check that
%        (Rad/2/pi)*diag(1/BesJ)*(Gj_dr*diag(BesJ)*PhiInc-Gj*diag(BesJ)*Phi
%        Inc_dr)= { PhiInc/4 (j=p or m)
%               = { 0        (j=til)   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tols(1) = 100; % - max num of terms in Green's fns
Tols(2) = 1e-1; % - cut off tol for terms in Green's fns
Tols(3) = 1e-2; % - the point at which the sings are dealt with anay in Logs
Tols(4) = 1e-4; % - accuracy of numerical integration
Tols(5) = 7.5e-4; %5e-2; % - tol on resonances
Tols(6) = 1e-1; % - tol on irreg freqs
           
% INITIALISE

theta_vec=linspace(0,2*pi,res+1);
Gqp_vals = zeros(res); Gqpmod_vals = Gqp_vals; 
gsk_vals = zeros(res); gskmod_vals = gsk_vals; % - the skipped values!!!

% - x is field; x0 is source - %

xx=cc(1)+Rad*cos(theta_vec); yy=cc(2)+Rad*sin(theta_vec);
x0=c0(1)+Rad0*cos(theta_vec); y0=c0(2)+Rad0*sin(theta_vec);

% GREEN'S FUNCTION %

if isempty(skip)==1
 for loop_N=1:Dim   
  vars=[nan,width,kk(loop_N),Rad];
  for loop=1:length(theta_vec)
   for loop2=1:length(theta_vec)
    [Gqp_vals(loop2,loop,loop_N)] = fn_GreensFnQPmod...
       ([xx(loop),yy(loop)],[x0(loop2),y0(loop2)],vars,skip,Tols);
    [Gqpmod_vals(loop2,loop,loop_N)] = fn_GreensFnQPmod...
       ([xx(loop),yy(loop)],[x0(loop2),-y0(loop2)],vars,skip,Tols);
   end
  end
 end
else
 for loop=1:length(theta_vec)
  for loop2=1:length(theta_vec)
   [Gqp_vals(loop2,loop),gsk_vals(loop,loop2)] = fn_GreensFnQPmod...
       ([xx(loop),yy(loop)],[x0(loop2),y0(loop2)],vars,skip,Tols);
   [Gqpmod_vals(loop2,loop),gskmod_vals(loop,loop2)] = fn_GreensFnQPmod...
       ([xx(loop),yy(loop)],[x0(loop2),-y0(loop2)],vars,skip,Tols);
  end
 end
end

% - Find their inner-prods with Fourier exps - %

Gqp = fn_FillMat_NoSymm(Gqp_vals, Dim, Az_Dim, theta_vec);
Gqpmod = fn_FillMat_NoSymm(Gqpmod_vals, Dim, Az_Dim, theta_vec);

if isempty(skip)==0
 gv = fn_FillMat_SymmD(gsk_vals, Dim, Az_Dim, theta_vec);
 gvmod = fn_FillMat_NoSymm(gskmod_vals, Dim, Az_Dim, theta_vec);
end

G = (Gqp + Gqpmod)/2/pi;  

if isempty(skip)~=1
 Gsk = gsktil + gv + gvmod;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IRREGULAR VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ASSUMING THEY ARE ONLY CAUSED BY REAL WAVENUMBER

Gdr_add = [];

if ~isempty(irreg_vals)
    
 Gqp_vals = zeros(length(irreg_vals(1,:)),res); Gqpmod_vals = Gqp_vals;
 vars=[nan,width,kk(1),Rad];
  
 for loop=1:length(irreg_vals(1,:))
  for loop2=1:length(theta_vec)
   [Gqp_vals(loop,loop2)] = fn_GreensFnQPmod...
     (irreg_vals(:,loop).',[x0(loop2),y0(loop2)],vars,skip,Tols);
   [Gqpmod_vals(loop,loop2)] = fn_GreensFnQPmod...
     (irreg_vals(:,loop).',[x0(loop2),-y0(loop2)],vars,skip,Tols);
  end   
 end
 
 G_add = zeros(length(irreg_vals(1,:)),Dim*(2*Az_Dim+1)); Gdr_add = G_add;
 v1 = 1:Dim:Dim*(2*Az_Dim+1);
 
 G_add(:,v1) = fn_FillVec(Gqp_vals+Gqpmod_vals, 1, Az_Dim, theta_vec);
 G = [G;G_add];
 
 Gqpdr_vals = zeros(length(irreg_vals(1,:)),res); Gqpmoddr_vals = Gqpdr_vals;
 for loop=1:length(irreg_vals(1,:))
  for loop2=1:length(theta_vec)
   [Gqpdr_vals(loop,loop2)] = fn_GreensFndr_til...
     ([x0(loop2),y0(loop2)],irreg_vals(:,loop).',c0,vars,Tols);
   [Gqpmoddr_vals(loop,loop2)] = fn_GreensFndr_til...
     ([x0(loop2),y0(loop2)],[1,-1].*irreg_vals(:,loop).',c0,vars,Tols);
   end   
  end
 
 Gdr_add(:,v1) = fn_FillVec(Gqpdr_vals+Gqpmoddr_vals, 1, Az_Dim, theta_vec);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Extra info for resonances, etc
ExtraOut = {[], res};
if isempty(skip)~=1
 SkipOut = {skipinfo, Gsk, Gsk_dr};
 ExtraOut{3} = SkipOut;
end

return
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Subfunctions - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% - Matrix Fills - %%

function Mat = fn_FillMat_NoSymm(Vals, Dim, Az_Dim, theta_vec)

% - No symmetry of matrix

Mat = zeros(2*Az_Dim+1);

v1 = 1:Dim; loop_Az1 = -Az_Dim; 

for loop1=1:2*Az_Dim+1
 v2 = 1:Dim; loop_Az2 = -Az_Dim;
 for loop2=1:2*Az_Dim+1   
 Mat(v1, v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
 
 loop_Az2 = loop_Az2+1;
 v2 = v2 + Dim; 
 end
 v1 = v1+Dim; loop_Az1 = loop_Az1+1;
end 

return

function Mat = fn_FillVec(Vals, Dim, Az_Dim, theta_vec)

m = length(Vals(:,1,1));

Mat = zeros(Dim*m,Dim*(2*Az_Dim+1));

v1 = 1:m*Dim; 
v2 = 1:Dim; loop_Az2 = -Az_Dim;
for loop2=1:2*Az_Dim+1   
 Mat(v1, v2) = get_G_mat_1D(loop_Az2,Vals,theta_vec,Dim);
 
 loop_Az2 = loop_Az2+1;
 v2 = v2 + Dim; 
end

return

%

function G = get_G_mat_1D(az, G0_vals, theta_vec, Dim)

m = length(G0_vals(:,1,1));
G = zeros(m*Dim,Dim);
res=length(theta_vec);

Exp_vec(1,:)=exp(1i*az*theta_vec);

v1 = 1:Dim;
for loop_m=1:m
 for loop_Dim=1:Dim
  dum_Mat = G0_vals(loop_m,:,loop_Dim).*Exp_vec;
  G(v1(loop_Dim),loop_Dim) = NumInt_1D(dum_Mat, res);        
 end
 v1 = v1+Dim;
end

return

%

function G = get_G_mat(az1, az2, G0_vals, theta_vec, Dim)

G = zeros(Dim);
res=length(theta_vec);

Exp_vec1 = zeros(res,1); Exp_vec2 = zeros(1,res);

Exp_vec1(:,1)=exp(1i*az1*theta_vec);
Exp_vec2(1,:)=exp(1i*az2*theta_vec);

for loop_Dim=1:Dim
 dum_Mat = G0_vals(:,:,loop_Dim).*(Exp_vec1*Exp_vec2);
 G(loop_Dim,loop_Dim) = NumInt_2D(dum_Mat, res);        
end

return

% 1D NUMERICAL INTEGRATION OVER [-PI,PI] 

function I = NumInt_1D(M, res)

dt = 2*pi/(res-1);

I = M(1) +  M(res);

I = I + 2*sum(M(2:res-1));

I = I*dt/2;

return 

% 2D NUMERICAL INTEGRATION OVER [-PI,PI]^2 

function I = NumInt_2D(M, res)

dt = 2*pi/(res-1);

I = M(1,1) + M(res,1) + M(1,res) + M(res,res);

I = I + 2*(sum(M(1,2:res-1)) + sum(M(2:res-1,1)) + sum(M(res,2:res-1)) + sum(M(2:res-1,res)));

I = I + 4*sum(sum(M(2:res-1,2:res-1)));

I = I*(dt^2)/4;

return 
