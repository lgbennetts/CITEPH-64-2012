function [G, Gdr_add, ExtraOut] = fn_getG(...
    Dim, Az_Dim, kk, Rad, width, c0, res, Tols, irreg_vals, skip, skipinfo)

% - 26.09.12 - modified from find_getG for irregular values
% - 21.09.12 - modified from find_GandGr_vTol_Resonx2_v3
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

% Tols(1) = 100; % - max num of terms in Green's fns
% Tols(2) = 1e-1; % - cut off tol for terms in Green's fns
Tols(3) = 1e-2; % - the point at which the sings are dealt with anay in Logs
Tols(4) = 1e-4; % - accuracy of numerical integration
Tols(5) = 7.5e-4; %5e-2; % - tol on resonances
Tols(6) = 1e-1; % - tol on irreg freqs
           
% - INITIALISE

theta_vec=linspace(-pi,pi,res+1);

xx=c0(1)+Rad*cos(theta_vec); yy=c0(2)+Rad*sin(theta_vec);

% - GREEN'S FN - %

% NB. source vals vary in dim 1 (vertically); field dim 2 (horizontally)
% BUT doesn't make a difference here (does in find_getGmix)

Gqp_vals = zeros(res,res,Dim); Gqpmod_vals = zeros(res,res,Dim); %Gtil_vals = Gqp_vals; S0_vals = zeros(res);
%gsk_vals = zeros(res,res,Dim); % - the skipped values!!!

for loop_N=1:Dim   
 vars=[nan,width,kk(loop_N),Rad];
 for loop=1:length(theta_vec)
  for loop2=1:length(theta_vec)
   [Gqp_vals(loop,loop2,loop_N)] = ...
       fn_GreensFnQuasiPer([xx(loop),yy(loop)],[xx(loop2),yy(loop2)],c0,vars,skip,Tols);
   [Gqpmod_vals(loop,loop2,loop_N)] = ...
       fn_GreensFnQPmod([xx(loop),yy(loop)],[xx(loop2),-yy(loop2)],vars,skip,Tols);
  end
 end
end
% else
%  for loop=1:length(theta_vec)
%   for loop2=1:length(theta_vec)
%    [Gqp_vals(loop,loop2),gsk_vals(loop,loop2)] = fn_GreensFnQuasiPer_Reson_v2...
%        ([xx(loop),yy(loop)],[xx(loop2),yy(loop2)],c0,vars,skip,Tols);
%    [Gqpmod_vals(loop,loop2),gskmod_vals(loop,loop2)] = fn_GreensFnQPmod_Reson...
%        ([xx(loop),yy(loop)],[xx(loop2),-yy(loop2)],vars,skip,Tols);
%   end
%  end
% end

% - Find their inner-prods with Fourier exps - %

Gqp = fn_FillMat_SymmD(Gqp_vals, Dim, Az_Dim, theta_vec);
Gqpmod = fn_FillMat_SymmE(Gqpmod_vals, Dim, Az_Dim, theta_vec);

clear Gqp_vals Gqpmod_vals

% if isempty(skip)==0
%  gv = fn_FillMat_SymmD(gsk_vals, Dim, Az_Dim, theta_vec);
%  gvmod = fn_FillMat_NoSymm(gskmod_vals, Dim, Az_Dim, theta_vec);
%  clear gsk_vals gskmod_vals
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    
Gqp_vals0 = zeros(res+1); Gqpmod_vals0 = Gqp_vals0; Ln_vals = Gqp_vals0;

for loop=1:length(theta_vec)    
 for loop2=1:length(theta_vec)
 [Gqp_vals0(loop,loop2)] = fn_GreensFnQuasiPer_LessBasic...
  ([xx(loop),yy(loop)],[xx(loop2),yy(loop2)],vars,skip,[Tols(1:2),1e-8,Tols(4:end)]);
 [Gqpmod_vals0(loop,loop2)] = fn_GreensFnQuasiPer_LessBasic...
  ([xx(loop),yy(loop)],[xx(loop2),-yy(loop2)],vars,skip,[Tols(1:2),1e-8,Tols(4:end)]);
 end
end

Gqp0 = fn_FillMat_NoSymm(Gqp_vals0+Gqpmod_vals0, Dim, Az_Dim, theta_vec);

for loop=1:length(theta_vec)
 for loop2=1:length(theta_vec)
 Ln_vals(loop,loop2) = fn_LogSing...
  ([xx(loop),yy(loop)],[xx(loop2),yy(loop2)],c0,-1,vars,[Tols(1:2),1e-8,Tols(4:end)]);
 end
end
Ln = fn_FillMat_NoSymm(Ln_vals, Dim, Az_Dim, theta_vec); 

Gqp1 = Gqp + Gqpmod + Ln;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Log_Mat = zeros(1,Dim*(2*Az_Dim+1));

v1 = 1:Dim; loop_Az=-Az_Dim; v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;
for loop=1:Az_Dim
 Log_Mat(1,v1) = 2*fn_LogSingmod(theta_vec,loop_Az);    
 Log_Mat(1,v1_flip) = Log_Mat(1,v1); 
 v1 = v1+Dim; loop_Az=loop_Az+1; v1_flip = v1_flip-Dim;
end
%Log_Mat(1,v1) = fn_LogSingmod(theta_vec,loop_Az);
Log_Mat(1,v1) = 4*pi*log(sqrt(2)*Rad*pi/width);

Log_Mat = Log_Mat - 2*pi*log(2);

Log_Mat = Log_Mat/2;

Gqp = Gqp+diag(Log_Mat);

%% -

%G = Gtil + Gqp + Gqpmod;  
G = (Gqp + Gqpmod)/2/pi;  

% if isempty(skip)~=1
%  Gsk = gv + gvmod; % (10.01.13) What was gsktil +  ????
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% IRREGULAR VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% ASSUMING THEY ARE ONLY CAUSED BY REAL WAVENUMBER

Gdr_add = [];

if ~isempty(irreg_vals)
    
 Gqp_vals = zeros(length(irreg_vals(1,:)),res); Gqpmod_vals = Gqp_vals;
 vars=[nan,width,kk(1),Rad];
  
 for loop=1:length(irreg_vals(1,:))
  for loop2=1:length(theta_vec)
   [Gqp_vals(loop,loop2)] = fn_GreensFnQPmod...
     (irreg_vals(:,loop).',[xx(loop2),yy(loop2)],vars,skip,Tols);
   [Gqpmod_vals(loop,loop2)] = fn_GreensFnQPmod...
     (irreg_vals(:,loop).',[xx(loop2),-yy(loop2)],vars,skip,Tols);
  end   
 end
 
 G_add = zeros(length(irreg_vals(1,:)),Dim*(2*Az_Dim+1)); Gdr_add = G_add;
 v1 = 1:Dim:Dim*(2*Az_Dim+1);
 
 G_add(:,v1) = fn_FillVec(Gqp_vals+Gqpmod_vals, 1, Az_Dim, theta_vec);
 G = [G;G_add];
 
 Gqpdr_vals = zeros(length(irreg_vals(1,:)),res); Gqpmoddr_vals = Gqpdr_vals;
 for loop=1:length(irreg_vals(1,:))
  for loop2=1:length(theta_vec)
%    [Gqpdr_vals0(loop,loop2)] = fn_Gpm_dr_v0...
%      (irreg_vals(:,loop).',[xx(loop2),yy(loop2)],c0,vars,Tols);
%    [Gqpmoddr_vals0(loop,loop2)] = fn_Gpm_dr_v0...
%      (irreg_vals(:,loop).',[xx(loop2),-yy(loop2)],c0,vars,Tols);   
   [Gqpdr_vals(loop,loop2)] = fn_GreensFndr_til...
     ([xx(loop2),yy(loop2)],irreg_vals(:,loop).',c0,vars,Tols);
   [Gqpmoddr_vals(loop,loop2)] = fn_GreensFndr_til...
     ([xx(loop2),yy(loop2)],[1,-1].*irreg_vals(:,loop).',c0,vars,Tols);
   end   
  end
 
 Gdr_add(:,v1) = fn_FillVec(Gqpdr_vals+Gqpmoddr_vals, 1, Az_Dim, theta_vec);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Extra info for resonances, etc
ExtraOut = []; %{[], res};
if isempty(skip)~=1
%  SkipOut = {skipinfo}; %, Gsk_dr (30.07.12) , Gsk (10.01.13)
%  ExtraOut{3} = SkipOut;
 ExtraOut = {skipinfo};
end

% if isempty(irregs)~=1
%  IrregOut = {irregs};
%  ExtraOut{4} = IrregOut;
% end

return
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - Subfunctions - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% - Matrix Fills - %%

% % - this uses analytic expressions
% 
% function Mat = fn_FillMat_SymmA_anlyt(funct, num_mats, Dim, Az_Dim, fn_vars)
% 
% % - nb funct = funct(Az1,Az2,fn_vars) !!!!!!!!!!!
% %
% % - num mats = the number of outputs given by funct
% 
% % - When M_{m,n} = M_{-m,n} = M_{-m,-n} = M_{m,-n} only!
% 
% Mat = zeros(2*Az_Dim+1,2*Az_Dim+1,num_mats);
% 
% v1 = 1:Dim; loop_Az1 = -Az_Dim; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;
% 
% for loop1=1:Az_Dim
%  v2 = 1:Dim; loop_Az2 = -Az_Dim; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;   
%  for loop2=1:Az_Dim  
%   Mat(v1,v2,:) = feval(funct,loop_Az2,-loop_Az1,fn_vars);  
%   % - use the ability to swap (m,n)->(-m,-n)->(-m,n)->(m,-n)
%   Mat(v1_flip,v2_flip,:) = Mat(v1,v2,:);
%   Mat(v1,v2_flip,:) = Mat(v1,v2,:);
%   Mat(v1_flip,v2,:) = Mat(v1,v2,:);
%   loop_Az2 = loop_Az2+1;
%   v2 = v2 + 1*Dim; v2_flip = v2_flip - 1*Dim;
%  end
%  Mat(v1,v2,:) = feval(funct,loop_Az2,-loop_Az1,fn_vars);   
%  % - use the ability to swap (m,n)->(-m,n)
%  Mat(v1_flip,v2,:) = Mat(v1,v2,:);
%  v1 = v1+1*Dim; loop_Az1 = loop_Az1+1;
%  v1_flip = v1_flip-1*Dim;
% end 
% % - last row
% v2 = 1:Dim; loop_Az2 = -Az_Dim; 
% v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; 
% for loop2=1:Az_Dim   
%  Mat(v1,v2,:) = feval(funct,loop_Az2,-loop_Az1,fn_vars);   
%  % - use the ability to swap (m,n)->(m,-n)
%  Mat(v1,v2_flip,:) = Mat(v1,v2,:);
%  loop_Az2 = loop_Az2+1;
%  v2 = v2 + 1*Dim; v2_flip = v2_flip - 1*Dim;
% end
% % - and finally the (m,n)=(0,0) entry 
% Mat(v1,v1,:) = feval(funct,loop_Az2,-loop_Az1,fn_vars); 
% 
% return
% 
% % - this uses numerical integration
% 
% function Mat = fn_FillMat_SymmA(Vals, Dim, Az_Dim, theta_vec)
% 
% % - When M_{m,n} = M_{-m,n} = M_{-m,-n} = M_{m,-n} only!
% 
% Mat = zeros(2*Az_Dim+1);
% 
% v1 = 1:Dim; loop_Az1 = -Az_Dim; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;
% 
% for loop1=1:Az_Dim
%  v2 = 1:Dim; loop_Az2 = -Az_Dim; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;   
%  for loop2=1:Az_Dim  
%   Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%   % - use the ability to swap (m,n)->(-m,-n)->(-m,n)->(m,-n)
%   Mat(v1_flip,v2_flip) = Mat(v1,v2);
%   Mat(v1,v2_flip) = Mat(v1,v2);
%   Mat(v1_flip,v2) = Mat(v1,v2);
%   loop_Az2 = loop_Az2+1;
%   v2 = v2 + 1*Dim; v2_flip = v2_flip - 1*Dim;
%  end
%  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%  % - use the ability to swap (m,n)->(-m,n)
%  Mat(v1_flip,v2) = Mat(v1,v2);
%  v1 = v1+1*Dim; loop_Az1 = loop_Az1+1;
%  v1_flip = v1_flip-1*Dim;
% end 
% % - last row
% v2 = 1:Dim; loop_Az2 = -Az_Dim; 
% v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; 
% for loop2=1:Az_Dim   
%  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%  % - use the ability to swap (m,n)->(m,-n)
%  Mat(v1,v2_flip) = Mat(v1,v2);
%  loop_Az2 = loop_Az2+1;
%  v2 = v2 + 1*Dim; v2_flip = v2_flip - 1*Dim;
% end
% % - and finally the (m,n)=(0,0) entry 
% Mat(v1,v1) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% 
% return
% 
% function Mat = fn_FillMat_SymmAB(Vals, Dim, Az_Dim, theta_vec)
% 
% % - When M_{m,n} = M_{-m,n} = M_{-m,-n} = M_{m,-n} 
% %   and M_{m,n} = (-1)^{n-m}M_{-m,-n} => M_{m,n}=0 (n-m \in 2\mathbb{Z}+1)
% 
% Mat = zeros(2*Az_Dim+1);
% 
% v1 = 1:Dim; loop_Az1 = -Az_Dim; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;
% 
% for loop1=1:(Az_Dim+1)/2
%  v2 = 1:Dim; loop_Az2 = -Az_Dim; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;   
%  for loop2=1:(Az_Dim+1)/2   
%   Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%   % - use the ability to swap n->-1 and m->-m
%   Mat(v1_flip,v2) = Mat(v1,v2); 
%   Mat(v1_flip,v2_flip) = Mat(v1,v2);
%   Mat(v1,v2_flip) = Mat(v1,v2);
%   loop_Az2 = loop_Az2+2;
%   v2 = v2 + 2*Dim; v2_flip = v2_flip - 2*Dim;
%  end
% %  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% %  % - use the ability to swap n and m
% %  Mat(v2,v1) = Mat(v1,v2); Mat(v2,v1_flip) = Mat(v1,v2); Mat(v1_flip,v2) = Mat(v1,v2);
%  v1 = v1+2*Dim; loop_Az1 = loop_Az1+2;
%  v1_flip = v1_flip-2*Dim;
% end 
% % - and finally the (m,n)=(0,0) entry 
% % Mat(v1,v1) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% 
% v1 = 1:Dim; v1 = v1+Dim; loop_Az1 = -Az_Dim+1; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; v1_flip = v1_flip-Dim;
% 
% for loop1=1:Az_Dim/2
%  v2 = 1:Dim; v2 = v2+Dim; loop_Az2 = -Az_Dim+1; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; v2_flip = v2_flip-Dim;   
%  for loop2=1:Az_Dim/2   
%   Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%   % - use the ability to swap n->-1 and m->-m
%   Mat(v1_flip,v2) = Mat(v1,v2); 
%   Mat(v1_flip,v2_flip) = Mat(v1,v2);
%   Mat(v1,v2_flip) = Mat(v1,v2);
%   loop_Az2 = loop_Az2+2;
%   v2 = v2 + 2*Dim; v2_flip = v2_flip - 2*Dim;
%  end
%  v1 = v1+2*Dim; loop_Az1 = loop_Az1+2;
%  v1_flip = v1_flip-2*Dim;
% end 
% 
% loop_Az1=0; loop_Az2=0;
% v1 = Az_Dim*Dim+1:(Az_Dim+1)*Dim; v1_flip=v1; v2=v1; v2_flip=v1;
% 
% % - the (m,n)=(0,0) entry 
% Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% v1 = v1-2*Dim; loop_Az1 = loop_Az1-2;
% v1_flip = v1_flip+2*Dim;
% for loop1=1:Az_Dim/2 % - work out from centre going up
%  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
%  % - and left
%  Mat(v2,v1) = get_G_mat(loop_Az1,-loop_Az2,Vals,theta_vec,Dim);
%  % - use the ability to swap (m,n) -> (-m,n) or (m,-n)
%  Mat(v1_flip,v2) = Mat(v1,v2); Mat(v2,v1_flip) = Mat(v2,v1); 
%  v1 = v1-2*Dim; loop_Az1 = loop_Az1-2;
%  v1_flip = v1_flip+2*Dim;
% end 
% 
% return
% 
% function Mat = fn_FillMat_SymmABC(Vals, Dim, Az_Dim, theta_vec)
% 
% % - When M_{m,n} = M_{-m,n} = M_{-m,-n} = M_{m,-n} 
% %   and M_{m,n} = M_{n,m}
% %   and M_{m,n} = (-1)^{n-m}M_{-m,-n} => M_{m,n}=0 (n-m \in 2\mathbb{Z}+1)
% 
% Mat = zeros(2*Az_Dim+1);
% 
% v1 = 1:Dim; loop_Az1 = -Az_Dim; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;
% 
% for loop1=1:(Az_Dim+1)/2
%  v2 = 1:Dim; loop_Az2 = -Az_Dim; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;   
%  for loop2=1:(Az_Dim+1)/2   
%   Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%   % - use the ability to swap n->-1 and m->-m
%   Mat(v1_flip,v2) = Mat(v1,v2); 
%   Mat(v1_flip,v2_flip) = Mat(v1,v2);
%   Mat(v1,v2_flip) = Mat(v1,v2);
%   loop_Az2 = loop_Az2+2;
%   v2 = v2 + 2*Dim; v2_flip = v2_flip - 2*Dim;
%  end
% %  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% %  % - use the ability to swap n and m
% %  Mat(v2,v1) = Mat(v1,v2); Mat(v2,v1_flip) = Mat(v1,v2); Mat(v1_flip,v2) = Mat(v1,v2);
%  v1 = v1+2*Dim; loop_Az1 = loop_Az1+2;
%  v1_flip = v1_flip-2*Dim;
% end 
% % - and finally the (m,n)=(0,0) entry 
% % Mat(v1,v1) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% 
% v1 = 1:Dim; v1 = v1+Dim; loop_Az1 = -Az_Dim+1; 
% v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; v1_flip = v1_flip-Dim;
% 
% for loop1=1:Az_Dim/2
%  v2 = 1:Dim; v2 = v2+Dim; loop_Az2 = -Az_Dim+1; 
%  v2_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim; v2_flip = v2_flip-Dim;   
%  for loop2=1:Az_Dim/2   
%   Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);  
%   % - use the ability to swap n->-1 and m->-m
%   Mat(v1_flip,v2) = Mat(v1,v2); 
%   Mat(v1_flip,v2_flip) = Mat(v1,v2);
%   Mat(v1,v2_flip) = Mat(v1,v2);
%   loop_Az2 = loop_Az2+2;
%   v2 = v2 + 2*Dim; v2_flip = v2_flip - 2*Dim;
%  end
%  v1 = v1+2*Dim; loop_Az1 = loop_Az1+2;
%  v1_flip = v1_flip-2*Dim;
% end 
% 
% loop_Az1=0; loop_Az2=0;
% v1 = Az_Dim*Dim+1:(Az_Dim+1)*Dim; v1_flip=v1; v2=v1; v2_flip=v1;
% 
% % - the (m,n)=(0,0) entry 
% Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
% v1 = v1-2*Dim; loop_Az1 = loop_Az1-2;
% v1_flip = v1_flip+2*Dim;
% for loop1=1:Az_Dim/2 % - work out from centre going up
%  Mat(v1,v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
%  % - use the ability to swap (m,n) -> (n,m) and (m,n) -> (-m,n) or (m,-n)
%  Mat(v2,v1) = Mat(v1,v2); Mat(v2,v1_flip) = Mat(v1,v2); Mat(v1_flip,v2) = Mat(v1,v2);
%  v1 = v1-2*Dim; loop_Az1 = loop_Az1-2;
%  v1_flip = v1_flip+2*Dim;
% end 
% 
% return

function Mat = fn_FillMat_SymmD(Vals, Dim, Az_Dim, theta_vec)

% - When M_{m,n} = M_{-m,-n} 
%   and M_{m,n} = M_{n,m}
%   and M_{m,n} = (-1)^{n-m}M_{-m,-n} => M_{m,n}=0 (n-m \in 2\mathbb{Z}+1)

Mat = zeros(Dim*(2*Az_Dim+1));

v1 = 1:Dim; loop_Az1 = -Az_Dim; 
v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;

for loop1=1:Az_Dim+1
 Mat(v1, v1) = get_G_mat(loop_Az1,-loop_Az1,Vals,theta_vec,Dim);
 if loop_Az1 ~= 0
  Mat(v1_flip, v1_flip)=Mat(v1, v1);
  Mat(v1, v1_flip) = get_G_mat(-loop_Az1,-loop_Az1,Vals,theta_vec,Dim);
  Mat(v1_flip, v1) = Mat(v1, v1_flip);
 end
  % - and in-between
 v2 = v1+2*Dim; loop_Az2 = loop_Az1+2;
 v2_flip = v1_flip-2*Dim;
 for loop2=1:(2*abs(loop_Az1)-1)/2
  Mat(v1, v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
  % - use to get LEFT TRI
  Mat(v2, v1) = Mat(v1, v2);
  % - use to get BOTTOM & RIGHT TRI'S
  Mat(v1_flip, v2_flip) =  ((-1)^(loop_Az1+loop_Az2))*Mat(v1, v2);
  Mat(v2_flip, v1_flip) =  ((-1)^(loop_Az1+loop_Az2))*Mat(v2, v1);
  
  loop_Az2 = loop_Az2+2;
  v2 = v2 + 2*Dim; v2_flip = v2_flip - 2*Dim;
 end
 v1 = v1+Dim; loop_Az1 = loop_Az1+1;
 v1_flip = v1_flip-Dim;
end 

return

function Mat = fn_FillMat_SymmE(Vals, Dim, Az_Dim, theta_vec)

% - When M_{m,n} = M_{-n,-m} 
%   and M_{m,n} = (-1)^{n-m}M_{n,m} 

Mat = zeros(Dim*(2*Az_Dim+1));

v1 = 1:Dim; loop_Az1 = -Az_Dim; 
v1_flip = 2*Az_Dim*Dim+1:(2*Az_Dim+1)*Dim;

for loop1=1:Az_Dim+1
 % - DIAG TOP LEFT TO BOTTOM RIGHT - %   
 Mat(v1, v1) = get_G_mat(loop_Az1,-loop_Az1,Vals,theta_vec,Dim);
 if loop_Az1 ~= 0
  Mat(v1_flip, v1_flip)=Mat(v1, v1);
 end
 % - DIAG TOP RIGHT TO BOTTOM LEFT - %   
 Mat(v1, v1_flip) = get_G_mat(loop_Az1,loop_Az1,Vals,theta_vec,Dim);
 if loop_Az1 ~= 0
  Mat(v1_flip, v1)=Mat(v1, v1_flip);
 end
 % - and in-between
 v2 = v1+Dim; loop_Az2 = loop_Az1+1;
 v2_flip = v1_flip-Dim;
 for loop2=1:2*abs(loop_Az1)-1
  Mat(v1, v2) = get_G_mat(loop_Az2,-loop_Az1,Vals,theta_vec,Dim);
  % - use to get LEFT TRI
  Mat(v2, v1) = ((-1)^(loop_Az1+loop_Az2))*Mat(v1, v2);
  % - use to get BOTTOM & RIGHT TRI'S
  Mat(v1_flip, v2_flip) =  Mat(v2, v1);
  Mat(v2_flip, v1_flip) =  Mat(v1, v2);
  
  loop_Az2 = loop_Az2+1;
  v2 = v2 + Dim; v2_flip = v2_flip - Dim;
 end
 v1 = v1+Dim; loop_Az1 = loop_Az1+1;
 v1_flip = v1_flip-Dim;
end 

return
%

function Mat = fn_FillMat_NoSymm(Vals, Dim, Az_Dim, theta_vec)

% - No symmetry of matrix

Mat = zeros(Dim*(2*Az_Dim+1));

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

%

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

%

function L = fn_LogSingmod(theta_vec,nn)

% INTEGRAL OF -LOG[1-COS[PSI]][1-COS[nn*PSI] dPSI PSI IN [0,PI]

res = (length(theta_vec)-1)/2;
theta_vec = theta_vec(res+1:end);

L_vals = zeros(1,res+1);

L_vals(2:res+1) = -log(1-cos(theta_vec(2:res+1))).*(1-cos(nn*theta_vec(2:res+1)));

dt = pi/res;

L = L_vals(1) + L_vals(res+1);

L = L + 2*sum(L_vals(2:res));

L = L*dt/2;

return

