function Main_Channel_freq(Np,lam_vec,scatyp,row_seps,file_marker,COMM,POOL)

if ~exist('POOL','var'); POOL=0; end
if ~exist('COMM','var'); COMM=1; end

if POOL
 matlabpool('open',POOL)
end

if ~exist('Np','var'); Np=1; end
if ~exist('extra_pts','var'); extra_pts=[]; end
if ~exist('file_marker','var'); file_marker=0; end
if ~exist('lam_vec','var'); lam_vec=2*pi/0.6; end %4.0317;1.0472;%2*pi-0.1:0.1:2*pi+0.1;0.95*pi/16:0.05*pi/16:1.05*pi/16;
if ~exist('scatyp','var'); scatyp='cyl'; end
if ~exist('row_seps','var'); Rows=1; row_seps = 3 + zeros(1,Rows-1); end

% if file_marker~=0
%  file_name = ['Main_Channel_freq_',scatyp,'_',file_marker,'.mat'];
%  dummy = 1.234;   
%  save(['../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/',file_name],...
%      'dummy'); clear dummy
% end

Rows = length(row_seps)+1;

fortyp = 'waveno';

if COMM
disp('%-----------------------------------------------%')
disp('%---------- START: Main_Channel_freq -----------%')
end

%% Creation of local paths
%CreatePaths_Mac

%% Wavelengths
%%% NB. There are m `tank modes' when k\in \pi(m,m+1)/width

%lam_vec = 1:0.025:3;

%% Parametrisation
% Note that the origin of the Cartesian system is located at the bottom 
% left corner of the wave tank and coincides with the free-surface at rest.

%%% Geometry of the wave tank: [length width height]
TankDim = [40 16 2]; %[15 10 0.5]; %

%%% Separation of the rows
SAME_ROWS = 1; % =0 DIFFERENT ROWS; =1 SAME ROWS-> NO NEED TO RECALC SCATS

%%% Location and geometry of disks: [x_c y_c Rad thick]
for loop=1:Np
 GeomDisks(loop,:)=[4 0+(TankDim(2)/2/Np)+(TankDim(2)/Np)*(loop-1) 0.5 1e-1];
end

%GeomDisks = [7.5 1.5 0.5 0.005; 7.5 8.5 0.5 0.005]; %[4 8 0.5 1e-10; 8 8 0.5 1e-10];

if ~SAME_ROWS
 disp('PUT SOMETHING HERE!!!!!')
 GeomDisks(:,:,1) = GeomDisks(:,:);
 for loop=2:Rows
  GeomDisks(:,:,loop) = GeomDisks(:,:,1);
 end
end

%%% Check disk/disk and disk/wall overlaps
overlap(GeomDisks, TankDim)

%%% Physical properties
Param = ParamDef3d(GeomDisks);

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy
Param = ModParam_def(Param,extra_pts);

parfor loop_lam=1:length(lam_vec)
         
 %%% Forcing
 Forcing = Force_def(Param.g(1), TankDim(3), fortyp, lam_vec(loop_lam));

 %% Solution - Disks in a channel
 % Scattering matrix is S = [Rm Tp;Tm Rp]
 % Wavenumbers matrices v_vec and u_vec are scaled by vertical wavenumbers.
 if SAME_ROWS
  [Rm,Tm,Rp,Tp,v_vec,u_vec,k0,weight_0,x_lim,reson_mkr(loop_lam)] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.g], Param.N, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).';GeomDisks(:,2).'], ...
    Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...
    [Param.res_green, Param.terms_green, Param.cutoff_green, Param.tolres], ...
    Param.extra_pts, scatyp, COMM);
 else
  [Rm,Tm,Rp,Tp,v_vec,u_vec,k0,weight_0,x_lim,reson_mkr(loop_lam)] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.g], Param.N, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4,1).', Param.draft.', GeomDisks(:,3,1).', ...
    [GeomDisks(:,1,1).';GeomDisks(:,2,1).'], ...
    Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...
    [Param.res_green, Param.terms_green, Param.cutoff_green, Param.tolres], ...
    Param.extra_pts, scatyp, COMM);
 end
 
 kv = reshape(diag(k0)*v_vec,1,size(v_vec,1)*size(v_vec,2));
 
 Ep0 = diag(exp(-1i*kv*(x_lim(1)-0) ));
 Rm = Rm*Ep0; Tm = Tm*Ep0;
 Ep = diag(exp( 1i*kv*(x_lim(2)-TankDim(1)) ));
 Rp = Rp*Ep; Tp = Tp*Ep;
 
 v1 = 1:size(v_vec,1)*size(v_vec,2); v2 = v1+size(v_vec,1)*size(v_vec,2);
 
 Full_S = zeros(2*length(v_vec),2*length(v_vec),Rows);
 
 Full_S(v1,v1,1)=Rm; Full_S(v2,v2,1)=Rp;
 Full_S(v1,v2,1)=Tp; Full_S(v2,v1,1)=Tm;
    
%  Full_S_rl(v1,v1,1)=Rm; Full_S_rl(v2,v2,1)=Rp;
%  Full_S_rl(v1,v2,1)=Tp; Full_S_rl(v2,v1,1)=Tm;
    
 II = eye(size(v_vec,1)*size(v_vec,2));
 for loop_Sc=2:Rows
  % - left-to-right method
  Ep = diag(exp(1i*kv*row_seps(loop_Sc-1)));
  
  if ~SAME_ROWS
   [Rm,Tm,Rp,Tp,~,~,~,~,x_lim] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.g], Param.N, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).';GeomDisks(:,2).'], ...
    Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...
    [Param.res_green, Param.terms_green, Param.cutoff_green, Param.tolres], ...
    Param.extra_pts, scatyp, COMM);
  else
    x_lim  = x_lim + row_seps(loop_Sc-1) + (x_lim(2)-x_lim(1));  
  end
  
  R0m = Full_S(v1,v1,loop_Sc-1); T0p = Full_S(v1,v2,loop_Sc-1);
  T0m = Full_S(v2,v1,loop_Sc-1); R0p = Full_S(v2,v2,loop_Sc-1);

  IV1 = II/(II - Rm*Ep*R0p*Ep); IV2 = II/(II - R0p*Ep*Rm*Ep);

  Full_S(v1,v1,loop_Sc) = R0m + T0p*Ep*IV1*Rm*Ep*T0m;
  Full_S(v1,v2,loop_Sc) = T0p*Ep*IV1*Tp;
  Full_S(v2,v2,loop_Sc) = Rp + Tm*Ep*IV2*R0p*Ep*Tp;
  Full_S(v2,v1,loop_Sc) = Tm*Ep*IV2*T0m;

%   if 0
%    % - right-to-left method
%    Ep = diag(exp(1i*kv*row_seps(Rows-loop_Sc+1)));
% 
%    R0m = RR; T0p = TT; T0m = TT; R0p = RR; 
%    Rm = Full_S_rl(v1,v1,loop_Sc-1); Tp = Full_S_rl(v1,v2,loop_Sc-1);
%    Tm = Full_S_rl(v2,v1,loop_Sc-1); Rp = Full_S_rl(v2,v2,loop_Sc-1);
% 
%    IV1 = II/(II - Rm*Ep*R0p*Ep); IV2 = II/(II - R0p*Ep*Rm*Ep);
% 
%    Full_S_rl(v1,v1,loop_Sc) = R0m + T0p*Ep*IV1*Rm*Ep*T0m;
%    Full_S_rl(v1,v2,loop_Sc) = T0p*Ep*IV1*Tp;
%    Full_S_rl(v2,v2,loop_Sc) = Rp + Tm*Ep*IV2*R0p*Ep*Tp;
%    Full_S_rl(v2,v1,loop_Sc) = Tm*Ep*IV2*T0m;
%   end
 end
 
%  Rm = Rm/Ep0; Tm = Tm/Ep0;
%  Ep = diag(exp(-1i*kv*(x_lim(2)-TankDim(1)) ));
%  Rp = Rp*Ep; Tp = Tp*Ep;
 
 R_vec{loop_lam} = Full_S(v1,v1,Rows)/Ep0; 
 T_vec{loop_lam} = Full_S(v2,v1,Rows)/Ep0; 

 v_vecs{loop_lam} = v_vec(1,:); %k_vec(loop_lam) = k0(1);

end

if and(and(~file_marker, length(lam_vec)==1),length(v_vecs{1})<6)
 disp('Rm='); disp(R_vec{1})
 disp('|Rm|='); disp(abs(Rm))
 disp('Tm='); disp(T_vec{1})
 disp('|Tm|='); disp(abs(Tm))
end

if file_marker~=0
    
 w = TankDim(2);
 
 for loopi=1:4
  dum_str{loopi} = num2str(GeomDisks(1,loopi));    
  for loopii=2:Np   
    dum_str{loopi} = [ dum_str{loopi}, ', ' ...
        num2str(GeomDisks(loopii,loopi)) ];
  end
 end
 clear loopi loopii
 dum_str{5} = num2str(Param.E(1));    
  for loopii=2:Np   
    dum_str{5} = [ dum_str{5}, ', ' ...
        num2str(Param.E(loopii)) ];
  end
  clear loopii 

 description = {['Problem: ' scatyp], ...
     [fortyp ' :' num2str(lam_vec(1)) ' -> ' num2str(lam_vec(end)) ' (' ...
      int2str(length(lam_vec)) ')'] ...
     ['Tank: ' num2str(TankDim(1)) 'x' num2str(TankDim(2)) ...
     'x' num2str(TankDim(3))], ...
     [int2str(Np) ' Disk(s): x=' dum_str{1} ...
     '; y=' dum_str{2} ...
     '; r=' dum_str{3} ...
     '; h=' dum_str{4} ...
     '; E=' dum_str{5} ], ...
     ['Vert modes=' int2str(Param.N) ...
     ', extra points=' int2str(Param.extra_pts) ...
     ', angular res=' int2str(Param.res_green) ...
     ', Green fn terms=' int2str(Param.terms_green) ...
     ', Greens fn cutoff=' int2str(Param.cutoff_green) ...
     ', resonance tol=' int2str(Param.tolres)],
     ['Version 1.1: changed Greens fn cutoff to fixed truncation' ...
     'rather than size of terms']};
     
 file_name = ['Main_Channel_freq_',scatyp,'_',file_marker,'.mat'];

 save(['../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/',file_name],...
    'R_vec', 'T_vec', 'lam_vec', 'v_vecs', 'w', 'description','reson_mkr')

else
 beep; beep; beep
end

if POOL; matlabpool close; end

if COMM
disp('%----------- END: Main_Channel_freq ------------%')
disp('%-----------------------------------------------%')
end

return