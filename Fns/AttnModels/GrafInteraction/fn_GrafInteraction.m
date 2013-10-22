% function [D2_body,a2_body,Dt]=fn_GrafInteracion(N,pos,k0)
%
% DESCRIPTION: Transform scattered wavefield for one scatterer to incident
%              wavefield for second scatterer
%
% INPUTS:
%
% N   = number of azimuthal modes
% pos = [x1 y1 R1; x2 y2 R2 ]
%       (x,y) scatterer coordinate system origin
%       R     scatterer radius (enscribed circle)
% k   = free-space wavenumber
% DTM = diffraction transfer matrix for first scatterer
% D0  = amplitudes of incident wave (default is plane wave e^{ikx})
%
% OUTPUTS:
%
% A2I   = amplitudes wave amplitudes on 2nd scatterer are A2I*D0
%         where D0 is vector of ambient wave amplitudes
% A2I   = A2inc + A21
%         where A2inc is transform of incident wave to local coords
%               A21   is scattered wave 
% D2    = A2I*D0
%
% L Bennetts Sept 2013 / Adelaide
%
% REVISION HISTORY: 
% - see Directional_spectrum/Fns/Graf/interac_circ_scat_main
%   and Directional_spectrum/Fns/Graf/interaction_solver_proj2 
% - original code by F Montiel circa 2009                 

function [A2I,A2inc,A21,D2]=fn_GrafInteraction(N,k,pos,DTM,D0)

if ~exist('N','var');   N=2; end
if ~exist('k','var');   k=1; end
if ~exist('pos','var'); pos = [0,0,0.1;1,0,0.1]; end
if ~exist('D0','var');  D0= (1i.^[-N:N])*exp(abs(imag(k*pos(1,3)))); end
D0=reshape(D0,(2*N+1),1);
if ~exist('DTM','var'); DTM = fn_CircScat(N, {'k',k}, pos(1,3)); end

%% COORDINATE SYSYEMS

np=size(pos,1); % = 2 in this case                 

glob=[0,0];

%%% Polar coordinates in relation to global origin 

R0=zeros(1,np); psi0=zeros(1,np); 

for l=1:np
 [R0(l),psi0(l)]=coordinate_change(glob,pos(l,1:2));
end

%%% Mean center of scatterer l wrt local coordinates of scatterer j

R=zeros(np,np); psi=zeros(np,np);

for j=1:np
 for l=1:np
  [R(j,l),psi(j,l)]=coordinate_change(pos(j,:),pos(l,:));        
 end
end

%%% Check that plates don't overlap

for loop1=1:np-1
 for loop2=loop1+1:np
  sumRads = pos(loop1,3)+pos(loop2,3);
  if R(loop1,loop2)<sumRads
   disp(['plates ',num2str(loop1),' & ',num2str(loop2),' overlapping!!!!!'])
   D_body=[]; a_body=[];
   return
  end
  clear sumRads
 end
end

%% INCIDENT WAVE (AMBIENT)

Dt=zeros((2*N+1)*np,(2*N+1));
  
v1 = 1:(2*N+1); 
for l=1:np
 Dt(v1,:) = fn_TransInc(N, R0(l), psi0(l),k);
 v1 = v1+(2*N+1);
end

%% COORDINATE TRANSFORM

% Coordinate transformation matrices T_jl

T=Graf_H2J(k,R,psi,N);

%% INCIDENT WAVE FOR 2ND SCATTERER

v1=1:(2*N+1); v2=v1+(2*N+1); 

A2inc = Dt(v2,v1);

A21   = T(v1,v2)*DTM*Dt(v1,v1);

A2I   = A2inc + A21;

D2    = A2inc*D0;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% - SUBFUNCTIONS  - %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TRANSFORM GENERAL INC WAVE TO LOCAL POLAR COORD SYSTEM

function Mat = fn_TransInc(Az_Dim, dist, rot, k0)

% - Inc wave is phiI=sum Jm(k*r0)*Im*exp(i*m*th0) - %
% - Jn(k0*r0)*exp(i*n*th0) Transforms (using Graf) to
%   sum_m Jm(k*rj)*Amn*exp(i*m*thj)
% - where Amn = J{n-m}(k*R0j)*exp{i*ang0j*(n-m)}
% - so if inc wave for plate j is phiI_(j)=sum Jm(k*rj)*Im_(j)*exp(i*m*thj) - %
% - then I_(j) = {Amn}*I, where I_(j)=[I-M_(j),...,IM_(j)], I=[I-M,...,IM] - %

Vert_Dim=1;

Mat = zeros(Vert_Dim*(2*Az_Dim));
v1=1:Vert_Dim; 
for loop_Az1=-Az_Dim:Az_Dim
 v2=1:Vert_Dim;   
 for loop_Az2=-Az_Dim:Az_Dim     
     Mat(v1,v2)=besselj(loop_Az2-loop_Az1,k0*dist)*exp(1i*rot*(loop_Az2-loop_Az1));     
     v2=v2+Vert_Dim;
 end
 v1=v1+Vert_Dim;
end 

return
