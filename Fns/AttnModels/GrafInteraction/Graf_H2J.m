function T=Graf_H2J(k,dists,rots,N)

% - 01.06.10 modified from Fabien's coor_trans_mat by LB

% - M = vert dim
% - N = az_dim

% - In the below note that my M=Az_Dim!!!!!!!! - %

% - the incident wave for plate j is the sum of the corresponding ambient
% inc wave component and the scatterered waves from all other plates
% - Thus, if phiIj = sum_m Jm(k0*rj)*Dj_{m}*exp(i*m*thj)
% - we have Dj_{u} = til{D}j_{u} + sum_{p~=j}Tjp_{u}*ap
% - where ap = [ap_{-M},...ap_{M}] -> (length Vert_Dim*(2Az+1)) scattered waves for plate p
% ( assuming that phiSp=sum_m Hm(k*rp)*ap_{m}*exp(i*m*thp) )
% - Tjp_{u} is a matrix (size Vert_Dim x Vert_Dim*(2Az+1)) derived from Graf's formula
% - Tjp_{u} = [Tjp_{u,-M},...,Tjp_{u,M}]
% - Tjp_{u,m}=H_{m-u}(k*dist_pj)*exp(i*ang_pj*(m-u))

% - so say Dj=[Dj_{-M},...,Dj_{M}]=til{D}j + sum_{p~=j} Tjp*ap
%   Tjp = [Tjp_{-M};...;Tjp_{M}] is size Vert_Dim*(2Az+1)-squared

% - so D=[D1,..DNp]=til{D}+T*a=til{D}+T*[a1,...aNp]
%   T = [T11,T12,...,T1Np;...;TNp1,T12,...,TNpNp] is size Np*Vert_Dim*(2Az+1)-squared
%   and Tmm=0

if ~exist('M','var'); M=1; end

np=length(dists);
T=zeros(np*(2*N+1)*M,np*(2*N+1)*M);

v1 = 1:M*(2*N+1);
for loop_p1=1:np % - j
 v2 = 1:M*(2*N+1);
 for loop_p2=1:np % - p
  if loop_p1~=loop_p2 % - only needed if p~=j
   dist = dists(loop_p1,loop_p2);
   rot = rots(loop_p1,loop_p2);
   Tjp=zeros((2*N+1)*M);
    u1=1:M;
    for loop_Az1=-N:N % - u
     u2=1:M;   
     for loop_Az2=-N:N % - m
      Tjp(u1,u2)=besselh(loop_Az2-loop_Az1,k*dist)...
        *exp(1i*(loop_Az2-loop_Az1)*rot);
      u2=u2+M;            
     end
     u1=u1+M;
    end
    T(v1,v2)=Tjp;
  end
  v2 = v2+M*(2*N+1); 
 end
 v1 = v1+M*(2*N+1);
end

return