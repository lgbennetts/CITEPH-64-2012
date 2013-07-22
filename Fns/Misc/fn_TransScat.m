function Mat_Out = fn_TransScat(Mat_In, opt)

% function Mat_Out = TransScat(Mat_In, opt)
% - opt = -1 Scat -> Trans
% - opt =  1 Trans -> Scat

% - calc trans from scat and vice versa

N = length(Mat_In)/2;

if opt == -1 %- trans from scat

Rm = Mat_In(1:N, 1:N);
Tm = Mat_In(N+1:2*N,1:N);
Rp = Mat_In(N+1:2*N,N+1:2*N);
Tp = Mat_In(1:N,N+1:2*N); invTp = inv(Tp);

Mat_Out = zeros(2*N);

Mat_Out(1:N, 1:N) = Tm - Rp*invTp*Rm;
Mat_Out(N+1:2*N,1:N) = - invTp*Rm;
Mat_Out(N+1:2*N,N+1:2*N) = invTp;
Mat_Out(1:N,N+1:2*N) = Rp*invTp;

elseif opt == 1 % - scat from trans
    
P11 = Mat_In(1:N, 1:N);
P21 = Mat_In(N+1:2*N,1:N);
P22 = Mat_In(N+1:2*N,N+1:2*N); invP22 = inv(P22);
P12 = Mat_In(1:N,N+1:2*N);    

Mat_Out(1:N, 1:N) = -invP22*P21;
Mat_Out(N+1:2*N,1:N) = P11 - P12*invP22*P21;
Mat_Out(N+1:2*N,N+1:2*N) = P12*invP22;
Mat_Out(1:N,N+1:2*N) = invP22;

end

return