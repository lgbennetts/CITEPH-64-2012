function [X_new,D_new,q_new] = fn_ArrangeTrans(X,D)

tol = 1e-4;

q = log(diag(D))/1i; 

vR = find(abs(imag(q))<tol);

vI = find(abs(imag(q))>=tol);

[qR,iR]=sort(real(q(vR)));

[qI,iI]=sort(imag(q(vI)));

qR=flipud(qR); qI=flipud(qI); iR=flipud(iR); iI=flipud(iI);

qRp = qR(1:length(qR)/2); qRm = qR(1+length(qR)/2:end);
iRp = iR(1:length(iR)/2); iRm = iR(1+length(iR)/2:end);
qIp = qI(1:length(qI)/2); qIm = qI(1+length(qI)/2:end);
iIp = iI(1:length(iI)/2); iIm = iI(1+length(iI)/2:end);

q_new = [qRp;1i*qIm;flipud(qRm);flipud(1i*qIp)];

ind = [vR(iRm);vI(iIm);flipud(vR(iRp));flipud(vI(iIp))];

D_new = diag(exp(1i*q_new));

X_new = X(:,ind);

return
