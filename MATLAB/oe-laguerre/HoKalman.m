function [A,B,C,D]=HoKalman(h,s,tol)
%   HOKALMAN(h,s) returns a discrete-time, minimal, state-space realization 
%   (A,B,C,D) with McMillan degree of at most s, from an impulse response h.

    if nargin<3
        tol=0.01;
    end
    
    H = hankel(h(2:s+1),h(s+1:end));
    [U,S,V] = svd(H);

    nhat=sum(diag(S)>=tol*S(1));
    
    T=eye(s);

    Psi_o = U(:,1:s)*S(1:s,1:s)^(1/2)*T;
    Psi_c = inv(T)*S(1:s,1:s)^(1/2)*V(:,1:s)';

    A=pinv(Psi_o(1:end-1,1:nhat))*Psi_o(2:s,1:nhat);
    B=Psi_c(1:nhat,1);
    C=Psi_o(1,1:nhat);
    D=h(1);

end