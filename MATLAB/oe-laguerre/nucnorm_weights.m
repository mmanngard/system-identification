function [W1,W2] = nucnorm_weights(W1,W2,h,eps)
%NUCNORM_WEIGHTS updates the weights for the iterative nuclear norm
%minimization method (Mohan and Fazel, 2010).
global s m
    
    [U,S,V] = svd(W1*hankel(h(2:s+1),h(s+1:end))*W2,'econ'); S=S(:,1:s);
    Y = inv(W1)*U*S*U'*inv(W1);
    Z = inv(W2)*V*S*V'*inv(W2);
    W1 = real((Y+eps*eye(s))^-0.5);
    W2 = real((Z+eps*eye(m-s))^-0.5);  
end

