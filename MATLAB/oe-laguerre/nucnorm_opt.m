function [g,h] = nucnorm_opt(y,X,gamma,W1,W2)
global s m alpha

cvx_begin sdp quiet
    variables g(m,1) h(m,1)
        minimize( norm_nuc(W1*hankel(h(2:s+1),h(s+1:end))*W2) + gamma*sum_square(y - X*g) )
        subject to
        h(1) == alpha/sqrt(1-alpha^2)*g(1);
        for k=2:m
            h(k) == 1/sqrt(1-alpha^2)*g(k-1) + alpha/sqrt(1-alpha^2)*g(k);
        end
cvx_end

end

