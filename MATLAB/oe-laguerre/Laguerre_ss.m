% State-space matrices A and B of Laguerre-filter expansion
%
function [A_L,B_L]=Laguerre_ss(alpha,M)
%
beta=sqrt(1-alpha*alpha);
%
if M==1
    A_L=alpha;
    B_L=beta;
end
%
if M==2
    A_L=[alpha 0; 1-alpha*alpha alpha];
    B_L=[beta; -alpha*beta];
end
%
if M>2
    A_new=[1-alpha*alpha alpha];
    A_L=[alpha 0; A_new];
    B_new=-alpha*beta;
    B_L=[beta; B_new];
%
    for k=3:M
        A_new=[-alpha*A_new(1) A_new];
        A_L=[A_L zeros(size(A_L,1),1); A_new];
        B_new=-alpha*B_new;
        B_L=[B_L; B_new];
    end
end
%
% end