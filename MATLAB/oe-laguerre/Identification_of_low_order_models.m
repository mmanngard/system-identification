clear all, close all, clc
%Identification of low-order models using Laguerre basis function expansions.
%Presented at the 18th Symposium on System Identification, 2018, Stockholm.
%M. Manngård and H. T. Toivonen, Åbo Akademi University, 16.7.2018.

global s m alpha
warning('off','all')

%% Generate data
%y(t) = G(q)u(t) + F(q)e(t)
F = tf([1 -1.38 0.4],[1 -1.9 0.91],1,'variable','z^-1');
G = tf([0 0.2530 -0.9724 1.4283 -0.9493 0.2410],[1 -4.15 6.8831 -5.6871 2.3333 -0.3787],1,'variable','z^-1');

u = randn(1200,1);                                                          %input
e = sqrt(0.0025)*randn(1200,1);                                             %noise
y = lsim(G,u) + lsim(F,e);                                                  %output

%% Expand the state equations in terms of a Laguerre basis
m = 20;                                                                     %number of basis functions
alpha = 0.7;                                                                %pole
[AL,BL] = laguerre_ss(alpha,m);                                             %state-space matrices A and B of Laguerre-filter expansion
X = lsim(ss(AL,BL,eye(m),0,1),u);                                           %input to state


%% Nuclear norm minimization
gamma = 5;                                                                  %weight parameter
s = 10;                                                                     %maximum system order
W1=eye(s); W2=eye(m-s);                                                     %nuclear norm weights  

tic
for it=1:20                                                                 %iterative reweighting (Mohan and Fazel, 2010)   
    %cvx_solver mosek_2
    [g,h] = nucnorm_opt(y,X,gamma,W1,W2);                                   %optimization (requires CVX)
    [W1,W2] = nucnorm_weights(W1,W2,h,1e-4);                                %update weights
end
toc

%% Compute Laguerre-domain state-space  realization
[At,Bt,Ct,Dt] = HoKalman(h,s);                                              %Ho-Kalman algorithm (Ho and Kalman, 1966)

%% Inverse Laguerre transform
I=eye(size(At));

A = inv(I+alpha*At)*(alpha*I+At);
B = sqrt(1-alpha^2)*inv(I+alpha*At)*Bt;
C = sqrt(1-alpha^2)*Ct*inv(I+alpha*At);
D = Dt - alpha*Ct*inv(I+alpha*At)*Bt;

%identified model
sys = ss(A,B,C,D,1);

%%
% References:
% Ho, B. and Kalman, R.E. (1966). Effective construction of linear 
% 	state-variable models from input/output functions. Regelungstechnik,
% 	14(1-12), 545-548.
% Manng�rd, M. and Toivonen, H.T. (2018). Identification of low-order models
%   using Laguerre basis function expansions. IFAC-PapesOnline, XXX-XXX.
% Mohan, K. and Fazel, M. (2010). Reweighted nuclear norm minimization with 
% 	application to system identification. In Proceedings of the American 
% 	Control Conference (ACC), 2953-2959.
