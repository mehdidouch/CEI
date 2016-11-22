%% Application regression COX

n = 100;

T = zeros(n,1);

Z = randn(n, 3);

 beta = [10; -2; 3];
 
 lambda = exp(Z * beta);
 

for i=1:n 
    
    T(i) = exprnd(1/lambda(i));
    
end

t_censor = median(T);
bool_censor = T>t_censor;

% Value of beta given by fminsearch
[beta_opt,fval] = fminsearch(@(x) -vraisemblance_cox(x,T,Z,t_censor),[0;0;0]);

% We have implemented our own Newton-Raphson
beta_opt_NR = NR_COX(T, Z, t_censor, [10;-2;3], 0.01);

% Cox proportional hazards regression by Matlab
b = coxphfit(Z,T,'censoring',bool_censor);
