function V = vraisemblance_cox(beta, T, Z, t_censure)
%V returns the cox-likelihood for a given value of the input parameters% This function takes the following parameters:
%   - beta: a vector
%   - Z : covariables matrix
%   - T : event time 
%   - t_censure : limit of event time considered

%Generation of Y

L = length(T);

delta_censure = T>t_censure;

V = 1;

for i=1:L
    if ~delta_censure(i)
        V = V*exp(Z(i,:)*beta)/((T(i)<= T)'*exp(Z*beta));
    end
    
end

end
