function U = score_cox(beta, T, Z, t_censure)
%V returns the score of the cox-likelihood for a given value of the input parameters% 
%This function takes the following parameters:
%   - beta: a vector
%   - Z : covariables matrix
%   - T : event time 
%   - t_censure : limit of event time considered

%Generation of Y

L = length(T);

delta_censure = T > t_censure;

U = zeros(1,length(beta));

for i=1:L
    if ~delta_censure(i)
        U = U + Z(i,:) - ((Z.'* ((T(i)<= T).* exp(Z*beta))) ./ ((T(i)<= T)'*exp(Z*beta))).';
    end
    
end

U = U.';


end