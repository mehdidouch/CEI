function I = Fischer(beta, T, Z, t_censure)
%V returns the score of the cox-likelihood for a given value of the input parameters%
%This function takes the following parameters:
%   - beta: a vector
%   - Z : covariables matrix
%   - T : event time
%   - t_censure : limit of event time considered

%Generation of Y

L = length(T);

delta_censure = T > t_censure;

I = zeros(length(beta),length(beta));

ARG= Z*beta;

for k = 1 : length(beta)
    for l =  1 : length(beta)
        for i=1:L
            if ~delta_censure(i)
                
                terme_denominateur = ((T(i)<= T)'*exp(ARG));
                terme_1 = 0;
                terme_2_1 = 0;
                terme_2_2 = 0;
                
                for j=1:L
                    t = T(i)<= T;
                    terme_1 = terme_1 + t(j) * Z(j,k) * Z(j,l) * exp(ARG(j))*terme_denominateur;
                    terme_2_1 = terme_2_1 +  t(j) * Z(j,l) * exp(ARG(j));
                    terme_2_2 = terme_2_2 +  t(j) * Z(j,k) * exp(ARG(j));
                end
                
                I(k,l) = I(k,l) - (terme_1 - terme_2_1 * terme_2_2) / (terme_denominateur^2);
            end
        end
        
    end  
    I = - I;
end