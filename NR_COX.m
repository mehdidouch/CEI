function NR_COX = NR_COX(T, Z, t_censure, beta_start, stop_criterion)

%

beta = beta_start;

compteur = 0;

% stop_criterion < sum(square(score_cox(beta, T, Z, t_censure))) && compteur <  10
while compteur <  10
    
    beta = beta + (diffScore(beta,T,Z, t_censure)^-1)*score_cox(beta,T,Z, t_censure);
    compteur = compteur + 1;
    
end

NR_COX = beta;

end