% Calculate the big term from rho_0 equation.
% Big term is the second term in the equation for reactivity required for
% steady-state (rho_0). See ORNL-MSR-67-102 page 4. MATLAB function
% originally from University of Tennessee paper 
% Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.
function bterm=bigterm(bet,lam,t_L,t_C)
    bterm = 0;
    for i = 1:6
        bterm = bterm + bet(i)/(1.0 + ((1.0-exp(-lam(i)*t_L))/(lam(i)*t_C)));
    end
end