%Populate lookup table for Reheater Pressure. 
%Function uses XSteam to take an array of possible densities and enthalpies
%of the fluid in the reheater. XSteam provides pressure given an enthalpy 
%and density. This lookup table will be used to find the pressure in the 
%reheater vaporization of fluid in reheater when running sim.

%Function based off of the same function used in the University of
%Tennessee paper with differences to use XSteam for values
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.

% density of steam in the reheater in kg/m^3
rho_rh_table = [1, 2.5, 5, 7.5, 10, 15, 20];

% Enthalpy of steam in the reheater
H_rh_inp = [100, 200, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500]; %kJ/kg
H_rh_table2 = H_rh_inp./1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_rh_table);
n = length(H_rh_inp);

% matrix of zeros to store result
P_rh_table = zeros(n, m);

for i=1:m
    for j=1:n
        % Call XSteam function 'p_hrho' for each value of enthalpy
        % while holding density constant
        % Save the result in each row of the result matrix
        P_rh_table(j,i) = XSteam('p_hrho', H_rh_inp(j), rho_rh_table(i))*1e-1; % bar to MPa
    end
end