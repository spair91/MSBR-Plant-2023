%Populate lookup table for Nozzle chest pressure. Function uses XSteam.m to
%take an array of possible densities in the nozzle chest and an array of
%enthalpies in the nozzle chest. XSteam.m provides the pressure at the
%density and entalpy given. This lookup table will be used to find nozzle
%chest pressure from the calculated enthalpy and density when running sim.

%Function based off of the same function used in the University of
%Tennessee paper with differences to use XSteam for values
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.

% density of steam in the nozzle chest in kg/m^3
rho_table = [10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 125, 150];

% Enthalpy of steam in the nozzle chest
% in kJ/kg
H_inp = [1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400];
H_table = H_inp/1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_table);
n = length(H_inp);

% matrix of zeros to store result
P_nc_table = zeros(n, m);

for i=1:m
    %iterate through densities 
    for j=1:n
        %iterate through enthalpies with density constant
        %call XSteam function 'p_hrho' for each value of enthalpy
        %save the result in each row of the result matrix
        P_nc_table(j,i) = XSteam('p_hrho', H_inp(j), rho_table(i))*1e-1; % bar to MPa
    end
end