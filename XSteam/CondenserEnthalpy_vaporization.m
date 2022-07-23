%Populate lookup table for Condenser enthalpy of vaporization. 
%Function uses XSteam.m to take an array of possible pressure in the 
%condenser. XSteam provides enthalpy of sat fluid and sat vapor
%at given pressures and the function uses this to find enthalpy difference. 
%This lookup table will be used to find the enthalpy of vaporization of 
%the condenser when running sim.

%Function is based off of University of Tennessee paper for entire-plant
%model with differences to use XSteam for most values.
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.

% Pressure of saturated steam in MPa
P_con_table = [1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 5e-2, 1e-1];

% length of each array
m = length(P_con_table);

% matrix of zeros to store result
H_convap_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'hV_p' and 'hL_p' for each value of enthalpy
    % with given pressure in bar. Find difference between hV_p and hL_p.
    % Convert kJ/kg to MJ/kg and save the result in each row of the result matrix
    H_convap_table(1,i) = (XSteam('hV_p', P_con_table(i).*1e1) - XSteam('hL_p', P_con_table(i).*1e1))./1e3; % Pressure from MPa to bar, Enthalpy from kJ/kg to MJ/kg
end