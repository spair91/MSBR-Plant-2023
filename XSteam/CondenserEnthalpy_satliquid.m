%Populate lookup table for Condenser Saturated fluid enthalpy. 
%Function uses XSteam.m to take an array of possible pressure in the 
%condenser. XSteam provides enthalpy of sat fluid at given pressures. 
%This lookup table will be used to find the sat. fluid enthalpy of 
%the condenser when running sim.


%Function based off of the same function used in the University of
%Tennessee paper with differences to use XSteam for values
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.


% Pressure of saturated steam in MPa
P_con_table = [1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 5e-2, 1e-1];

% length of each array
m = length(P_con_table);

% matrix of zeros to store result
H_con_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'hL_p' for each value of enthalpy
    % with condenser pressure input in bars. 
    % Save the result in each row of the result matrix
    H_con_table(1,i) = XSteam('hL_p', P_con_table(i).*1e1)./1e3; % Pressure from MPa to bar, Enthalpy from kJ/kg to MJ/kg
end