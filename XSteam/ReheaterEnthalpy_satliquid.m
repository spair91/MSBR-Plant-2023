%Populate lookup table for Reheater enthalpy of saturated fluid. 
%Function uses XSteam.m to take an array of possible pressure in the 
%reheater. XSteam provides enthalpy of reheater sat fluid at given pressure
%This lookup table will be used to find reheater enthalpy
%when running sim.

%Function based off of the same function used in the University of
%Tennessee paper with differences to use XSteam for values
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.

% Pressure of saturated steam in MPa
P_rh_table1 = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0, 12.5, 15.0];

% length of each array
m = length(P_rh_table1);

% matrix of zeros to store result
H_rh_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'hL_p' for each value of enthalpy
    % while holding density constant
    % Save the result in each row of the result matrix
    H_rh_table(1,i) = XSteam('hL_p', P_rh_table1(i).*1e1)./1e3; % Pressure from MPa to bar, Enthalpy from kJ/kg to MJ/kg
end