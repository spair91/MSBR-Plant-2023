%Populate lookup table for Reheater enthalpy of vaporization. 
%Function uses XSteam.m to take an array of possible pressure in the 
%reheater. XSteam provides enthalpy of sat vapor and sat fluid 
%at given pressure. Difference between the enthalpies provides the enthalpy
%of vaporization. This lookup table will be used to find the enthalpy of 
%vaporization of fluid in reheater when running sim.

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
H_rhvap_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'hV_p' and 'hL_p' for enthalpy of sat vapor and
    % fluid, respectively. Values are found in XSteam by pressure in P_rh_table.
    % Find the difference of the sat vapor enthalpy and sat fluid enthalpy.
    % Save the result in each row of the result matrix
    H_rhvap_table(1,i) = (XSteam('hV_p', P_rh_table1(i).*1e1) - XSteam('hL_p', P_rh_table1(i).*1e1))./1e3; % Pressure from MPa to bar, Enthalpy from kJ/kg to MJ/kg
end