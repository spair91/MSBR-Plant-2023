%Populate lookup table for BoP Steam Enthalpy input. Function uses 
%XSteam to take an array of possible steam temperatures and an array of
%possible steam pressures. XSteam provides the enthalpy of steam at the
%given pressure and temperature. This lookup table will be used to find 
%Steam Enthalpy for the steam going into the nozzle chest.

%Function based off of the same function used in the University of
%Tennessee paper with differences to use XSteam for values
%Vikram Singh, Alexander M. Wheeler, Belle R. Upadhyaya, Ondřej Chvála, 
% and M. Scott Greenwood. 2020. Plant-level dynamic modeling of a 
% commercial-scale molten salt reactor system. Nucl. Eng. Des. 360, 
% (Apr, 2020), 110457. DOI: https://doi.org/ 10.1016/j.nucengdes.2019.110457.

temp_table = [275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 800, 900, 1000, 1500, 2000, 2500]; % temp in deg-C
pres_inp = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200, 250, 300, 350, 400, 450, 500, 550]; % pressure in bar
pres_table = pres_inp*1e-1; % convert to MPa

m = length(temp_table);
n = length(pres_table);

% matrix of zeros to store result
H_s_table = zeros(m, n); % Enthalpy of steam

for i=1:m
    for j=1:n
        H_s_table(i, j) = XSteam('h_pT', pres_inp(j), temp_table(i))/1e3; % in MJ/kg
    end
end

% % Correction for error in the data in XSteam
% H_s_table(1 ,5) = H_s_table(1 ,5) + 1.4;
% H_s_table(1 ,6) = H_s_table(1 ,6) + 1.3;