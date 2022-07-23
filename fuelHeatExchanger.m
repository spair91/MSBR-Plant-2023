function [fhxVolume,fhxArea] = fuelHeatExchanger(numTubeDown,numTubeUp)

%function used to define the fuel heat exchanger properties based off of
%the inputs of number of tubes (Down, and Up)

%Volume for the Fuel Salt (Primary), Coolant Salt (Secondary), and tube
%metal are calculated and provided. Values of Tube Thickness, outer
%diameter, and lengths are kept constant. Only tube number is varied.

%Purpose of this function is to be used with Magic Systems of Systems
%Analyst to perform Parametric Analysis. Magic Systems of Systems Analyst
%will be run over-and-over to provide new tube numbers to this function,
%use this function to determine heat exchanger volume and area properties,
%and then provide the input to the function that runs the entire plant
%simulation.

%% Fuel Heat Exchanger Primary Volume, and Area Parameters 
% Table 3.11 of ORNL-3996
OD_f = 0.009525; % (m), 0.375 in Outer diameter of tubes
t_f = 0.000889; % (m), 0.035 in tube thickness
l_fdn = 3.56616; % (m), 11.7 ft downflow tube length
l_fup = 4.17576; % (m), 13.7 ft upflow tube length
n_fdn = numTubeDown; % number of tubes in downflow
n_fup = numTubeUp; % number of tubes in upflow

% Volume of fuel salt in the downflow tubes
V_fueldn = (pi()*((OD_f-2*t_f)^2)/4)*(n_fdn*l_fdn); %m^3
% Volume of fuel salt in the upflow tubes
V_fuelup = (pi()*((OD_f-2*t_f)^2)/4)*(n_fup*l_fup); %m^3

% Fuel Primary Heat Exchanger Area of Heat Transfer
A_fuelhx_dn = pi()*(OD_f-2*t_f)*n_fdn*l_fdn; % m^2,
% fuel (primary, inner side of fuel heat exchanger tubes, downflow section)
A_fuelhx_up = pi()*(OD_f-2*t_f)*n_fup*l_fup; % m^2,
% fuel (primary, inner side of fuel heat exchanger tubes, upflow section)

%% Fuel Heat Exchanger Tube Volume and Area Properties
%Cross-sectional Area of Hastelloy N Tubes
A_tfhx = pi()*((OD_f)^2)/4 - pi()*((OD_f-2*t_f)^2)/4; %m^2
%Tube metal Volume, downflow tubes
V_tfuelhxdn = A_tfhx*l_fdn*n_fdn; %m^3
%Tube metal Volume, upflow tubes
V_tfuelhxup = A_tfhx*l_fup*n_fup; %m^3

%% Fuel Heat Exchanger Secondary Volume and Area Properties
A_sfuelhx_dn = pi()*(OD_f)*n_fdn*l_fdn; %m^2, Calculated for different tube numbers.

A_sfuelhx_up = pi()*(OD_f)*n_fup*l_fup; %m^2, Calculated for different tube numbers.

% Shell ID in Fuel downflow area (annular section) ORNL-3996 Table 3.11
d_i_dn = 1.6891; % m, 66.5 in
% Shell ID in Fuel upflow area (center section) ORNL-3996 Table 3.11
d_i_up = 1.02108; % m, 40.2 in
% Shell thickness
t_shell = 0.0127; %m, 0.5 inch ORNL-3996 Table 3.11
%Volume of annulus in the fuel downflow (annular section)
V_fhxdn = ((pi()*(d_i_dn^2)/4)-(pi()*((d_i_up+2*t_shell)^2)/4))*l_fdn; %m^3
%Volume of center in the fuel upflow (center section) 
V_fhxup = (pi()*(d_i_up^2)/4)*l_fup; %m^3

%Volume of Secondary fluid in the fuel downflow (annular section)
%Equals previously calculated total volume of annulus minus tube volume and
%fuel volume in the downflow section.
V_sfuelhxdn = (V_fhxdn - V_fueldn) - V_tfuelhxdn; %m^3
%Volume of Secondary fluid in the fuel upflow (center section)
%Equals previously calculated total volume of center minus tube volume and
%fuel volume in the upflow section.
V_sfuelhxup = (V_fhxup - V_fuelup) - V_tfuelhxup; %m^3


fhxVolume = [V_fueldn,V_fuelup,V_tfuelhxdn,V_tfuelhxup,V_sfuelhxdn,V_sfuelhxup];
fhxArea = [A_fuelhx_dn,A_fuelhx_up,A_sfuelhx_dn,A_sfuelhx_up];
end