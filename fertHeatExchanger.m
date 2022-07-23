function [bhxVolume,bhxArea] = fertHeatExchanger(numTube)

%function used to define the fertile heat exchanger properties based off of
%the inputs of number of tubes (once-through heat exchanger)

%Volume for the Fertile Salt (Blanket), Coolant Salt (Secondary), and tube
%metal are calculated and provided. Values of Tube Thickness, outer
%diameter, and lengths are kept constant. Only tube number is varied.

%Purpose of this function is to be used with Magic Systems of Systems
%Analyst to perform Parametric Analysis. Magic Systems of Systems Analyst
%will be run over-and-over to provide new tube numbers to this function,
%use this function to determine heat exchanger volume and area properties,
%and then provide the input to the function that runs the entire plant
%simulation.

%% Fertile Heat Exchanger Fertile Salt (Blanket) Volume and Area
OD_b = 0.009525; % (m), 0.375 in Outer diameter of tubes
t_b = 0.000889; % (m), 0.035 in tube thickness
l_b = 2.5146; % (m), 8.25 ft tube length
n_b = numTube; % total number of tubes
% Volume of fuel salt in the tubes
V_fertb = (pi()*((OD_b-2*t_b)^2)/4)*(n_b*l_b); %m^3

A_pblankhx = pi()*(OD_b-2*t_b)*n_b*l_b; % m^2, area of heat transfer for
% fertile salt (primary, inner side of fertile heat exchanger tubes)

%% Fertile Heat Exchanger Tube Cross Sectional Area and Volume Properties

%Cross-sectional Area of Hastelloy N Tubes
A_tbhx = pi()*((OD_b)^2)/4 - pi()*((OD_b-2*t_b)^2)/4; %m^2
%Tube metal Volume, fertile heat exchanger
V_tblankhx = A_tbhx*l_b*n_b; %m^3

%% Fertile Heat Exchanger Secondary Coolant Area and Volume Properties

%A_sbhx = 123.561; %m^2, 1330 ft^2 from Table 3.12 of ORNL-3996
A_sblankhx = pi()*(OD_b)*n_b*l_b; %m^2, Calculated for different tube number

% Shell ID for Fertile HX ORNL-3996 Table 3.12
d_bi = 0.9271; % m, 36.5 in
%total Volume of fertile HX using shell ID 
V_bhx = (pi()*(d_bi^2)/4)*l_b; %m^3

%Volume of Secondary fluid in the fertile heat exchanger
%Volume of secondary fluid is equal to volume of heat exchanger minus the
%fertile salt and metal tube volume.
V_sblankhx = (V_bhx - V_fertb) - V_tblankhx; %m^3

bhxVolume = [V_fertb,V_tblankhx,V_sblankhx];
bhxArea = [A_pblankhx,A_sblankhx];
end