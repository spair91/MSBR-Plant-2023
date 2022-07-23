%% MSBR Plant Parameters File
% Author: Stephen Pair
% Advised by: Dr. Micheal Jones
% 
% MATLAB to support SAM 2022 conference paper:
% Stephen Pair, Micheal Jones. 2022. Model-based systems engineering and 
% simulation of a molten salt reactor power plant for requirements analysis. 
%
% Note that this MATLAB model is for a single-reactor MSBR plant proposed by
% ORNL:
% 
% Paul R. Kasten, E. S. Bettis, and Roy C. Robertson. 1966. Design Studies 
% of 1000-Mw(e) Molten-Salt Breeder Reactors. Oak Ridge National Laboratory
% Document ORNL-3996. Oak Ridge National Laboratory, Oak Ridge, TN.
% 
% The modeling techniques used here is a a modification of the previous 
% modeling work performed by The University of Tennessee and ORNL: 
% 
% V. Singh, A.M. Wheeler, B.R. Upadhyaya, O. Chvála, M.S. Greenwood, 
% Plant-level dynamic modeling of a commercial-scale molten salt reactor 
% system, Nuclear Engineering and Design, Vol. 360, 2020.
% 
% &
%
% Vikram Singh, Matthew R. Lish, Alexander M. Wheeler, Ondřej Chvála, and 
% Belle R. Upadhyaya. 2018. Dynamic modeling and performance analysis of a 
% two-fluid molten-salt breeder reactor system. Nucl. Technol. 202, 1, 
% (Feb, 2018), 177-193. DOI: https://doi.org/10.1080/00295450.2017.1416879.
%
% &
%
% Thomas W. Kerlin. 1967. Preliminary Dynamics Model of the MSBR. 
% Oak Ridge National Laboratory Document ORNL-MSR-67-102. 
% Oak Ridge National Laboratory, Oak Ridge, TN.
%
% Other references for values and methods are provided throughout. Refer to
% references section of the 2022 paper for other references.

%% Core Mass Paramters
%MSBR Core used is the Single-Core MSBR Plant proposed in ORNL-3996

core_vol = 27.8; %(m^3), 982 ft^3 Table 3.1, total core volume

fuel_core_vol_frac = 0.169; % Table 3.1, volume fraction of fuel in core
fert_core_vol_frac = 0.0735; % Table 3.1, volume fraction of fertile salt in core
grap_core_vol_frac = 0.7575; % Table 3.1, volume fraction of graphite in core

core_fuelvol = fuel_core_vol_frac*core_vol; % core fuel volume m^3
core_fertvol = fert_core_vol_frac*core_vol; % core fertile volume m^3
core_grapvol = grap_core_vol_frac*core_vol; % core grphite volume m^3

den_f = 2034; % kg/m^3, 127 lb/ft^3, Table 3.2
den_b = 4437; % kg/m^3, 277 lb/ft^3, Table 3.2
den_g = 2002; % kg/m^3, 125 lb/ft^3, Table 3.2

mass_fuel = core_fuelvol*den_f; % core fuel mass [kg]
mass_fert = core_fertvol*den_b; % core fertile salt mass [kg]
mass_grap = core_grapvol*den_g; % core graphite mass [kg]

mass_fn = mass_fuel/4; % mass of each of the four fuel lumps in core
mass_bn = mass_fert/2; % mass of each of the two fertile lumps in core

%% Core Heat Deposition Fractions
% Note that total Gross Thermal heat is generated in the core and blanket
% as shown on Table 3.1 of ORNL-3996 (2114 + 111) MW. We'll consider 
% the gross thermal to be generated in the core only. Note that the 
% dynamics analysis in ORNL-MSR-67-102 also does this. See page 11 of
% ORNL-MSR-67-102. Compare the fraction of power in the blanket and
% annulus and compare that value with the total power fraction of adding
% the Kf1, Kf2, Kf3, Kf4, Kg1, Kg2, Kg3, Kb1, and Kb2 (from
% ORNL-MSR-67-102) together.
therm_P = 2225; %[MW], ORNL-3996 Table 3.1
therm_P_blank = 111; %[MW], ORNL-3996 Table 3.1
K = (therm_P - therm_P_blank)/therm_P; % fraction of power in core only
K_Ba = (therm_P_blank)/therm_P; % fraction of power in blanket and annulus

K_G = 0.074; % total fraction of power in the graphite (ORNL-MSR-67-102 pg. 10)
K_B = 0.017; % total fraction of heat generated in the fertile salt within the core
% % ORNL-MSR-67-102 pg. 11.
K_fuel = 0.89; %total power in graphite removed by fuel stream (ORNL-MSR-67-102 pg. 5)
% % Note that there are 2 graphite lumps for the core.
k_g1 = ((K_G*K_fuel)/2)/K; % power gen in graphite lump 1 normalized
k_g2 = k_g1; % power generated in graphite lump 2
% % Note that there are 4 fuel lumps for the core
k_f1 = ((1 - K_G - K_B - K_Ba)/4)/K; % power gen in fuel lump 1 normalized
k_f2 = k_f1; % power gen in fuel lump 2 normalized
k_f3 = k_f1; % power gen in fuel lump 3 normalized
k_f4 = k_f1; % power gen in fuel lump 4 normalized

k_b1 = (K_B/2)/K; % power gen in blanket lump 1 normalized, ORNL-MSR-62-107 pg. 16
k_b2 = k_b1; % power gen in blanket lump 2 normalized

k_g3 = (K_G*(1-K_fuel))/K; % power gen in graphite lump 3 normalized

%% Core Total Heat Capacity Values of Lumps
% graphite specific heat, from ORNL-TM-0728 pg. 87
Cp_g = 1.758e-3; % MJ/kg/deg C, 0.42 Btu/lb/deg F
% fuel specific heat capacity, from ORNL-3996 pg. 41
Cp_f = 2.303e-3; % MJ/kg/deg C, 0.55 Btu/lb/deg F
% fertile salt specific heat capacity, from ORNL-3996 pg. 41
Cp_b = 0.921e-3; % MJ/kg/deg C, 0.22 Btu/lb/deg F

mcp_f1 = mass_fn*Cp_f; % total heat capacity of fuel salt lump 1
mcp_f2 = mass_fn*Cp_f; % total heat capacity of fuel salt lump 2
mcp_f3 = mass_fn*Cp_f; % total heat capacity of fuel salt lump 3
mcp_f4 = mass_fn*Cp_f; % total heat capacity of fuel salt lump 4

% Total heat capacity of graphite lump 1, see ORNL-MSR-67-102 pg. 10
mcp_g1 = (mass_grap*Cp_g)*K_fuel/2;
mcp_g2 = mcp_g1; % heat capacity of graphite lump 2

% Total heat capacity of graphite lump 3,
mcp_g3 = (mass_grap*Cp_g)*(1-K_fuel); % heat capacity of graphite 3 ORNL-MSR-67-102 pg. 15

% Fertile salt nodes total heat capacity
mcp_b1 = mass_bn*Cp_b; % heat capacity of fertile node 1, ORNL-MSR-67-102 pg. 16
mcp_b2 = mcp_b1; % heat capacity of fertile node 2

%% Fractions of Heat Between Fuel/Fertile Lumps and Graphite Lumps
% core upflow
k_1 = 0.5; %fraction of heat from the graphite node 1 which goes to the fuel lump 1.
k_2 = 0.5; %fraction of heat from the graphite node 1 which goes to the fuel lump 2.
% core downflow
k_3 = 0.5; %fraction of heat from the graphite node 2 which goes to the fuel lump 3.
k_4 = 0.5; %fraction of heat from the graphite node 2 which goes to the fuel lump 4.
% Fertile flow
k_1b = 0.5; %fraction of heat from the graphite node 3 which goes to the fertile lump 1.
k_2b = 0.5; %fraction of heat from the graphite node 3 which goes to the fertile lump 2.

%% Thermal Power
P = 2225; % MW

%% Heat Transfer Coefficient and Area for Fuel/Graphite Interfaces
h_fg = 8.165e-3; %(MW/m2/deg C); 1438 Btu/hr/ft^2/deg F ORNL-MSR-67-102 pg. 9

n_fuel_tubes = 534; %number of graphite fuel tubes ORNL-3996 Table 3.1
n_passages = 8; %number of fuel upflow passages per fuel tube ORNL-3996 pg. 26
D_fuel_up = 0.0135; % m, 0.53 in diameter of the fuel upflow passages ORNL-3996 pg. 26
L_core = 3.81; %m, 12.5 ft length of active core ORNL-3996 Table 3.1

A_fuel_up = pi()*D_fuel_up*n_fuel_tubes*n_passages*L_core; %m^2; area of fuel upflow heat transfer
hA_fg_up = h_fg*A_fuel_up; %MW/deg heat transfer coefficient of fuel upflow node

D_fuel_down = 0.0381; %m, 1.5 in diameter of the fuel downflow passages ORNL-3996 pg. 26
A_fuel_down = pi()*D_fuel_down*n_fuel_tubes*L_core; %m^2; area of fuel downflow heat transfer area
hA_fg_dn = h_fg*A_fuel_down; %MW/deg C, heat transfer coefficient of fuel downflow node

%% Heat Transfer Coefficient for Fertile /Graphite Interface
% Note that the Area for heat transfer in ORNL-MSR-67-102 is set equal to
% 4000 ft^2 without explanation. Dimensions for the estimation
% of surface area between the fertile/graphite interface are not available
% for the MSBR reactor in ORNL-3996. Consider scaling up the provided
% surface area of 4000 ft^2 for the fert/graphite interface by a factor
% of the A_fuel_down calculated here divided by the area fuel down
% calculated in ORNL-MSR-67-102. Also, consider scaling up by the areas of
% fuel upflow calculated here and in ORNL-67-102.
% 4000 ft^2 = 371.6122 m^2    A_fuel_down/(1319 ft^2) = 1.99
% A_fuel_up/(2285 ft^2) = 3.25. Scale the 4000 ft^2 area between these two
% numbers.

A_fert = (4000*.09290304)*2.62; %m^2, area of fertile/graphite heat exchange
% fertile salt/graphite heat transfer coefficient
h_bg = 2.839e-3; %(MW/m2/deg C); 500 Btu/hr/ft^2/deg F ORNL-MSR-67-102 pg. 15
hA_bg = A_fert*h_bg; %(MW/deg C); heat tranfer*area coefficient for fert/graphite node

%% Fuel and Fertile Core Lump Transit values
% Use the same deriviations as ORNL-MSR-67-102 with the appropriate fuel
% velocity and core height in ORNL-3996
fuel_vel = 4.572;% {m/s}, 15 ft/s fuel salt velocity in core ORNL-3996 pg. 26
tau_1 = (1/2)*L_core/fuel_vel; %s, fuel transit time across fuel lump 1
tau_2 = tau_1; %s, fuel transit time across fuel lump 2 upflow
tau_3 = tau_1; %s, fuel transit time across fuel lump 3 downflow
tau_4 = tau_1; %s, fuel transit time across fuel lump 4 downflow

% Fertile Lump transit Time
% Calculate Blanket lump transit time according to pg. 15/16 of
% ORNL-MSR-67-102
delta_T_fert = 676.67 - 621.11; % deg C, 100 deg F, Temperature differenace
% of fertile stream across the reactor.
W_fert = therm_P_blank/(Cp_b*delta_T_fert); %kg/s, fertile stream mass flow rate
tau_b1 = (0.5)*mass_fert/W_fert; %s, fertile stream transit time for lump b1
tau_b2 = tau_b1; %s, fertile stream transit time for lump b2

%% Fuel/Fertile/Graphite Node Initial Temperatures

% Assume initial temperatures at full power
% see table 3.1 of ORNL-3996 for F4 and B2 temperatures, other
% temperatures are between the fuel/fertile inlet and outlet temperature.
%%% NOTE THAT ACTUAL TEMPERATURE OF THE VARIABLES ARE DIFFERENT FROM deg F TEMPS
%%% THIS IS FROM SEEING STEADY-STATE REACTION OF THE CORE.
T0_f1  = 578; % in �C  1075 deg F approx. in ORNL-3996
T0_f2  = 619; % in �C  1150 deg F approx. in ORNL-3996
T0_f3  = 660; % in �C, 1225 deg F, approx. in ORNL-3996
T0_f4  = 702; %C, 1300 deg F Fuel Outlet Temperature ORNL-3996 Fig. 3.7
% approx. in ORNL-3996
T0_b1  = 638; % in �C, 1200 deg F approx. in ORNL-3996
%Fertile Salt Input Temperature (fertile salt from the core)
T0_b2 = 653; %C, 1250 deg F ORNL-3996 Figure 3.7 approx. in ORNL-3996

% Overall, graphite temps are between temps of lower and upper fuel and 
% blanket nodes. 
% Use initial temperatures of 19 deg F larger than initial fuel node
T0_g1  = 591; % in �C, 
T0_g2  = 700; % in �C, 

%Fertile Salt Graphite lump
T0_g3  = 645; % in �C, 

% Inlet temp Tf_in 1000 degF, approx. see ORNL-3996 Table 3.1
Tf_in  = 536; % in �C

% Inlet temp T_b_in 1150 degF, approx.  see ORNL-3996 Table 3.1
Tb_in = 623; % in �C

%% Core Neutronics Parameters
W_f = 5508.627; % kg/s, 43720000 lb/hr fuel salt flow rate ORNL-3996 Table 3.1

tau_c = mass_fuel/W_f; %s Fuel core transit time

%% Calculate fuel heat exchanger volumes first and then calculated fuel transit time
Wp_fhx1 = 1360.777; %kg/s, 1.08e7 lb/hr mass flow rate in Fuel Hx 1
Wp_fhx2 = Wp_fhx1; %kg/s, 1.08e7 lb/hr mass flow rate in Fuel Hx 2
Wp_fhx3 = Wp_fhx1; %kg/s, 1.08e7 lb/hr mass flow rate in Fuel Hx 3
Wp_fhx4 = Wp_fhx1; %kg/s, 1.08e7 lb/hr mass flow rate in Fuel Hx 4
W_f1 = Wp_fhx1; %kg/s, Flow out of fuel primary HX 1 and into Pre-Core Mixing
W_f2 = Wp_fhx2; %kg/s, Flow out of fuel primary HX 2 and into Pre-Core Mixing
W_f3 = Wp_fhx3; %kg/s, Flow out of fuel primary HX 3 and into Pre-Core Mixing
W_f4 = Wp_fhx4; %kg/s, Flow out of fuel primary HX 4 and into Pre-Core Mixing
OD_f = 0.009525; % (m), 0.375 in Outer diameter of tubes
t_f = 0.000889; % (m), 0.035 in tube thickness
l_fdn = 3.56616; % (m), 11.7 ft downflow tube length
l_fup = 4.17576; % (m), 13.7 ft upflow tube length
%%Original Tube Numbers
n_fdn = 4167; % number of tubes in downflow
n_fup = 3624; % number of tubes in upflow
% Volume of fuel salt in the downflow tubes
V_fdn = (pi()*((OD_f-2*t_f)^2)/4)*(n_fdn*l_fdn); %m^3
% Volume of fuel salt in the upflow tubes
V_fup = (pi()*((OD_f-2*t_f)^2)/4)*(n_fup*l_fup); %m^3
%Fuel heat exchanger mixing plenum transit time.
tau_mp = 0.5; % s, ORNL-MSR-67-102 pg. 21 assumes 1 s. We'll assume 0.5 s
%Volume of the mixing plenum for each HX.
V_mp = tau_mp*Wp_fhx1/den_f;
%Total Volume of Fuel in Heat Exchangers
V_hx = 4*(V_fdn+V_fup+V_mp);

%% Back to Core Neutronics

pipe_fv = 4.163 + 2.77611 + V_hx; % m^3, (147 + 345) ft^3 volume of fuel salt in plena, heat 
% exchangers, and piping, ORNL-3996 Table 3.1 (when using the number of
% tubes of original heat exchangers for fuel)
pipe_fm = pipe_fv*den_f; % kg mass of fuel in external pipes and plena
tau_l = pipe_fm/W_f; %s Fuel loop transit time (loop being outside of the core)
Lam = 3.3e-4; %s, mean neutron generation time, ORNL-MSR-67-102 pg. 5
% lam and beta values from ORNL-MSR-67-102 pg 6, beta is Betai (MSBR)
lam = [1.260e-02, 3.370e-02, 1.390e-01, 3.250e-01, 1.130e+00, 2.500e+00]; 
beta = [2.290e-04, 8.320e-04, 7.100e-04, 8.520e-04, 1.710e-04, 1.020e-04];
beta_t = sum(beta); % total delayed neutron fraction MSBR
% rho_0 equation from ORNL-MSR-67-102 pg 4, reacitvity required for
% steady-state
rho_0 = beta_t - bigterm(beta,lam,tau_l,tau_c); % reactivity change in
% going from stationary to circulating fuel

n_frac0 = 1; % initial fractional neutron density (i.e. full power)

% initial fractional delayed-neutron densities for groups 1..6
C0(1)   = ((beta(1))/Lam)*(1.0/(lam(1) - (exp(-lam(1)*tau_l) - 1.0)/tau_c));
C0(2)   = ((beta(2))/Lam)*(1.0/(lam(2) - (exp(-lam(2)*tau_l) - 1.0)/tau_c));
C0(3)   = ((beta(3))/Lam)*(1.0/(lam(3) - (exp(-lam(3)*tau_l) - 1.0)/tau_c));
C0(4)   = ((beta(4))/Lam)*(1.0/(lam(4) - (exp(-lam(4)*tau_l) - 1.0)/tau_c));
C0(5)   = ((beta(5))/Lam)*(1.0/(lam(5) - (exp(-lam(5)*tau_l) - 1.0)/tau_c));
C0(6)   = ((beta(6))/Lam)*(1.0/(lam(6) - (exp(-lam(6)*tau_l) - 1.0)/tau_c));

%% Temperature Feedback Coefficients
a_f = -8.172e-5; % fuel temperature feedback coefficient in drho/ deg C
% -4.54e-5 drho/ deg F from ORNL-MSR-67-102
a_g = 2.016e-5; % fuel temperature feedback coefficient in drho/ deg C
% 1.12e-5 drho/ deg F from ORNL-MSR-67-102
a_b = 1.656e-5; % fuel temperature feedback coefficient in drho/ deg C
% 9.2e-6 drho/ deg F from ORNL-MSR-67-102

%% Xenon Reactivity Parameters
% Yield of I-135 and Xe-135 for U-233 fission (pulled from Univ of 
% Tennessee MSDR Model). 
% https://www-nds.iaea.org/wimsd/fpyield.htm
% IAEA website helps corroborate these numbers for U-233 yield
gamma_I  = 5.1135e-2; % weighted yield for I-135
gamma_Xe = 1.1628e-2; % weighted yield for Xe-135

% Reference to use for decay constant of Iodine-135, ORNL-4037 pg. 164
lam_I    = 2.875e-5;  % decay constant for I-135 (s^-1)
% Reference to use for decay constant of Xenon, ORNL-4191 pg. 91.
lam_Xe   = 2.0916e-5; % decay constant for Xe-135 (s^-1)

% Ridley, Chvala 2017 paper pg. 272
lam_bubl = 2.0e-2; % effective bubbling out constant (s^-1)

% Microscopic cross-section of Xe-135
% Used value in the MSDR 2019 Powerplant model by Vikram et. al.
sig_Xe   = 2.66449e-18; % (cm^2) microscopic cross-section for Xe (n,gamma) reaction
% Molecular weight, see ORNL-3996 Table 3.2. Used molecular percentages
% of U-233 and U-235 in ORNL-MSR-67-102 pg. 6
molc_wt  = .683*(7.016+18.998)+.312*(9.012+2*18.998)+.005*(4*18.998+(233.0396*0.934 + 235.0439*0.066)); % (g/mol)
molc_den = 0.001*den_f/molc_wt;  % (mol/cm^3)

% Uranium density per cm^3, 0.5% mole % UF4 per fuel molecule
U_den    = .005*molc_den*6.022e23;  % (#U/cm^3)

% Uranium-233 thermal neutron fission cross section
% Used BNL National Nuclear Data Center value
% https://www.nndc.bnl.gov/atlas/atlasvalues.html
% Note that ORNL-MSR-67-102 has the fuel as mostly U-233.
U_sig    = 5.291e-22; % (cm^2)

Sigma_f_msbr = U_den*U_sig; % (cm^-1)
% Steady state thermal neutron flux
% Ridley/Chvala 2017 paper pg. 269
% 3.04414e-17*1e6 is 190 MeV converted.
% 190 MeV is energy released per fission of U-233
% See Nuclear Reactor Analysis by Duderstadt and Hamilton pg.  67
phi_0 = P/(3.04414e-17*1e6*Sigma_f_msbr);  % neutrons cm^-2 s^-1

% initial concentration of iodine, see 
I_0 = gamma_I*Sigma_f_msbr*phi_0/lam_I;

% initial concentration of Xenon,
% note that lam_bubl is the removal of Xenon through bubbling/sparging
Xe_0 = (gamma_Xe*Sigma_f_msbr*phi_0 + lam_I*I_0)/(lam_Xe + sig_Xe*phi_0 + lam_bubl);

Xe0_og = lam_bubl*Xe_0/(lam_Xe); % initial Xe conc. in off-gas system

% assumed absorption cross section of the MSDR power plant core.
% See MSDR Power Plant 2019 Pape by Vikram et. al.
Sig_a = 1.02345; % (cm^-1) macroscopic absorption cross-section for core

% steady-state reactivity due to Xe-135, See MSDR Power Plant paper 2019
% by Vikram et. al.
rhoXe_0 = (-sig_Xe/Sig_a)*(gamma_I+gamma_Xe)*Sigma_f_msbr*phi_0/(lam_Xe + sig_Xe*phi_0 + lam_bubl);

%% Fuel Primary Heat Exchanger Flow, Mass, Heat Capacity, Heat Transfer Parameters 

% Mass of fuel salt in HX downflow tubes
M_tdn = den_f*V_fdn; % kg
% Mass of fuel salt in HX upflow tubes
M_tup = den_f*V_fup; % kg
tau_tdn = M_tdn/Wp_fhx1; % s Transit time in downflow tubes of Fuel Hx
tau_tup = M_tup/Wp_fhx1; % s Transit time in upflow tubes of Fuel Hx

mp_fhx1m = tau_mp*Wp_fhx1; % kg, mass of fluid in mixing plenum of Fuel HX
mp_fhx2m = tau_mp*Wp_fhx2; % kg, mass of fluid in mixing plenum of Fuel HX
mp_fhx3m = tau_mp*Wp_fhx3; % kg, mass of fluid in mixing plenum of Fuel HX
mp_fhx4m = tau_mp*Wp_fhx4; % kg, mass of fluid in mixing plenum of Fuel HX

% Fuel Primary Heat Exchanger fuel mass lumps.
mp_fhx1 = [M_tdn/2,M_tdn/2,M_tup/2,M_tup/2];
mp_fhx2 = mp_fhx1;
mp_fhx3 = mp_fhx1;
mp_fhx4 = mp_fhx1;

% Fuel Primary Heat Exchanger Total Heat Capacity of Lumps
mcp_p_fhx1 = mp_fhx1*Cp_f;
mcp_p_fhx2 = mp_fhx2*Cp_f;
mcp_p_fhx3 = mp_fhx3*Cp_f;
mcp_p_fhx4 = mp_fhx4*Cp_f;

% Fuel Primary Heat Exchanger Area of Heat Transfer
A_pfhx_dn = pi()*(OD_f-2*t_f)*n_fdn*l_fdn; % m^2,
% fuel (primary, inner side of fuel heat exchanger tubes, downflow section)
A_pfhx_up = pi()*(OD_f-2*t_f)*n_fup*l_fup; % m^2,
% fuel (primary, inner side of fuel heat exchanger tubes, upflow section)

% heat transfer coefficients from ORNL-MSR-67-102 pg. 21
% heat transfer coefficient from primary to tubes
h_p = 1.306e-02; % Mw/m^2C, 2300 Btu/hr/ft^2/deg F

%Mw/deg C, h*A of primary side of fuel heat exchanger
hA_p_fhx1 = [h_p*A_pfhx_dn/2,h_p*A_pfhx_dn/2,h_p*A_pfhx_up/2,h_p*A_pfhx_up/2];
hA_p_fhx2 = hA_p_fhx1;
hA_p_fhx3 = hA_p_fhx1;
hA_p_fhx4 = hA_p_fhx1;

%% Fuel Heat Exchanger Tube Lump Mass Properties
%Cross-sectional Area of Hastelloy N Tubes
A_tfhx = pi()*((OD_f)^2)/4 - pi()*((OD_f-2*t_f)^2)/4; %m^2
%Tube metal Volume, downflow tubes
V_tfhxdn = A_tfhx*l_fdn*n_fdn; %m^3
%Tube metal Volume, upflow tubes
V_tfhxup = A_tfhx*l_fup*n_fup; %m^3

%Hastelloy N tube density from ORNL-TM-0728 pg. 20
den_t = 8774.528; %kg/m^3, 0.317 lb/in^3
%Hastelloy N tube specific heat capacity, ORNL-TM-0728 pg. 20
Cp_t  =  577.78e-6; %Mw*s/kg/C, 0.138 Btu/lb/deg F
%Mass*Cp of the tube nodes for the fuel heat exchanger
mcp_t_fhx1 = [den_t*V_tfhxdn*Cp_t, den_t*V_tfhxup*Cp_t];
mcp_t_fhx2 = mcp_t_fhx1;
mcp_t_fhx3 = mcp_t_fhx1;
mcp_t_fhx4 = mcp_t_fhx1;

%% Fuel Heat Exchanger Secondary Coolant Flow, Mass, Heat Capacity, Heat Transfer Parameters

%Fuel heat exchanger coolant flow rate (ORNL-3996 Table 3.11)
Ws_fhx1 = 2116.764; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fuel Hx 1
Ws_fhx2 = Ws_fhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fuel Hx 2
Ws_fhx3 = Ws_fhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fuel Hx 3
Ws_fhx4 = Ws_fhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fuel Hx 4

% heat transfer coefficients from ORNL-MSR-67-102 pg. 21
% heat transfer coefficient from secondary to tubes
h_s = 2.555e-02; % Mw/m^2C, 4500 Btu/hr/ft^2/deg F

%A_sfhx_dn = 445.01; %m^2, 4790 ft^2 from Table 3.11 of ORNL-3996
A_sfhx_dn = pi()*(OD_f)*n_fdn*l_fdn; %m^2, Calculated for different tube numbers.
% Note that annular section heat transfer area is used for secondary
% downflow area. The term downflow is with respect to the fuel salt flow in
% annular region.
%A_sfhx_up = 452.90; %m^2, 4875 ft^2 from Table 3.11 of ORNL-3996
A_sfhx_up = pi()*(OD_f)*n_fup*l_fup; %m^2, Calculated for different tube numbers.
% Note that center section heat transfer area is used for secondary
% upflow area. The term upflow is with respect to the fuel salt flow in
% center region.

%Mw/deg C, h*A of secondary side of fuel heat exchanger
% Note that upflow is used for the first two secondary nodes and downflow
% is used for the last two. This is because Lump 1 and 2 of the secondary
% is the secondary fluid transferring heat with with the tube and fuel
% upflow lumps. See ORNL-MSR-67-102 pg. 18
hA_s_fhx1 = [h_s*A_sfhx_up/2,h_s*A_sfhx_up/2,h_s*A_sfhx_dn/2,h_s*A_sfhx_dn/2];
hA_s_fhx2 = hA_s_fhx1;
hA_s_fhx3 = hA_s_fhx1;
hA_s_fhx4 = hA_s_fhx1;

%Calculate Secondary lumped masses by finding the center and annulus
%volumes of the heat exchanger and subtract the tube metal volume and fuel
%volume. Assume the length of the center and annulus volumes equal the 
% associated tube lengths. Baffles, disks, donuts, and tubesheets are not 
% considered as they have minimal effect on overall volume.

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
V_sfhxdn = (V_fhxdn - V_fdn) - V_tfhxdn; %m^3
%Volume of Secondary fluid in the fuel upflow (center section)
%Equals previously calculated total volume of center minus tube volume and
%fuel volume in the upflow section.
V_sfhxup = (V_fhxup - V_fup) - V_tfhxup; %m^3

%density of the coolant salt ORNL-MSR-67-102 pg. 23
den_s = 2002.308; %kg/m^3, 125 lb/ft^3

%mass of secondary fluid in the fuel downflow (annular) section
m_sdn = V_sfhxdn*den_s; %kg
%mass of secondary fluid in the fuel upflow (center) section
m_sup = V_sfhxup*den_s; %kg

%Fuel heat exchanger coolant mass lumps. Assumed same mass for each lump 1
%and 2, and same mass for each lump 3 and.
ms_fhx1 = [m_sup/2,m_sup/2,m_sdn/2,m_sdn/2]; %kg
ms_fhx2 = ms_fhx1;
ms_fhx3 = ms_fhx1;
ms_fhx4 = ms_fhx1;
%coolant salt specific heat, ORNL-3996 Table 3.2
Cp_sec = 1.716588e-3; %Mw*s/kg/deg C, 0.41 Btu/lb*deg F
%coolant salt m*Cp values for fuel heat exchanger
mcp_s_fhx1 = ms_fhx1*Cp_sec; %Mw*s/deg C
mcp_s_fhx2 = ms_fhx2*Cp_sec; %Mw*s/deg C
mcp_s_fhx3 = ms_fhx3*Cp_sec; %Mw*s/deg C
mcp_s_fhx4 = ms_fhx4*Cp_sec; %Mw*s/deg C

%% Fertile Heat Exchanger Flow, Mass, Heat Capacity, Heat Transfer Table 3.12 of ORNL-3996
Wp_bhx1 = 541.7909; %kg/s, 4.3e6 lb/hr mass flow rate in blanket Hx 1
Wp_bhx2 = Wp_bhx1; %kg/s, 4.3e6 lb/hr mass flow rate in blanket Hx 2
Wp_bhx3 = Wp_bhx1; %kg/s, 4.3e6 lb/hr mass flow rate in blanket Hx 3
Wp_bhx4 = Wp_bhx1; %kg/s, 4.3e6 lb/hr mass flow rate in blanket Hx 4
W_b1 = Wp_bhx1; %kg/s, Flow out of blanket HX 1 and into Pre-Core Mixing
W_b2 = Wp_bhx2; %kg/s, Flow out of blanket HX 2 and into Pre-Core Mixing
W_b3 = Wp_bhx3; %kg/s, Flow out of blanket HX 3 and into Pre-Core Mixing
W_b4 = Wp_bhx4; %kg/s, Flow out of blanket HX 4 and into Pre-Core Mixing
OD_b = 0.009525; % (m), 0.375 in Outer diameter of tubes
t_b = 0.000889; % (m), 0.035 in tube thickness
l_b = 2.5146; % (m), 8.25 ft tube length
%original tube number
n_b = 1641; % total number of tubes
%n_b = 975; % total number of tubes
% Volume of fuel salt in the tubes
V_fb = (pi()*((OD_b-2*t_b)^2)/4)*(n_b*l_b); %m^3
% Total Mass of blanket salt in HX tubes
M_tb = den_b*V_fb; % kg
tau_bhx = M_tb/Wp_bhx1; % s Fertile Transit time in tubes of blanket HX
% no mixing plenum transit time for blanket HX

% Fertile Primary Heat Exchanger fuel mass lumps.
mp_bhx1 = [M_tb/2,M_tb/2];
mp_bhx2 = mp_bhx1;
mp_bhx3 = mp_bhx1;
mp_bhx4 = mp_bhx1;

%Fertile Primary Heat Exchanger M*Cp values for lumps
mcp_p_bhx1 = mp_bhx1*Cp_b;
mcp_p_bhx2 = mp_bhx2*Cp_b;
mcp_p_bhx3 = mp_bhx3*Cp_b;
mcp_p_bhx4 = mp_bhx4*Cp_b;

A_pbhx = pi()*(OD_b-2*t_b)*n_b*l_b; % m^2, area of heat transfer for
% fertile salt (primary, inner side of fertile heat exchanger tubes)

% heat transfer coefficients from ORNL-MSR-67-102 pg. 26
% Note that ORNL-MSR-67-102 uses the same h on both sides by an assumption.
% It may be worthwhile to explore ORNL-1545 for other heat transfer
% coefficient values.
% heat transfer coefficient from fertile primary to tubes
h_pb = 1.1584e-02; % Mw/m^2C, 2040 Btu/hr/ft^2/deg F

%Mw/deg C, h*A of primary side of fertile heat exchanger
hA_p_bhx1 = [h_pb*A_pbhx/2,h_pb*A_pbhx/2];
hA_p_bhx2 = hA_p_bhx1;
hA_p_bhx3 = hA_p_bhx1;
hA_p_bhx4 = hA_p_bhx1;

%% Fertile Heat Exchanger Tube Lump Mass Properties

%Cross-sectional Area of Hastelloy N Tubes
A_tbhx = pi()*((OD_b)^2)/4 - pi()*((OD_b-2*t_b)^2)/4; %m^2
%Tube metal Volume, fertile heat exchanger
V_tbhx = A_tbhx*l_b*n_b; %m^3

%Mass*Cp of the tube nodes for the fertile heat exchanger
mcp_t_bhx1 = den_t*V_tbhx*Cp_t;
mcp_t_bhx2 = mcp_t_bhx1;
mcp_t_bhx3 = mcp_t_bhx1;
mcp_t_bhx4 = mcp_t_bhx1;

%% Fertile Heat Exchanger Secondary Coolant Flow, Mass, Heat Capacity, Heat Transfer Parameters

%Fertile primary heat exchanger coolant flow rate (ORNL-3996 Table 3.12)
Ws_bhx1 = 2116.764; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fertile Hx 1
Ws_bhx2 = Ws_bhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fertile Hx 2
Ws_bhx3 = Ws_bhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fertile Hx 3
Ws_bhx4 = Ws_bhx1; %kg/s, 1.68e7 lb/hr mass flow rate coolant in Fertile Hx 4
W_s1 = Ws_bhx1; %kg/s, secondary Flow out of blanket HX 1 and into Pre-BoP Mixing
W_s2 = Ws_bhx2; %kg/s, secondary Flow out of blanket HX 2 and into Pre-BoP Mixing
W_s3 = Ws_bhx3; %kg/s, secondary Flow out of blanket HX 3 and into Pre-BoP Mixing
W_s4 = Ws_bhx4; %kg/s, secondary Flow out of blanket HX 4 and into Pre-BoP Mixing

%A_sbhx = 123.561; %m^2, 1330 ft^2 from Table 3.12 of ORNL-3996
A_sbhx = pi()*(OD_b)*n_b*l_b; %m^2, Calculated for different tube number

% heat transfer coefficients from ORNL-MSR-67-102 pg. 26
% Note that ORNL-MSR-67-102 uses the same h on both sides by an assumption.
% It may be worthwhile to explore ORNL-1545 for other heat transfer
% coefficient values.
% heat transfer coefficient from secondary to tubes in fertile primary HX
h_sb = 1.1584e-02; % Mw/m^2C, 2040 Btu/hr/ft^2/deg F

%Mw/deg C, h*A of secondary side of fertile heat exchanger
hA_s_bhx1 = [h_sb*A_sbhx/2,h_sb*A_sbhx/2];
hA_s_bhx2 = hA_s_bhx1;
hA_s_bhx3 = hA_s_bhx1;
hA_s_bhx4 = hA_s_bhx1;

%Calculate Secondary lumped masses by finding the volume of the heat 
%exchanger and subtract the tube metal volume and fertile salt
%volume. Assume the length of the heat exchanger volume equals the
%tube length. Baffles, disks, donuts, and tubesheets are not considered
%as they have minimal effect on overall volume.

% Shell ID for Fertile HX ORNL-3996 Table 3.12
d_bi = 0.9271; % m, 36.5 in
%total Volume of fertile HX using shell ID 
V_bhx = (pi()*(d_bi^2)/4)*l_b; %m^3

%Volume of Secondary fluid in the fertile heat exchanger
%Volume of secondary fluid is equal to volume of heat exchanger minus the
%fertile salt and metal tube volume.
V_sbhx = (V_bhx - V_fb) - V_tbhx; %m^3

%mass of secondary fluid in the fertile heat exchanger
m_sb = V_sbhx*den_s; %kg

%fertile primary heat exchanger secondary coolant mass lumps. 
% Assumed same mass for each lump
ms_bhx1 = [m_sb/2,m_sb/2]; %kg
ms_bhx2 = ms_bhx1;
ms_bhx3 = ms_bhx1;
ms_bhx4 = ms_bhx1;

%coolant salt m*Cp values for fertile heat exchanger
mcp_s_bhx1 = ms_bhx1*Cp_sec; %Mw*s/deg C
mcp_s_bhx2 = ms_bhx2*Cp_sec; %Mw*s/deg C
mcp_s_bhx3 = ms_bhx3*Cp_sec; %Mw*s/deg C
mcp_s_bhx4 = ms_bhx4*Cp_sec; %Mw*s/deg C

%% Fuel/Fertile/Coolant Pre-core/Pre-BoP mixing masses and times
%Set the mixing section of the Pre-Core and Pre-BoP
time_mixf = 1; %s say 1 seconds of time is spent for fuel salt in the mixing
time_mixb = 1; %s say 1 seconds of time is spent for fertile salt in the mixing
time_mixs = 1; %s say 1 seconds of time is spent for coolant salt in the mixing
m_fm = (W_f1+W_f2+W_f3+W_f4)*time_mixf; %kg, mass of fuel in mixing node
m_bm = (W_b1+W_b2+W_b3+W_b4)*time_mixb; %kg, mass of blanket in mixing node
m_sm = (W_s1+W_s2+W_s3+W_s4)*time_mixs; %kg, mass of coolant in mixing node

%% Component-to-Component Transit Times
% With the calculated transition times in fuel HX and fertile HX,
% use the transition times for fuel in external loop and use volumes
% of the fuel or fertile salt outside the core to calculate transition
% times between the core and heat exchanger.
tau_c_fhx = (tau_l - (tau_tdn+tau_tup+tau_mp))/2; %s transit time from core to 
% fuel heat exchanger. Note that because tau_l includes the plena, the 
% tau_c_fhx includes the time across the plena + the piping. Also divide
% the transit time by 2 since the return and supply piping to core is same
% length. This is an assumption as little information is provided on the
% connecting piping to the core.
pipe_bv = 2.053685 + 4*V_fb; %m^3, 100 ft^3 blanket volume in pipe+blanket hx ORNL-3996 Table 3.1 (when using the standard tube number)
pipe_bm = pipe_bv*den_b; %kg, mass of blanket in pipe+blanket hx
tau_lb = pipe_bm/W_fert; %s, blanket external loop transit time
tau_c_fehx = (tau_lb - tau_bhx)/2; %s, blanket transit time from core to 
% fertile heat exchanger

% Time for return of fuel from fuel HX to core
% with the mixing node time removed. Fuel mix node takes some time
% and is accounted for in the Pre-Core mixing node.
time_fhx_c = (tau_l - (tau_tdn+tau_tup+tau_mp))/2 - time_mixf; %s
% set value of taus for each fuel heat exchanger loop
tau_fhx_c = [time_fhx_c, time_fhx_c, time_fhx_c, time_fhx_c];

% Time for return of fertile salt from fertile HX to core
% with the mixing node time removed. Fertile mix node takes some time
% and is accounted for in the Pre-Core mixing node.
time_fehx_c = tau_c_fehx - time_mixb; %s,
% set value of taus for each Fertile heat exchanger loop
tau_fehx_c = [time_fehx_c,time_fehx_c,time_fehx_c,time_fehx_c]; %s

% Coolant Transition time between heat exchanger and OTSG components.
%Assume 1 sec transition from Fertile Heat Exchanger to Boiler/OTSG and
% Boiler/OTSG to Fuel Heat Exchanger.
%Similar to the MSDR Plant provided in Vikram et. al 2019 paper.
% Note that the transition between boiler and fertile heat exchanger
% in ORNL-MSR-67-102 is 13.5 s and the transition between boiler
%and fuel heat exchanger is 4.2 s.
tau_fehx_b = [1, 1, 1, 1];
tau_b_fhx = [1, 1, 1, 1];

%Coolant transition time between heat exchangers.
%Assume 1 sec transition from Fuel Heat Exchanger to Fertile Heat Exchanger
%Similar to the MSDR Plant provided in Vikram et. al 2019 paper.
% Note that the transition between Fuel and Fertile heat exchanger
% in ORNL-MSR-67-102 is 7.88 s
tau_fhx_fehx = [1, 1, 1, 1];

%% Fuel/Fertile Pumps Flow Fractions and Corresponding Heat Transfer (h*A) Fractions
% % Simple lookup table that takes in a Fuel Flow Fraction
% % and finds a Heat Transfer Fraction 
% % See Vikram's Thesis section 3.4.3.1, and ORNL-3014.
% % Note that primary info on Fuel Flow Fraction 
% % and heat transfer coefficient relationship is in personal
% % correspondence. Used digitized Heat Transfer Coefficient lookup table
% % from Vikram et. al 2019 MSDR Plant Paper
xdat=[0.00483288502733
0.024472562101481
0.044096021440641
0.058773071607534
0.083294286913988
0.102852875313184
0.127341655149656
0.146900243548852
0.161577293715744
0.18113588211494
0.19578049681185
0.215339085211046
0.229999917642947
0.254472479744428
0.274031068143625
0.293589656542821
0.313132027207026
0.337588371573515
0.371824010139604
0.406075866440683
0.44029528727178
0.47456336130785
0.513745409046206
0.543083291645
0.58717931308564
0.626328925354014
0.660629434860066
0.694897508896137
0.734079556634493
0.773229168902867
0.822239164045793
0.856507238081864
0.900587041787514
0.949564601460459
0.983865110966511
1.02306337643986
1.05736388594591
1.08675042174968
1.12598112269301
1.16029784993405
1.21422181877926
1.26326424939217
1.3024949503355
1.3466558427161
1.39571449106401
1.43497762747732
1.50863857980663];

ydat=[0.013508209475996
0.024765050705993
0.038273260181989
0.054032837903984
0.072043783871979
0.094557466331972
0.117071148791966
0.139584831251959
0.155344408973954
0.177858091433948
0.198120405647942
0.220634088107936
0.23864503407593
0.263410084781923
0.285923767241916
0.308437449701909
0.333202500407902
0.360218919359894
0.398492179541882
0.434514071477873
0.47503869990586
0.50880922359585
0.544831115531841
0.578601639221831
0.616874899403819
0.657399527831807
0.686667315029798
0.720437838719788
0.756459730655778
0.796984359083766
0.837508987511754
0.871279511201745
0.911804139629733
0.95683150454972
0.98609929174771
1.0198698154377
1.04913760263569
1.07615402158768
1.10542180878568
1.13243822773767
1.17521422441165
1.21123611634765
1.24050390354564
1.26977169074363
1.30354221443362
1.32830726513961
1.3688318935676];

%% Secondary Coolant Heat Transfer Coefficient Lookup Table
% Heat Exchanger Test module used a different Heat Transfer Coefficient
% Lookup Table that was sourced from Vikram Singh's University of Tennessee
% Thesis section 3.4.3.1 and ORNL-3014. The digitized table from V_Singh
% et. al 2017 MSRE Operational Anomaly paper was used. Note that the
% Coolant Salt in the MSRE was LiF-BeF_2. The MSBR uses Na-F-NaBF_4. The
% MSDR 2019 Plant simulation uses the same coefficient lookup tablesfor the
% fuel, secondary, and tertiary salt in the plant. Since the information is 
% lacking for the MSBR coolant salt, it appears to be a good compromise to
% use the same xdat and ydat information. Therefore, set the coolant salt
% heat transfer coefficient lookup table to be the same.

xdatc = xdat;
ydatc = ydat;

%% Initial Temperatures for Heat Exchanger Lumps and Input Temperatures

% Note that these temperatures are selected after viewing the steady-state
% behavior of the reactor plant and heat exchangers. Due to modeling
% choices, the simulation does not match the plant in ORNL-3996 exactly and
% if coolant salt temperatures are sufficiently low enough this causes an
% overcooling of the reactor and a significant steady-state power
% production above 100%.

%Fuel Salt Initial Heat Exchanger Lump Temperatures
%To make sure that the initial temperature leaving the Fuel HX is 1000 deg F
%initialize the fuel lump temperatures in the fuel heat exchanger
%with a temp difference of 75 deg F between each lump.
%T0_p_fhx1 = [T0_f4-41.667, T0_f4-2*41.667, T0_f4-3*41.667, T0_f4-4*41.667];
T0_p_fhx1 = [658, 615, 575, 536];
T0_p_fhx2 = T0_p_fhx1;
T0_p_fhx3 = T0_p_fhx1;
T0_p_fhx4 = T0_p_fhx1;

%Initial Temp of the Fuel Primary Heat Exchanger Plenum Node:
T0_p_fhx1m = T0_p_fhx1(2);
T0_p_fhx2m = T0_p_fhx2(2);
T0_p_fhx3m = T0_p_fhx3(2);
T0_p_fhx4m = T0_p_fhx4(2);

%Fertile Salt Initial Heat Exchanger Lump Temperatures
%To make sure that the initial temperature leaving the Fertile HX is 1150 
%deg F initialize the fertile temperatures with a temp difference of 50 
%deg F between each lump.
%T0_p_bhx1 = [T0_b2-27.78, T0_b2-2*27.78];
T0_p_bhx1 = [638, 623];
T0_p_bhx2 = T0_p_bhx1;
T0_p_bhx3 = T0_p_bhx1;
T0_p_bhx4 = T0_p_bhx1;

%Coolant Salt Initial Fuel Heat Exchanger Lump Temperatures
%To make sure that the coolant initial temperature leaving the Fuel HX is 
%1111 deg F initialize the secondary temperatures with a temp difference of 
%65.25 deg F between each lump.
%T0_s_fhx1 = [T_p6+36.25, T_p6+2*36.25, T_p6+3*36.25, T_p6+4*36.25];
T0_s_fhx1 = [503, 537, 575, 613];
T0_s_fhx2 = T0_s_fhx1;
T0_s_fhx3 = T0_s_fhx1;
T0_s_fhx4 = T0_s_fhx1;
T0_sbhx_in = T0_s_fhx1(4);

%Fuel Heat Exchanger Tube Lumps Initial Temperatures
%Initial Temps for secondary and fuel salt are for full power and for heat
%exchange to occur from fuel to the secondary salt. Therefore, initialize
%tube temperatures to be a temperature between the lowest temperature of
%the fuel lumps it's associated with and the temperature of the highest
%secondary lump it's associated with (e.g. Tube Lump 2 initial temp should
%be between the primary lump 4 initial temp of 1000 deg F and the secondary
%lump 2 initial temp of approx. 980.5 deg F)
%T0_t_fhx1 = [(T0_p_fhx1(2) + T0_s_fhx1(4))/2 , (T0_p_fhx1(4) + T0_s_fhx1(2))/2];
T0_t_fhx1 = [599, 524];
%T0_t_fhx2 = [(T0_p_fhx2(2) + T0_s_fhx2(4))/2 , (T0_p_fhx2(4) + T0_s_fhx2(2))/2];
T0_t_fhx2 = [599, 524];
%T0_t_fhx3 = [(T0_p_fhx3(2) + T0_s_fhx3(4))/2 , (T0_p_fhx3(4) + T0_s_fhx3(2))/2];
T0_t_fhx3 = [599, 524];
%T0_t_fhx4 = [(T0_p_fhx4(2) + T0_s_fhx4(4))/2 , (T0_p_fhx4(4) + T0_s_fhx4(2))/2];
T0_t_fhx4 = [599, 524];

%Coolant Salt Initial Fertile Heat Exchanger Lump Temperatures
%To make sure that the initial temperature leaving Fertile HX is 1125 deg F
%initialize the secondary temperatures with a temp difference of 7 deg F
%between each lump.
%T0_s_bhx1 = [T0_sbhx_in+3.889, T0_sbhx_in+2*3.889];
T0_s_bhx1 = [615, 617];
%T0_s_bhx2 = T0_s_bhx1;
T0_s_bhx2 = [615, 617];
%T0_s_bhx3 = T0_s_bhx1;
T0_s_bhx3 = [615, 617];
%T0_s_bhx4 = T0_s_bhx1;
T0_s_bhx4 = [615, 617];

%Fertile Heat Exchanger Tube Lumps Initial Temperatures
%Initial Temps for secondary and fertile salt are for full power and for heat
%exchange to occur from fertile to the secondary salt. Therefore, initialize
%tube temperatures to be a temperature between the lowest temperature of
%the fertile lumps it's associated with and the temperature of the highest
%secondary lump it's associated with (e.g. Tube Lump tb initial temp should
%be between the blanket lump 2 initial temp of 1150 deg F and the secondary
%lump 2 initial temp of 1125 deg F)
%T0_t_bhx1 = (T0_p_bhx1(2) + T0_s_bhx1(2))/2;
T0_t_bhx1 = 625;
%T0_t_bhx2 = (T0_p_bhx2(2) + T0_s_bhx2(2))/2;
T0_t_bhx2 = 625;
%T0_t_bhx3 = (T0_p_bhx3(2) + T0_s_bhx3(2))/2;
T0_t_bhx3 = 625;
%T0_t_bhx4 = (T0_p_bhx4(2) + T0_s_bhx4(2))/2;
T0_t_bhx4 = 625;

%Set Initial Temperature of Pre-core Fuel Mixing Lump to the weighted
%average of the product of fuel flow rates and temperature.
T0f_in = (W_f1*T0_p_fhx1(4)+W_f2*T0_p_fhx2(4)+W_f3*T0_p_fhx3(4)+W_f4*T0_p_fhx4(4))/(W_f1+W_f2+W_f3+W_f4);

%Set Initial Temperature of Pre-core Blanket Mixing Lump to the weighted
%average of the product of fuel flow rates and temperature.
T0b_in = (W_b1*T0_p_bhx1(2)+W_b2*T0_p_bhx2(2)+W_b3*T0_p_bhx3(2)+W_b4*T0_p_bhx4(2))/(W_b1+W_b2+W_b3+W_b4);

%Set Initial Temperature of Pre-BoP Coolant Mixing Lump to the weighted
%average of the product of coolant flow rates and temperature.
T0s_in = (W_s1*T0_s_bhx1(2)+W_s2*T0_s_bhx2(2)+W_s3*T0_s_bhx3(2)+W_s4*T0_s_bhx4(2))/(W_s1+W_s2+W_s3+W_s4);

%% OTSG Parameters
%OTSG parameters have been changed from the test as the core was being overcooled in
%steady-state conditions and caused the power to be produced above 100%.
%Since the OTSG in this plant is not designed specifically, we can change
%parameters how we want for our first plant run.

T_fw  = 212; % [�C] Feedwater Input Temperature
P_set = 13.0; % [MPa] Steam Pressure Setpoint
P_low = P_set-2.0; %MPa Lower Pressure for establishing Pressure/Temp/Enthalpy variance
% heat capacities
% % Cp_p is specific heat capacity of the secondary coolant salt
% % Salt in Table 3.2 of ORNL-TM-3996.
Cp_p = Cp_sec; % [[MJ/(kg/�C)] primary fluid of the OTSG

% % Cp_w is the specific heat capacity of the OTSG Tubes
% % 0.112; [BTU/F-lbm] Hastelloy N (tube), reference matweb at 752 deg F
% % At initial temperatures, initial tube wall is around 400 deg C on
% average. Use specific heat at about 400 deg C/ 752 deg F
Cp_w = 468.992e-6; % [[MJ/(kg/�C)]

Cp_fw = XSteam('Cp_pT',P_set*10,T_fw)/1000; % [[MJ/(kg/�C)] Specific heat of feedwater, 
% at inlet temp and Pressure setpoint
T_satcalc = XSteam('Tsat_p',P_set*10); %deg C, Saturation temp at pressure
T_sc2 = (T_satcalc + T_fw)/2; % temperature of subcooled lump set to halfway 
% between feedwater temp and saturation temp
Cp_sc = XSteam('Cp_pT',P_set*10,T_sc2)/1000; % [[MJ/(kg/�C)] Specific heat 
% of subcooled lump 

% % rho_p is the density of coolant salt at approx. 1000 deg F.
rho_p = den_s; % [kg/m^3], density of the primary fluid (coolant salt)                  
rho_w = den_t; % 8774.528 kg/m^3, 0.317 lb/in^3 density of Hastelloy N Tube                    
rho_fw = XSteam('rho_pT',P_set*10,T_fw); % [kg/m^3] density of feedwater
rho_f = XSteam('rhoL_p',P_set*10); % [kg/m^3] boiling water density 
% (sat liquid water at 13.0 MPa)          
rho_sc = (rho_fw + rho_f)/2; % [kg/m^3] subcooled fluid density        
rho_b = XSteam('rhoV_p',P_set*10); % [kg/m^3] boiling lump fluid density 
% (vapor at 13.0 MPa)               

P_s = P_set-0.5; % [MPa] (Initial Steam Pressure)
T_s = (540 + 600)/2; % deg C (Steam Temperature)
Cp_s  = XSteam('Cp_pT',P_s*10,T_s)/1000; % [MJ/kg-C] specific heat of superheated steam
rho_s = XSteam('rho_pT',P_s*10,T_s); % [kg/m^3] density of superheated steam

%%%%%%%%%%%%%%%%%%%%
%%% SG TUBE SIDE %%%
%%%%%%%%%%%%%%%%%%%%

%N = 6546; % Number of Tubes, V_Singh Thesis pg. 123
N = 6400;
%L = 20.0; % [m] Active Tube Length in m,
L = 22.0;
% % Initialize the subcooled/boiling/superheated tube lengths as percentages of
% % full active tube length.
% % Note that the initial superheated length is about 55% of total length,
% % the boiling length is about 28%, and the subcooled length is about
% % and subcooled length is about 17% of the total length in V_Singh MSDR Power Plant model.
%L_b = L*0.28; % [m] Boiling Length
L_b = 1.9;% [m] Boiling Length
%L_s = L*0.55; % [m] Superheated Steam Length
L_s = 16.53; % [m] Superheated Steam Length
L_sc = L-L_b-L_s; % [m] Subcooled Length 
% 5/8 tubes with 0.034" thickness. V_Singh MSDR Power Plant Model and Chen Thesis pg. 23 
D_ot = 0.015875; % [m] 0.0520833 [ft] outer tube diam.
T_th = 0.000864; % [m] 0.002833333 [ft]
D_it = D_ot-2*T_th; % [m] internal tube diameter
R_it = D_it/2; % [m] inner tube radius
R_ot = D_ot/2; % [m] outer tube radius
A_sit = pi()*D_it*L*N; % [m^2] total surface area inner tube 
A_sot = pi()*D_ot*L*N; % [m^2] total surface area outer tube 

%%%%%%%%%% Area %%%%%%%%%%
D_OTSG_i = 4.00; %m, Inner Diameter of OTSG Shell.
A_p_OTSG = N*pi()*(D_it)^2/4; %m^2, total primary fluid cross sectional area in the OTSG
A_tube = (pi()*((D_ot)^2)/4) - (pi()*((D_it)^2)/4); %m^2, cross sectional area of one tube in OTSG
A_w_OTSG = A_tube*N; %m^2, total cross sectional area of tubes in OTSG
A_s = pi()*((D_OTSG_i^2)/4) - A_w_OTSG - A_p_OTSG; % [m^2] Cross sectional secondary flow (say inner diameter of the shell of OTSG is 70.6 in or 1.79324 m wide)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INITIAL STATE CALCULATIONS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_stm = 0.018; % [kg/mol] Molar weight of steam

P_table = P_low:0.1:P_set; % Pressure for linearizing Psat vs. Tsat relationship;

Ts_avg = (540 + 600)/2; % C average steam temperature, 
T_table = [];
Hfg_table = [];
hs_table = [];
rho_bV_table = [];
T_fsat_table = []; %note that Tfsat is calculated with a 1.0 MPa to Set Pressure MPa
% this is done to match previous MSDR work K_1 values
rhosc_table = []; %subcooled density table, also done to match previous
% MSDR K_sc values

for PPP = P_low:0.1:P_set

T_sat = XSteam('Tsat_p',PPP*10); %deg C, saturated temp at pressure
hf = XSteam('hL_p',PPP*10)/1000; %Mj/kg, sat fluid enthalpy
hg = XSteam('hV_p',PPP*10)/1000; %Mj/kg, sat vapor enthalpy
kf = XSteam('tcL_p',PPP*10)*1000; %W/m deg C, sat fluid conductivity
kg = XSteam('tcV_p',PPP*10)*1000; %W/m deg C, sat vapor conductivity

hss = XSteam('h_pT',PPP*10,Ts_avg)/1000; %Mj/kg, superheat steam enthalpy
rho_bV = XSteam('rho_pT',PPP*10,Ts_avg); %kg/m^3; density of fluid over boiling region
T_table = [T_table, T_sat];
hs_table = [hs_table, hss];
Hfg_table = [Hfg_table, hg-hf];
rho_bV_table = [rho_bV_table, rho_bV];

end

PQ_table = 1.0:0.1:P_set;

for PPQ = 1.0:0.1:P_set

T_fsat = XSteam('Tsat_p',PPQ*10);
T_fsat_table = [T_fsat_table, T_fsat];
rhosc = XSteam('rhoL_p',PPQ*10);
rhosc_table = [rhosc_table, rhosc];

end

a = polyfit(P_table, T_table, 1);
X_5 = a(2); K_5 = a(1);
b = polyfit(P_table, Hfg_table, 1);
X_4 = b(2); K_4 = b(1);
c = polyfit(P_table, hs_table, 1);
dHs_dPs = c(1);
d = polyfit(P_table, rho_bV_table, 1);
K_b = d(1);
e = polyfit(PQ_table, T_fsat_table, 1);
K_1 = e(1);
f = polyfit(PQ_table, rhosc_table, 1);
K_sc = f(1);

% % Initial temperatures for the OTSG nodes
%%%%%%%%%% Temperature %%%%%%%%%%
% For primary coolant, calculate the overall delta T of the OTSG
% for the primary coolant in the MSBR. Then calculate the overall delta
% T calculated for the MSDR. Figure out the individual delta Ts between
% subsequent nodes (e.g. T_p1-T_p2) in the MSDR. Calculate a percent
% difference of the overall delta T for the MSDR on each node. Use the
% percent difference to determine the appropriate initial temperatures for
% the MSBR.
%T_p1 = 590.24; % [�C] Primary Coolant Temperature node 1 1094 deg F
T_p1 = 608; %deg C Primary Coolant Temperature Lump 1 
%T_p2 = 560.40; % [�C] Primary Coolant Temperature node 2 1041 deg F
T_p2 = 576; %deg C Primary Coolant Temperature Lump 2 
%T_p3 = 521.70; % [�C] Primary Coolant Temperature node 3 971 deg F
T_p3 = 538; %deg C Primary Coolant Temperature Lump 3 
%T_p4 = 489.78; % [�C] Primary Coolant Temperature node 4 914 deg F
T_p4 = 508; %deg C Primary Coolant Temperature Lump 4 
%T_p5 = 474.66; % [�C] Primary Coolant Temperature node 5 886 deg F
T_p5 = 491; %deg C Primary Coolant Temperature Lump 5 
%T_p6 = 454.44; % [�C] Primary Coolant Temperature node 6 850 deg F
T_p6 = 470; %deg C Primary Coolant Temperature Lump 6 
% 
%T_w1 = 533.07; % [�C] Temperature for wall node 1 992 deg F
T_w1 = 606; % [�C] Temperature for wall node 1 
%T_w2 = 489.88; % [�C] Temperature for wall node 2 914 deg F
T_w2 = 568; % [�C] Temperature for wall node 2 
%T_w3 = 383.52; % [�C] Temperature for wall node 3 722 deg F
T_w3 = 463; % [�C] Temperature for wall node 3 
%T_w4 = 374.75; % [�C] Temperature for wall node 4 707 deg F
T_w4 = 444; % [�C] Temperature for wall node 4 
%T_w5 = 385.17; % [�C] Temperature for wall node 5 725 deg F
T_w5 = 449; % [�C] Temperature for wall node 5 
%T_w6 = 345.83; % [�C] Temperature for wall node 6 654 deg F
T_w6 = 419; % [�C] Temperature for wall node 6 
% use same feedwater and temperature initial conditions as V_Singh's
% Thesis.
%T_s1  = 482.63; % [�C] Temperature for superheated node 901 deg F
T_s1  = 597; % [�C] Temperature for superheated node
%T_s2  = 427.66; % [�C] Temperature for superheated node 802 deg F
T_s2  = 538; % [�C] Temperature for superheated node

%%%%%%%%%% Pressure %%%%%%%%%%
P_ss0 = P_set; % [MPa] Pressure superheated steam lump
P_s   = P_set; % [MPa] Pressure superheated steam lump
P_sc  = P_set; % [MPa] Pressure subcooled initial

deltaP = 0.5;
P_sat = (P_s + deltaP);
%X5=402.94; K5=0.14;  %Tsat~Psat
T_sat = X_5 + K_5*P_sat;
%Tsat=546.6;  %Exit temperature=317C and Degree of superheat is 43.4;
H_fg = X_4 + K_4*P_sat;


%%%%%%%%%% Mass Flow Rate %%%%%%%%%%
% Coolant salt flow rate at full power ORNL-3996 Fig. 1. Since there is
% no Reheater in the model, assume full flow from Coolant Salt Loop goes to
% the OTSG, 130 ft^3/s + 20 ft^3/s = 150 ft^3/s = 4.247528 m^3/s
W_p0 = W_s1+W_s2+W_s3+W_s4; %kg/s

W_fw = 896; % [kg/s]Mass flow rate of feedwater, see note at beginning of
% OTSG module regarding increased flow rate.

% Compressibility factor. See Chen's Thesis pg. 55 for discussion on how
% to obtain. Based on initial steady state conditions.
Z_ss = 0.76634; % Compressibility factor at 376.5 C, 12.5 MPa

R = 8.314462E-6; % [MJ/mol-�C] % Universal gas constant 

%%%%%%%%%% HEAT TRANSFER COEFFICIENTS %%%%%%%%%%
% heat exchanger tubes for the OTSG are assumed to be made out of Hastelloy
% N material as stated in ORNL-3996. Therefore, heat transfer coefficient
% of the coolant to the Hastealloy N tubes in the OTSG should be similar
% to the heat transfer coefficient of the coolant salt in the primary and
% secondary heat exchangers.
h_pw = h_s; %[MW/m^2-�C] 

%OTSG Original heat transfer coefficients. Used original heat transfer
%coefficients of the OTSG in MSDR 2019 Plant Paper by V_Singh et. al.
h_ws = 5732.378387768E-6; % [MW/m^2-�C] 0.6672; % htc wall to steam node [BTU/s-ft^2-F]
h_wb = 13349.334646671E-6; % [MW/m^2-�C] 2.16647; % htc wall to boil node [BTU/s-ft^2-F]
h_wsc = 8385.005556375E-6; % [MW/m^2-�C] 1.18; %0.147 htc wall to subcooled node [BTU/s-ft^2-F]

%% Steam Enthalpy
%Populate lookup table for BoP Steam Enthalpy input. Function uses 
%XSteam to take an array of possible steam temperatures and an array of
%possible steam pressures. XSteam provides the enthalpy of steam at the
%given pressure and temperature. This lookup table will be used to find 
%Steam Enthalpy for the steam going into the nozzle chest.
% temp_table (deg C)
% pres_table (MPa)
% H_s_table (MJ/kg)
run('Xsteam/SteamEnthalpy.m');

%% Main Steam Valve
%Valve Flow Coefficient of Main Steam Valve, directs percent of steam to
%nozzle chest. (1-C_v) is the amount of flow that goes to the reheater
%C_v = 0.84; % ratio of total flow that goes to the nozzle chest
C_v = 0.84; % ratio of total flow that goes to the nozzle chest
%% Nozzle Chest

H_nc = XSteam('h_pT', P_s*10, T_s)/1000;%MJ/kg, enthalpy of steam from OTSG
%Density of steam from OTSG in nozzle chest
rho_nc = XSteam('rho_pT', P_s*10, T_s);%kg/m^3,
%Nozzle Chest Volume, Naghedolfeizi's Thesis pg. 132 for a 3436 MW(t) plant
V_nc = 5.6634; %m^3, 200 ft^3
%Populate lookup table for Nozzle chest pressure. Function uses XSteam.m to
%take an array of possible densities in the nozzle chest and an array of
%enthalpies in the nozzle chest. XSteam.m provides the pressure at the
%density and entalpy given. This lookup table will be used to find nozzle
%chest pressure from the calculated enthalpy and density when running sim.
% P_nc_table (MPa)
% H_table (MJ/kg)
% rho_table (kg/m^3)
run('Xsteam/NozzleChestPressure.m');
%Populate lookup table for Nozzle chest specific volume of saturated vapor. 
%Function uses XSteam.m to take an array of possible pressure in the nozzle 
%chest. XSteam.m provides the specific volume of saturated vapor at pressure
%given. This lookup table will be used to find nozzle chest specific volume
%when running sim.
% SV_nc_table (m^3/kg)
% P_rh_table1 (MPa)
run('Xsteam/NozzleChestSpecVol.m');
%Populate lookup table for Nozzle chest specific volume of saturated fluid. 
%Function uses XSteam.m to take an array of possible pressure in the nozzle 
%chest. XSteam.m provides the specific volume of saturated fluid at pressure
%given. This lookup table will be used to find nozzle chest specific volume
%when running sim.
% SVf_nc_table (m^3/kg)
% P_rh_table1 (MPa)
run('Xsteam/NozzleChestfSpecVol.m');
%Populate lookup table for Reheater enthalpy of saturated fluid. 
%Function uses XSteam.m to take an array of possible pressure in the 
%reheater. XSteam provides enthalpy of reheater sat fluid at given pressure
%This lookup table will be used to find reheater enthalpy
%when running sim.
% H_rh_table (MJ/kg)
% P_rh_table1 (MPa)
run('Xsteam/ReheaterEnthalpy_satliquid.m');
%Populate lookup table for Reheater enthalpy of vaporization. 
%Function uses XSteam.m to take an array of possible pressure in the 
%reheater. XSteam provides enthalpy of sat vapor and sat fluid 
%at given pressure. Difference between the enthalpies provides the enthalpy
%of vaporization. This lookup table will be used to find the enthalpy of 
%vaporization of fluid in reheater when running sim.
% H_rhvap_table (MJ/kg)
% P_rh_table1 (MPa)
run('Xsteam/ReheaterEnthalpy_vaporization.m');

%Initial Specific volume of sat fluid and vapor in the nc
nu_f0 = XSteam('vL_p', P_s*1e1);% [m^3/kg], input pressure converted from MPa to bar
nu_g0 = XSteam('vV_p', P_s*1e1);% [m^3/kg], input pressure converted from MPa to bar
%Initial Fluid enthalpy and enthalpy of vaporization for the nc
H_frh0 = XSteam('hL_p', P_s*1e1)/1e3;% [MJ/kg] input pressure converted from MPa to bar
H_fgrh0 = (XSteam('hV_p', P_s*1e1) - XSteam('hL_p', P_s*1e1))/1e3;% [MJ/kg] input pressure convereted form MPa to bar

C_nc_X_A_nc = 12.8; %product of C_nc and A_nc (flow coefficient and area of
% nozzle chest outlet). 12.8 is the original value used in the MSDR Plant.
% Naghedolfeizi's Thesis pg. 159 uses a COF of 3.0326 (ft^3/s/psi) which is
% about 12.45 (m^3/s/MPa) with an equation that does not square root the 
% pressure and density products.

%% High Pressure Turbine

% HP turbine steam enthalpy at isentropic endpoints
% Terms of the equation (36) on the Shankar Simulation paper page 272
% converted from terms of Btu, lbm, W*s, and psi to MJ, kg, and MPa
% Also see Naghedolfeizi's Thesis pg. 151
term1 = 1080.3*1054.35e-6/0.45359237; % [MJ/kg]
p1 = 200*6894.76/1e6; % [MPa]
p2 = 1000*6894.76/1e6; % [MPa]
ce1 = 1.2471e-07; % [m^3/kg]
ce2 = 5.3816e-14; % [m^6/kg^2]
ce3 = 3.3721e-08; % [m^3/kg]
%isentropic efficiency of a High Pressure Turbine. Value from Shankar pg.
%273
eta_hp_isen = 0.861;
% fraction of steam that bleeds out of the HPT.
% Value from Naghedolfeizi's Thesis pg. 159
K_bhp = 0.1634;
%HPT Residence time, Value from Naghedolfeizi's Thesis pg. 159
tau_hpt = 2; %s
%Initial Mass Flow Rate out of the High Pressure Turbine
%Took number provided by the MSDR Plant paper and increased it by
%multiplying 2.8
w_hp_out = 225.96; %[kg/s]

w_bhp_out = 44.134; % [kg/s], flow from HP Turbine to the High Pres FWH
% number is the initial flow rate used in the MSDR plant paper increased by
% multiplying 2.8.

%Enthalpy of isentropic endpoints of HPT and enthalpy of the fluid leaving
%HPT. Used original values from the MSDR Plant Paper
H_hpt_out = 2.6201; % [MJ/kg]
H_hp_isen0 = 2.5111; % [MJ/kg]

%% Moisture Seperator

%No values to provide as the moisture seperator simply take flow from the
%HPT, removes the liquid water and sends it to the HP feedwater heater, and
%sends the dry steam to the reheater.

%% Reheater
% Volume of reheater, Naghedolfeizi's Thesis pg. 160 uses 20000 ft^3 
% MSDR Plant used 1/10th of Naghedolfeizi volume and doubled it.
V_rh = 56.634*2;% [m^3] (20000 ft^3)
%Reheater residence time, 3 s is what is used in Naghedolfeizi's Thesis pg. 159
%MSDR Plant Paper Model appears to have doubled the residence time based
%on the increased volume.
tau_rht = 3*2; %s (used Naghedolfeizi residence time of 3 s)
%time required for reheater heat transfer, 10 s is what is used in Naghedolfeizi's Thesis pg. 159
%MSDR Plant Paper Model appears to have doubled the residence time based
%on the increased volume.
tau_rhs = 10*2; %s (used Naghedolfeizi residence time of 10 s)
C_rh_X_A_rh =33; % product of C_rh and A_rh (flow coefficient and area of reheater nozzle).
% 33 is the value used in the MSDR Model.
% Using values from Naghedolfeizi's Thesis this product is 16.67 
% (m^(3/2)kg^(1/2)/MPa^(1/2)) for the plant modeled there. Original value 
% in Naghedolfeizi is 431.35 (ft^(3/2)lb^(1/2)/psi^(1/2)). Consider using that value.

%Isobaric specific heat capacity of reheater fluid, Steam @ 277 deg C
%and 1.26 MPa. This value is from the MSDR Plant. Naghedolfeizi uses a
%specific heat of 22.589 btu/lbm*deg F which is for a much higher
%temperature (700 deg f) and pressure (3090 psi) for water vapor.
Cp_rh = 0.0021520; %[MJ/kg-�C]

%Populate lookup table for Reheater Pressure. 
%Function uses XSteam to take an array of possible densities and enthalpies
%of the fluid in the reheater. XSteam provides pressure given an enthalpy 
%and density. This lookup table will be used to find the pressure in the 
%reheater vaporization of fluid in reheater when running sim.
% rho_rh_table (kg/m^3)
% H_rh_table2 (MJ/kg)
% P_rh_table (MPa)
run('Xsteam/ReheaterPressure.m');
%Populate lookup table for Reheater enthalpy of saturated vapor. 
%Function uses XSteam.m to take an array of possible pressure in the 
%reheater. XSteam provides enthalpy of reheater sat vapor at given pressure
%This lookup table will be used to find reheater enthalpy
%when running sim.
% H_rh_table1 (MJ/kg)
% P_rh_table1 (MPa)
run('Xsteam/ReheaterEnthalpy_satvapor.m');

rho_rh = 13.254; % [kg/m^3], initial reheater density
%density initial value is from MSDR Plant.
%Initial enthalpy of the fluid in the Reheater
H_rh = 2.8115; % [MJ/kg] %enthalpy initial value is from MSDR Plant.
%
P_rh0 = 2.8582;% [MPa] initial pressure of the reheater fluid
%Heat transfer rate in the reheater, Original value from 2019 MSDR Plant
%model is 1.745 MJ/s. Used the number from Naghedolfeizi's Thesis (pg. 159) for the
%initial heat transfer rate of the reheater
Q_rh = 1.745; %[MJ/s]
%Output reheat steam flow rate, this value was defined but not used in the
%MSDR plant model.
w_rh_out = 33.7926; % [kg/s]
%Initial Output reheat steam flow rate, value below is the value from the 
%Simulink portion of MSDR Plant model multiplied by 2.8 since the plant we are
%modeling is larger. Consider using a different value
%(like the one above) if needed or using a value close to a value in
%Naghedolfeizi's Thesis
%w_rh_out = 203.11; %[kg/s]
%Initial enthalpy of vaporization of the steam leaving the reheater
%Value taken from the MSDR Plant Simulink definition
H_grh0 = 2.8029; %MJ/kg
%Initial flow rate from the reheater to the HP FWH. Value taken from MSDR
%Plant paper and increased by 2.8 due to plant size. Consider increasing 
%using a value from Naghedolfeizi's Thesis.
w_vrh_out = 51.447; %[Kg/s]
%Initial Temperature of the steam leaving the reheater. Value taken from
%the MSDR Plant Paper
T_rh0 = 466.87;%deg C

%% Low Pressure Turbine
%Isentropic efficiency of LP turbine. Naghedolfeizi's Thesis pg. 160
eta_lp_isen = 0.861; %
%Enthalpy of steam at the isentropic end point of the low pressure turbine
%Value taken from Naghedolfeizi's Thesis pg. 159 and converted to MJ/kg.
%This matches the original MSDR number.
H_lp_isen = 958.4*1054.35e-6/0.45359237; % [MJ/kg]
%Fraction of Bleed steam from the LP Turbine that goes to 
% the low pressure feedwater heater. Naghedolfeizi's Thesis pg. 159
K_blp = 0.2174; % fraction of bleed steam from LP turbine
%Residence time in the low pressure turbine, Naghedolfeizi's Thesis pg. 159
tau_lpt = 10; % time constant for LP turbine
%Initial Mass Flow Rate out of the Low Pressure Turbine
%Used number provided in 2019 MSDR Plant model and increased it based 
%off of the larger power and flow rate of our plant (MSDR_value*2.8).
w_lp_out = 158.95; % [kg/s]
%Initial bleed flow rate from the LP Turbine to the LP FWH. Used the MSDR
%Plant's value and increased by 2.8
w_blp_out = 44.156; % [kg/s]
%Initial enthalpy of the fluid leaving the LP Turbine. Directly from the
%MSDR Plant paper.
H_lpt_out = 2.3089; %[MJ/kg]

%% Low Pressure Feedwater Heater (Feedwater heater 1)
%Time constant for Low Pressure Feedwater Heater (heater 1)
%Naghedolfeizi's Thesis pg. 159 uses a time of 65 seconds
%100s is the value from MSDR Plant
tau_fwh1 = 100; %s
%constant heat transfer parameter for FWHl, Original 2019 MSDR Plant
%Parameter is here. Consider adjusting since Plant Thermal MWt, feedwater
%flow rate, and feedwater temperature is larger
% Naghedolfeizi uses a value of 1.10485 MJ/kg (475 Btu/lb)
h_fwh1 = 0.5927; % [MJ/kg]
% Initial Enthalpy of the Low Pressure Feedwater Heater (heater 1)
%MSDR Plant Value, Naghedolfeizi value 0.65221 MJ/kg = 280.4 Btu/lb
H_fwh10 = 0.46009; % [MJ/kg]

Q_fwh10 = 96.368;%[MW] Feedwater Heater 1 initial heat transfer rate. 
% MSDR Plant Value 96.368 MW. Value used here is 2.8 times the MSDR number
% as the mass flow rate for other values used to determine this heat
% transfer has been increased by 2.8. Naghedolfeizi value is 549.5441 MW =
% 520867.42 Btu/s

%% High Pressure Feedwater Heater (Feedwater Heater 2)
%Time constant for High Pressure Feedwater Heater (heater 2) heat transfer
%between the tubes and the hot flows in the shell.
%Naghedolfeizi's Thesis pg. 159 appears to use a time of 22 seconds
tau_fwh2 = 65; %s MSDR Value 65s
%constant heat transfer parameter for FWH2, Original 2019 MSDR Plant
%Parameter is here. Consider adjusting since Plant Thermal MWt, feedwater
%flow rate, and feedwater temperature is larger
%Naghedolfeizi uses a value of 2.009106 MJ/kg (863.76 Btu/lb)
%h_fwh2 = 0.5927; % [MJ/kg] MSDR Value is 0.5927
h_fwh2 = 1.230803; %% [MJ/kg] Naghedolfeizi value
% Initial Enthalpy of the High Pressure Feedwater Heater (heater 2)
%original MSDR sim value. Naghedolfeizi uses a value of 529.15 Btu/lb =
%1.230803 MJ/kg
H_fwh20 = 0.6784; %[MJ/kg] 387.0721 enthalpy of heater 2 [BTU/lbm]
T_fw = 212; %deg C, Feedwater temperature

tau_rh2 = 10; %s time constant for heater 2 flow, used MSDR Plant Value
%Naghedolfeizi's Thesis pg. 159 appears to use a time of 10 seconds

%Original value for the w_fwh in the MSDR Plant paper is 97.0688 kg/s.
%It is suspected to be the flow from the High Pressure Feedwater heater to
%the Low pressure feedwater heater as it matches closely with the initial
%value specified in the model. MSDR Simulink used a value of 118.44 kg/s
%This value has equals the MSDR Simulink model value increase by 2.8.
% Naghedolfeizi uses a value of 1217.8 lb/s = 552.3848 kg/s
w_fwh = 118.44; %[kg/s] mass flow rate from heater 2

Q_fwh20 = 70.197;%[MW] Feedwater Heater 2 initial heat transfer rate. 
% MSDR Plant Value 70.197 MW. The value used here has taken the MSDR value
% and increased it by 2.8 as the flow rates used to determine this value
% have been increased by 2.8. Naghedolfeizi uses a value of 1084.25 MW

%% Condenser
%Populate lookup table for Condenser Saturated fluid enthalpy. 
%Function uses XSteam.m to take an array of possible pressure in the 
%condenser. XSteam provides enthalpy of sat fluid at given pressures. 
%This lookup table will be used to find the sat. fluid enthalpy of 
%the condenser when running sim.
% H_con_table (MJ/kg)
% P_con_table (MPa)
run('Xsteam/CondenserEnthalpy_satliquid.m');
%Populate lookup table for Condenser enthalpy of vaporization. 
%Function uses XSteam.m to take an array of possible pressure in the 
%condenser. XSteam provides enthalpy of sat fluid and sat vapor
%at given pressures and the function uses this to find enthalpy difference. 
%This lookup table will be used to find the enthalpy of vaporization of 
%the condenser when running sim.
% H_convap_table (MJ/kg)
% P_con_table (MPa)
run('Xsteam/CondenserEnthalpy_vaporization.m');

% Condensate and condenser hot well parameters

M_hw0 = 10000; %[kg] initial mass of the liquid in the condenser hotwell
%MSDR Plant value, consider increasing for larger plant size
%Naghedolfeizi has a constant mass of liquid water in the condenser at
%41422.9 lbm = 18789.11 kg
W_co0 = 159; %[kg/s] initial condensate flow rate from condenser to 
% feedwater heater. Initial value found by taken the MSDR Plant value of 159 kg/s and
% multiplying by 2.8. Naghedolfeizi's thesis has feedwater flow rate out of
% condenser at 470.0986 kg/s
V_co = 25; % Volume of condenser [m^3] Original value from the MSDR Simulink model
% Consider using 4000 ft^3 = 113.2674 (original value from the MSDR Matlab params
% file) or 56.78118 m^3 (condenser hotwell volume for the plant in
% Naghedolfeizi's Thesis pg. 132)
H_co0 = 0.149; % [MJ/kg] initial enthalpy of the condenser. MSDR Plant Value
H_f_co0 = 0.149; % [MJ/kg] initial enthalpy of sat fluid in condenser. MSDR Plant Value
% Naghedolfeizi uses a value of 0.162215 MJ/kg
H_fg_co0 = 2.35; % [MJ/kg] initial enthalpy of vaporization of fluid in condenser
%MSDR Plant Value
% Naghedolfeizi uses a value of 2.41 MJ/kg
%P_co0 = 1735.9; % condenser pressure [psf] %MSDR MATLAB Value is approx. 1
%psi
P_co0 = 0.006; %[MPa] initial condenser pressure value from Simulink of MSDR
% Naghedolfeizi uses approx. 1 psi
T_co0 = 90; %deg C initial condensate fluid temperature, MSDR Simulink value

alpha_co = 316.2; % deg C/MPa, coefficient for the linear relationship of
% condenser temperature and pressure, see Cao's paper, value here same as
% in the MSDR simulation.
beta_co = 68.0958; % deg C, constant for the linear relationship of the 
% condenser temperature and pressure, see Cao's paper, value here same as 
% in the MSDR simulation.
hA_co0 = 140; %kW/deg C initial condenser overall heat transfer coefficient
%original MSDR Plant value from Simulink model. Note the units are in kW.
Q_co0 = 375; %MW, initial heat duty of the condenser, original MSDR Plant Value is 375 MW
% Note that the 1190.1 value is equal to the MSBR Net thermal MW at
% steady-state (2225 MW) minus the gross power generation on Fig. 1 of
% ORNL-3996 (1034.9 MW(e)). The gross generation is the power extracted by
% the turbines and the remaining heat in the fluid needs to be dissipated
% through the condenser.

%Cooling Water Parameters
T_cwout0 = 40; %deg C initial temperature of the cooling water leaving the condenser
%Plant in MSDR simulink has 40 deg C. Plant in MSDR Matlab Params uses 69.78939 deg C. Note that the value in the MSDR
%MATLAB Parameter file is different.
%Tcwo = 157.6209; % temperature of cooling water out [F]
T_cwi = 21.11; % [C] temperature of inlet cooling water into condenser, 
% original MSDR Plant Value of 70 deg F
W_cwi = 40000;%[kg/s] cooling water mass flow rate
% original MSDR Code Wcw0 = 4000;%1397.7; % cooling water mass flow rate [lbm/s]
% MSDR Simulink uses a value of 2552 kg/s cooling water
C_pcw = 0.0042; %[MJ/kg deg C] Cooling water heat capacity. From Cao's
% Steam Condenser model paper, as implemented in the MSDR Plant model.

M_cw = 65000; % kg cooling water hold up, Original value from MSDR Simulink
wilson_ratio = M_cw/6500; %ratio between the mass of the cooling water and 
% that used in Cao's paper for the Wilson Plot Method relationship.
a_1 = 8.7292e-2/wilson_ratio;% a_1 coefficient for the Wilson Plot Method 
% as shown in Cao's paper. Coefficient represents the heat transfer 
% resistance of the inner fluid in the condenser tubes
a_2 = 7.3787e-4/wilson_ratio;% a_2 constant for the Wilson Plot Method 
% as shown in Cao's paper. Constant represents the heat transfer resistance
% of the outer fluid and tube wall of the condenser.

%% Fuel/Fertile/Coolant Pump Control
% Fuel Pump (Note that first index in the array controls the first pump).
% Subsequent indices control other pumps.
tau_p_f = [10, 10, 10, 10];
therm_conv_fuel = [0.05, 0.05, 0.05, 0.05];
trip_time_f = [11000, 11000, 11000, 11000];

% Blanket Pump (the number in the array corresponds to the pump/heat
% exchanger loop number).
tau_p_b = [10, 10, 10, 10];
therm_conv_blank = [0.05, 0.05, 0.05, 0.05];
trip_time_b = [11000, 11000, 11000, 11000];

%Coolant Pump (the number in the array corresponds to the pump/heat
%exchanger loop number).
tau_p_c = [10, 10, 10, 10];
therm_conv_cool = [0.05, 0.05, 0.05, 0.05];
trip_time_c = [11000, 11000, 11000, 11000];

%% Reactivity Control
%Ramp reactivity insertion
rampstart = 750; %time to begin ramped reactivity insertion
rampend = 760; %time to end ramped reactivity insertion
t = (rampstart:0.1:rampend); 
endramp = 2e-3; %value of reactivity inserted at end of ramp
slope = endramp/(rampend-rampstart);
ramp = ((slope*t-rampstart*slope).*(t>=rampstart)-(slope*t-rampstart*slope).*(t>=rampend)) + (rampend-rampstart)*slope*(t>=rampend);

downstart = 1000; %time to begin ramped reactivity withdrawal
downend = 1010; %time to end ramped reactivity withdrawal
td = (downstart:0.1:downend);
dendramp = 2e-3; %value of reactivity withdrawn at end of ramp
sloped = dendramp/(downend-downstart);
downramp = ((dendramp-sloped*td+downstart*sloped).*(td>=downstart)+(sloped*td-downstart*sloped).*(td>=downend)) - (downend-downstart)*sloped*(td>=downend);

reacttime = [0, rampstart:0.1:rampend, rampend , downstart:0.1:downend];
reactdata = [0, ramp                  , endramp , downramp];

%Step reactivity insertion
%reactdata is 20 pcm oscillation.
%reactdata = [0 0      0     0     0     ];
%reacttime = [0 1500   1800  2700  3600  ];

react = timeseries(reactdata,reacttime);

sim_time = 2000; %time to conduct simulation
ts_max = 1e-1; %max step size
ReltTol = 1e-4;

outputTemps = sim('MSBR_Plant.slx');

%Get a subset of the time data at the second coolant pump trip for each
%Fuel Salt, Fertile Salt, and Coolant Salt Lump Node

