%Purpose of this program is to plot the Coolant Salt Temperatures over time
% for every Coolant Salt node that is in the simulation. Run this program
% post-simulation to see the temperatures and compare values in different
% components in the system. The x-axis is limited based off of time of the
% transient accident. Change the x-axis limit or delete to see entire plant
% simulation.

% Secondary Coolant Temperature Plot Script

%% Plot Coolant Temps from Fuel heat exchangers

%Setup the entire Figure placement tiles
tiledlayout(3,3)
%Goto first tile in the Coolant Temperature Figure
nexttile
%Fuel HX 1 Coolant Temperatures of each lump of the HX over the transient
plot(Fuel_Heat_Ex_1_Cool_Temps)
xlim([1000 3500])
%plot Coolant Liquidus Temp in deg C
yline(371,'-.b','Coolant Liquidus')
%plot Fertile High Temp Limit in deg C
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fuel HX 1 Temps')
legend('Inlet','Lump 1','Lump 2','Lump 3','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

%Goto next tile and plot other items. Other items to plot are the Coolant
%Temperatures for the Coolant lumps in the Fuel Heat Exchangers over the transient
nexttile
plot(Fuel_Heat_Ex_2_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fuel HX 2 Temps')
legend('Inlet','Lump 1','Lump 2','Lump 3','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fuel_Heat_Ex_3_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fuel HX 3 Temps')
legend('Inlet','Lump 1','Lump 2','Lump 3','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fuel_Heat_Ex_4_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fuel HX 4 Temps')
legend('Inlet','Lump 1','Lump 2','Lump 3','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

%% Plot Coolant Temps from Fertile heat exchangers

nexttile
plot(Fert_Heat_Ex_1_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fert HX 1 Temps')
legend('Inlet','Lump 1','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_2_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fert HX 2 Temps')
legend('Inlet','Lump 1','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_3_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fert HX 3 Temps')
legend('Inlet','Lump 1','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_4_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt Fert HX 4 Temps')
legend('Inlet','Lump 1','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

%% Plot Coolant Temps from Oncce-Through Steam Generator

nexttile
plot(OTSG_Primary_Cool_Temps)
xlim([1000 3500])
yline(371,'-.b','Coolant Liquidus')
yline(780,'-.r','High Temp Control')
title('Coolant Salt OTSG Temps')
legend('Inlet','Lump 1','Lump 2','Lump 3','Lump 4','Lump 5','Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')