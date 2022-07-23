%Purpose of this program is to plot the Fuel Salt Temperatures over time
% for every Fuel Salt node that is in the simulation. Run this program
% post-simulation to see the temperatures and compare values in different
% components in the system. The x-axis is limited based off of time of the
% transient accident. Change the x-axis limit or delete to see entire plant
% simulation.

%% Fuel Temperature Plot Script
%Setup the entire Figure placement tiles
tiledlayout(2,3)
%Goto first tile in the Fuel Temperature Figure
nexttile([1 2])
%Plot Core Fuel Temperatures of each lump of the core over the transient
plot(core_fuel_temp)
xlim([600 1200])
%plot Fuel Liquidus Temp in deg C
yline(450,'-.b','Fuel Liquidus')
%plot Fuel High Temp Limit in deg C
yline(780,'-.r','High Temp Control')
title('Fuel Core Temps')
legend('Fuel Inlet Temp','Lump 1','Lump 2','Lump 3','Fuel Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

%Goto next tile and plot other items. Other items to plot are the Fuel
%Temperatures for the fuel lumps in the Heat Exchangers over the transient
nexttile
plot(Fuel_Heat_Ex_1_Fuel_Temps)
xlim([600 1200])
yline(450,'-.b','Fuel Liquidus','FontSize',14)
yline(780,'-.r','High Temp Control','FontSize',14)
title('Fuel HX 1 Temps')
legend('HX Inlet','Lump 1','Lump 2','Lump 3','Fuel Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fuel_Heat_Ex_2_Fuel_Temps)
xlim([600 1200])
yline(450,'-.b','Fuel Liquidus')
yline(780,'-.r','High Temp Control')
title('Fuel HX 2 Temps')
legend('HX Inlet','Lump 1','Lump 2','Lump 3','Fuel Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fuel_Heat_Ex_3_Fuel_Temps)
xlim([600 1200])
yline(450,'-.b','Fuel Liquidus')
yline(780,'-.r','High Temp Control')
title('Fuel HX 3 Temps')
legend('HX Inlet','Lump 1','Lump 2','Lump 3','Fuel Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fuel_Heat_Ex_4_Fuel_Temps)
xlim([600 1200])
yline(450,'-.b','Fuel Liquidus')
yline(780,'-.r','High Temp Control')
title('Fuel HX 4 Temps')
legend('HX Inlet','Lump 1','Lump 2','Lump 3','Fuel Exit')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')