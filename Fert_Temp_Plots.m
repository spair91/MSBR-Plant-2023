%Purpose of this program is to plot the Fertile Salt Temperatures over time
% for every Fertile Salt node that is in the simulation. Run this program
% post-simulation to see the temperatures and compare values in different
% components in the system. The x-axis is limited based off of time of the
% transient accident. Change the x-axis limit or delete to see entire plant
% simulation.

%% Fertile Temperature Plot Script
%Setup the entire Figure placement tiles
tiledlayout(2,3)
%Goto first tile in the Fertile Temperature Figure
nexttile([1 2])
%Plot Core Fertile Temperatures of each lump of the core over the transient
plot(core_fert_temp)
xlim([600 1200])
%plot Fertile Liquidus Temp in deg C
yline(560,'-.b','Fertile Liquidus','FontSize',14)
%plot Fertile High Temp Limit in deg C
yline(780,'-.r','High Temp Control','FontSize',14)
title('Fertile Core Temps')
legend('Inlet Temp','Lump 1','Lump 2 (Core Exit)','FontSize',14)
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

%Goto next tile and plot other items. Other items to plot are the Fertile
%Temperatures for the Fertile lumps in the Heat Exchangers over the transient
nexttile
plot(Fert_Heat_Ex_1_Fert_Temps)
xlim([600 1200])
yline(560,'-.b','Fertile Liquidus')
yline(780,'-.r','High Temp Control')
title('Fert HX 1 Temps')
legend('Inlet Temp','Lump 1','Lump 2 (HX Exit)')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_2_Fert_Temps)
xlim([600 1200])
yline(560,'-.b','Fertile Liquidus')
yline(780,'-.r','High Temp Control')
title('Fert HX 2 Temps')
legend('Inlet Temp','Lump 1','Lump 2 (HX Exit)')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_3_Fert_Temps)
xlim([600 1200])
yline(560,'-.b','Fertile Liquidus')
yline(780,'-.r','High Temp Control')
title('Fert HX 3 Temps')
legend('Inlet Temp','Lump 1','Lump 2 (HX Exit)')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')

nexttile
plot(Fert_Heat_Ex_4_Fert_Temps)
xlim([600 1200])
yline(560,'-.b','Fertile Liquidus')
yline(780,'-.r','High Temp Control')
title('Fert HX 4 Temps')
legend('Inlet Temp','Lump 1','Lump 2 (HX Exit)')
ylabel(['Temperature (' char(176) 'C)'])
xlabel('Time(s)')