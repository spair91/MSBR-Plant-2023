%3D Plots for Parametric Analysis (average values)
%Note, this script requires the following variables to be defined
%From the Fuel Pump Accident Parametric Analysis (trip of Fuel Pump 1):
%Min_Fuel_Salt a matrix of fuel HX/fertile HX tubes and fuel temps
% - This is the min/average fuel salt temperature in the Fuel HX 1
%Min_Fert_Salt a matrix of fuel HX/fertile HX tubes and fert temps
% - This is the min/average fertile salt temperature in the Fertile HX 1
%From the Reactivity Accident Parametric Analysis (large Reactivity increase in core):
%Max_Fuel_Core a matrix of fuel HX/fertile HX tubes and fuel core temps
% - This is the max/average fuel salt temperaturs in the core
%Max_Fuel_HX a matrix with fuel HX/fertile HX tubes and fuel HX temps
% - This is the max/average fuel salt temperatures in the fuel HX 1

% store number of fuel hx tubes to x
x_minfuelhx1 = Min_Fuel_Salt(:,1).';
%store number of fert hx tubes to y
y_minfuelhx1 = Min_Fuel_Salt(:,2).';
%store min/max/avg temps in z
z_minfuelhx1 = Min_Fuel_Salt(:,3).';
%define linear space of fuel hx tube numbers
x_v_minfuelhx1 = linspace(min(x_minfuelhx1),max(x_minfuelhx1),291*2);
%define linear space of fert hx tube numbers
y_v_minfuelhx1 = linspace(min(y_minfuelhx1),max(y_minfuelhx1),291*2);
%create x,y mesh of the linear spaces
[X_minfuelhx1,Y_minfuelhx1] = meshgrid(x_v_minfuelhx1,y_v_minfuelhx1);
%interpolate surface grid points to the scattered data in the x,y,z vectors
Z_minfuelhx1 = griddata(x_minfuelhx1,y_minfuelhx1,z_minfuelhx1,X_minfuelhx1,Y_minfuelhx1);
figure(1)
%plot surface plot
surf(X_minfuelhx1,Y_minfuelhx1,Z_minfuelhx1)
%take interpolated shading to remove points and grids
shading interp
%set vertices and faces of the patch for the fuel/fert liquidus or max
%limit
v_minfuel = [3000, 1800, 450;5000, 1800, 450;5000, 400, 450; 3000, 400, 450];
f_minfuel = [1,2,3,4];
hold on
%plot patch
patch('Faces',f_minfuel,'Vertices',v_minfuel,'FaceColor','blue','FaceAlpha',.2);
%labels
xlabel('Number of Fuel HX Tubes');
ylabel('Number of Fertile HX Tubes');
zlabel(['Fuel Salt Temperatures (' char(176) 'C)']);
title('Minimum Fuel Temperatures in Fuel Pump Accident');
%add labeling to show what happens in fuel salt when temperatures are too
%cold or hot.
text(4600,500,444,'\downarrow Fuel Precipitation','Color','b','FontSize',20);
text(4600,500,454,'\uparrow Liquid Fuel','Color','k','FontSize',20);
view(-78,22);

%repeat code as much as needed. be sure to change variable names and figure
%numbers 

x_minferthx1 = Min_Fert_Salt(:,1).';
y_minferthx1 = Min_Fert_Salt(:,2).';
z_minferthx1 = Min_Fert_Salt(:,3).';
x_v_minferthx1 = linspace(min(x_minferthx1),max(x_minferthx1),291*2);
y_v_minferthx1 = linspace(min(y_minferthx1),max(y_minferthx1),291*2);
[X_minferthx1,Y_minferthx1] = meshgrid(x_v_minferthx1,y_v_minferthx1);
Z_minferthx1 = griddata(x_minferthx1,y_minferthx1,z_minferthx1,X_minferthx1,Y_minferthx1);
figure(2)
surf(X_minferthx1,Y_minferthx1,Z_minferthx1)
shading interp
v_minfert = [3000, 1800, 560;5000, 1800, 560;5000, 400, 560; 3000, 400, 560];
f_minfert = [1,2,3,4];
hold on
patch('Faces',f_minfert,'Vertices',v_minfert,'FaceColor','blue','FaceAlpha',.2);
xlabel('Number of Fuel HX Tubes');
ylabel('Number of Fertile HX Tubes');
zlabel(['Fertile Salt Temperatures (' char(176) 'C)']);
title('Minimum Fertile Temperatures in Fuel Pump Accident');
text(4600,500,530,'\downarrow Fertile Solidification','Color','b','FontSize',20);
text(4600,500,579,'\uparrow Liquid Fertile Salt','Color','k','FontSize',20);
view(-78,22);

x_maxfuelcore = Max_Fuel_Core(:,1).';
y_maxfuelcore = Max_Fuel_Core(:,2).';
z_maxfuelcore = Max_Fuel_Core(:,3).';
x_v_maxfuelcore = linspace(min(x_maxfuelcore),max(x_maxfuelcore),291*2);
y_v_maxfuelcore = linspace(min(y_maxfuelcore),max(y_maxfuelcore),291*2);
[X_maxfuelcore,Y_maxfuelcore] = meshgrid(x_v_maxfuelcore,y_v_maxfuelcore);
Z_maxfuelcore = griddata(x_maxfuelcore,y_maxfuelcore,z_maxfuelcore,X_maxfuelcore,Y_maxfuelcore);
figure(3)
surf(X_maxfuelcore,Y_maxfuelcore,Z_maxfuelcore)
shading interp
v_maxfuel = [3000, 1800, 780;5000, 1800, 780;5000, 400, 780; 3000, 400, 780];
f_maxfuel = [1,2,3,4];
hold on
patch('Faces',f_maxfuel,'Vertices',v_maxfuel,'FaceColor','red','FaceAlpha',.2);
xlabel('Number of Fuel HX Tubes');
ylabel('Number of Fertile HX Tubes');
zlabel(['Fuel Salt Core Temperatures (' char(176) 'C)']);
title('Maximum Fuel Core Temperatures in Reactivity Accident');
text(4600,500,770,'\downarrow Below Fuel Limit','Color','k','FontSize',20);
text(4600,500,789,'\uparrow Exceeds Fuel Limit','Color','r','FontSize',20);
view(-78,22);